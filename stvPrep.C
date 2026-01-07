#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "THStack.h"
#include "TLegend.h"
#include "TF1.h"
#include "TF2.h"
#include "TLine.h"
#include "TMath.h"
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <array>


// get everything into a single Tree (TriggerTree)
// compute extra flags (tutorial_bin_config.txt will apply the selection cuts)

// run via:
// -> ./start_singularity.sh
// -> source setup_stv.sh
// -> root -l stvPrep.C


// ******************************
// helpful functions
//

// FV determination
bool FidVol(double x, double y, double z, bool reco){
        double radius   = 100.;  //cm
        double y_min    = -150.; //cm
        double y_max    = 150.;  //cm
        double z_center = 168.1; //cm
        double y_offset = 14.46; //cm

		// reco vertices are in [m], with no offsets needed
		if (reco) {
        	x *= 100.0;   // m to cm
        	y *= 100.0;
        	z *= 100.0;

			return (y > y_min && y < y_max && radius > std::sqrt(z*z + x*x));
    	}
		
		// truth vertices in [cm] and need offsets applied
        else {
			return (y+y_offset > y_min && y+y_offset < y_max && radius > std::sqrt((z - z_center)*(z - z_center) + x*x));
        }
}

// Flag for whether events originated outside the water tank (External)
bool External(double x, double y, double z){
        double radius   = 152.4;  //cm
        double y_min    = -198.; //cm
        double y_max    = 198.;  //cm
        double z_center = 168.1; //cm
        double y_offset = 14.46; //cm

        bool inside = (y+y_offset > y_min && y+y_offset < y_max && radius > std::sqrt((z - z_center)*(z - z_center) + x*x));

		return !inside;   // will return True if external
}

// bunch cutting (whether the time was within the bunches)
bool bunch_cutting(double cluster_time) {
    const double bunch_sigma = 3.10;    // determined based on MC bunch width (update to include BRF fitting)
    const double avg_diff = 0.0588;     // weighted avg between NCQE and CCinc bunch positions

	// CCinc centers extracted from bunch fitting (update if needed)
    static const double CC_centers[] = {
        30.252862036804867, 49.105181438904836, 67.80980088829278,
        86.76183849558346, 105.93922181966119, 124.93311941790408,
        143.7700212101684, 162.74865224853755, 181.70948297717655,
        200.68638425864887, 219.68223380220599, 238.56628356910397,
        257.5515436174118, 276.3551973029019, 295.2815361852104,
        314.0061564628827, 333.1094964976512, 351.92585490124424,
        371.1708384264726, 389.940088303305, 409.01620112483363,
        427.9521719492367, 446.7455147480343, 465.86427294746875,
        484.66211690272115, 503.5070533671551, 522.6163201579387,
        541.2982051398702, 560.2045692165119, 579.1811607183519,
        598.3018221115304, 617.0865686436197, 636.1698744814963,
        655.123547467177, 673.9199027401905, 692.7613636037937,
        711.8237422436717, 730.5638916091223, 749.919228564482,
        768.5699692095528, 787.5428263094714, 806.4940847845247,
        825.3183895637842, 844.5902725026449, 863.361165573104,
        882.1651463832068, 901.2162145453196, 920.1687853634885,
        938.8960929169189, 957.9092898626772, 976.8941158136465,
        995.9592593646615, 1015.0336713179253, 1033.7734506838333,
        1052.6480088471567, 1071.7670838139609, 1090.8310119667751,
        1109.6275205410964, 1128.5324628090907, 1147.400000303218,
        1166.274220214323, 1185.3048892154845, 1204.1996876180758,
        1223.0746752109535, 1241.8634871861623, 1261.0549285487182,
        1279.9332745578495, 1298.7671384166624, 1317.7638511749599,
        1336.463999839768, 1355.6829568044154, 1374.3838629549405,
        1393.5433775292236, 1412.1745975355464, 1431.3917931289625,
        1450.2161654845772, 1469.2769758317038, 1488.1973920412447,
        1507.1931834101015, 1526.112741833664, 1545.1331477778976
    };

    const int N = sizeof(CC_centers) / sizeof(double);

    for (int i = 0; i < N; i++) {
        double diff = fabs((CC_centers[i] - avg_diff) - cluster_time);
        if (diff < bunch_sigma) {
            return true;
        }
    }
    return false;
}

// charge vector cut
std::array<double,3> pairwise_relative_direction(
        const std::vector<double>& hitX,
        const std::vector<double>& hitY,
        const std::vector<double>& hitZ,
        const std::vector<double>& hitPE,
        const std::vector<double>& hitT)
{
    const size_t N = hitX.size();
    if (N < 2) {
        return {0.0, 0.0, 0.0};
    }

    // --- Combine hits ---
    struct Hit {
        double x, y, z, pe, t;
    };
    std::vector<Hit> hits;
    hits.reserve(N);

    for (size_t i = 0; i < N; i++) {
        hits.push_back({hitX[i], hitY[i], hitZ[i], hitPE[i], hitT[i]});
    }

    // --- Sort by time ---
    std::sort(hits.begin(), hits.end(),
              [](const Hit& a, const Hit& b){ return a.t < b.t; });

    // --- Get earliest hit time ---
    double t0 = hits.front().t;

    // --- Filter hits: keep only those within 13 ns of earliest ---
    std::vector<Hit> filtered;
    filtered.reserve(N);
    for (auto& h : hits) {
        if ((h.t - t0) <= 13.0) {
            filtered.push_back(h);
        }
    }

    const size_t M = filtered.size();
    if (M < 2) {
        return {0.0, 0.0, 0.0};
    }

    // --- Compute cumulative weighted pairwise direction ---
    double vx = 0.0, vy = 0.0, vz = 0.0;

    for (size_t i = 0; i < M; i++) {
        for (size_t j = i + 1; j < M; j++) {

            double dx = filtered[j].x - filtered[i].x;
            double dy = filtered[j].y - filtered[i].y;
            double dz = filtered[j].z - filtered[i].z;

            double r = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (r == 0) continue;

            double w = filtered[i].pe * filtered[j].pe;

            vx += w * (dx / r);
            vy += w * (dy / r);
            vz += w * (dz / r);
        }
    }

    // --- Normalize ---
    double norm = std::sqrt(vx*vx + vy*vy + vz*vz);
    if (norm == 0.0) {
        return {0.0, 0.0, 0.0};
    }

    return {vx / norm, vy / norm, vz / norm};
}

// find if the calculated charge vector is within the cuts
bool is_inside_rotated_parabola(double Z_val,
                                double X_val,
                                double a = -0.45,
                                double c = 0.4,
                                double b = 0.0,
                                double theta_deg = 275.0)
{
    const double theta = theta_deg * M_PI / 180.0;

    // inverse rotation
    double u =  Z_val * std::cos(theta) + X_val * std::sin(theta);
    double v = -Z_val * std::sin(theta) + X_val * std::cos(theta);

    double parabola_val = a*u*u + b*u + c;

    return (v < parabola_val);
}

// reconstructed energy
double reco_energy(double pe) {
    return (pe + 2.90) / 5.90;
}


// ******************************
// main function
//

void stvPrep(){

	// ********************************
	// defintions
	//

	// MC adjusted prompt window and spill end times
	const double bunch_time_cutoff = 1560.0;
	const double time_to_prompt_end = 258.0;   // based on data
	const double prompt_window_end = bunch_time_cutoff + time_to_prompt_end;


	// ********************************
	// load libraries (modify this path to any compiled toolanalysis directory)
	gSystem->Load("/exp/annie/app/users/doran/NC_CC_Nov_6_2025/lib/libDataModel.so");
	gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");


	// loop for executing over N files
	for (int i = 0; i < 1000; ++i) {       // change 1000 to N to match N input files
		std::string file =
			"/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/MC_" + 
			std::to_string(i) + ".root";   // modify path if needed

		std::cout << "********************************************** " << std::endl;

		if (gSystem->AccessPathName(file.c_str())) {
			std::cout << "Skipping missing file: " << file << std::endl;
			continue;
		}

		std::cout << "Running over file: " << file << std::endl;
        
		//Open file and trees
        TFile *f = new TFile(file.c_str(),"UPDATE");
		if (!f || f->IsZombie()) {
        	std::cerr << "ERROR: Could not open file" << std::endl;
        	continue;
    	}
        
		TTree *tTrig = (TTree*)f->Get("phaseIITriggerTree");
		TTree *tCluster = (TTree*)f->Get("phaseIITankClusterTree");
		if (!tTrig || !tCluster) {
        	std::cerr << "ERROR: Could not find one or both TTrees." << std::endl;
        	f->Close();
        	delete f;
        	continue;
    	}


		// *********************************
		// Current branches (what you need to calculate new quantities / branches)
		//

		int evNTrig, clusterEN, targetZ, nuPDG, clusterNumber, tankMRDCoinc, noveto, vetoHit, isNC, isQEL, isMEC;
		unsigned int clusterHits;
		double clusterPE, clusterTime, clusterCB, bunchTime, recoVtX, recoVtY, recoVtZ, nuvtxx, nuvtxy, nuvtxz;		
		std::vector<double>* hitT    = nullptr;
		std::vector<double>* hitPE   = nullptr;
		std::vector<double>* hitID   = nullptr;
		std::vector<double>* hitX    = nullptr;
		std::vector<double>* hitY    = nullptr;
		std::vector<double>* hitZ    = nullptr;
		std::vector<double>* MRDhitT = nullptr;

		// all the existing weight branches
		std::vector<double>* All0_weight = nullptr;
		std::vector<double>* All1_weight = nullptr;
		std::vector<double>* All2_weight = nullptr;
		std::vector<double>* All3_weight = nullptr;
		std::vector<double>* All4_weight = nullptr;
		std::vector<double>* All5_weight = nullptr;

		std::vector<double>* flux_horncurrent = nullptr;
		std::vector<double>* flux_expskin = nullptr;
		std::vector<double>* flux_piplus = nullptr;
		std::vector<double>* flux_piminus = nullptr;
		std::vector<double>* flux_kplus = nullptr;
		std::vector<double>* flux_kminus = nullptr;
		std::vector<double>* flux_kzero = nullptr;
		std::vector<double>* flux_pionine = nullptr;
		std::vector<double>* flux_pionqe = nullptr;
		std::vector<double>* flux_piontot = nullptr;
		std::vector<double>* flux_nucine = nullptr;
		std::vector<double>* flux_nucqe = nullptr;
		std::vector<double>* flux_nuctot = nullptr;



		// trig branches
        tTrig->SetBranchAddress("eventNumber",&evNTrig);
        tTrig->SetBranchAddress("TankMRDCoinc",&tankMRDCoinc);
        tTrig->SetBranchAddress("NoVeto",&noveto);
		tTrig->SetBranchAddress("vetoHit",&vetoHit);
		tTrig->SetBranchAddress("MRDhitT",&MRDhitT);
        tTrig->SetBranchAddress("trueNC",&isNC);
        tTrig->SetBranchAddress("trueQEL",&isQEL);
		tTrig->SetBranchAddress("trueMEC",&isMEC);
		tTrig->SetBranchAddress("trueTargetZ",&targetZ);        
		tTrig->SetBranchAddress("trueNuPDG",&nuPDG);
		tTrig->SetBranchAddress("trueNuIntxVtx_X",&nuvtxx);
        tTrig->SetBranchAddress("trueNuIntxVtx_Y",&nuvtxy);
        tTrig->SetBranchAddress("trueNuIntxVtx_Z",&nuvtxz);
		
		//Weights
    	tTrig->SetBranchAddress("weight_All0_UBGenie",&All0_weight);
    	tTrig->SetBranchAddress("weight_All1_UBGenie",&All1_weight);
    	tTrig->SetBranchAddress("weight_All2_UBGenie",&All2_weight);
    	tTrig->SetBranchAddress("weight_All3_UBGenie",&All3_weight);
    	tTrig->SetBranchAddress("weight_All4_UBGenie",&All4_weight);
    	tTrig->SetBranchAddress("weight_horncurrent_FluxUnisim",&flux_horncurrent);
    	tTrig->SetBranchAddress("weight_expskin_FluxUnisim",&flux_expskin);
    	tTrig->SetBranchAddress("weight_pioninexsec_FluxUnisim",&flux_pionine);
    	tTrig->SetBranchAddress("weight_pionqexsec_FluxUnisim",&flux_pionqe);
    	tTrig->SetBranchAddress("weight_piontotxsec_FluxUnisim",&flux_piontot);
    	tTrig->SetBranchAddress("weight_nucleoninexsec_FluxUnisim",&flux_nucine);
    	tTrig->SetBranchAddress("weight_nucleonqexsec_FluxUnisim",&flux_nucqe);
    	tTrig->SetBranchAddress("weight_nucleontotxsec_FluxUnisim",&flux_nuctot);
    	tTrig->SetBranchAddress("weight_piplus_PrimaryHadronSWCentralSplineVariation",&flux_piplus);
    	tTrig->SetBranchAddress("weight_piminus_PrimaryHadronSWCentralSplineVariation",&flux_piminus);
    	tTrig->SetBranchAddress("weight_kminus_PrimaryHadronNormalization",&flux_kminus);
    	tTrig->SetBranchAddress("weight_kzero_PrimaryHadronSanfordWang",&flux_kzero);
    	tTrig->SetBranchAddress("weight_kplus_PrimaryHadronFeynmanScaling",&flux_kplus);


		// cluster branches
		tCluster->SetBranchAddress("eventNumber",&clusterEN);
		tCluster->SetBranchAddress("clusterNumber",&clusterNumber);
		tCluster->SetBranchAddress("clusterPE",&clusterPE);
		tCluster->SetBranchAddress("clusterChargeBalance",&clusterCB);
		tCluster->SetBranchAddress("clusterHits",&clusterHits);
		tCluster->SetBranchAddress("clusterTime",&clusterTime);
		tCluster->SetBranchAddress("bunchTimes",&bunchTime);
		tCluster->SetBranchAddress("hitX",&hitX);
		tCluster->SetBranchAddress("hitY",&hitY);
		tCluster->SetBranchAddress("hitZ",&hitZ);
		tCluster->SetBranchAddress("hitT",&hitT);
		tCluster->SetBranchAddress("hitPE",&hitPE);
		tCluster->SetBranchAddress("hitDetID",&hitID);
		tCluster->SetBranchAddress("recoLeastSqVtxX",&recoVtX);
		tCluster->SetBranchAddress("recoLeastSqVtxY",&recoVtY);
		tCluster->SetBranchAddress("recoLeastSqVtxZ",&recoVtZ);


		// NCQE-like will include 2p2h (MEC) and be flagged PRIOR to pion absorption --> i.e. if a pion is produced but later absorbed by the nucleus (RES) it will be added to NCother
		// GENIE does not include an NC 2p2h model, so we will not include a + MEC flag on our true events


		// ********************************
		// Pre-scan cluster tree and build the selection map
		// 

		// we need to first scan the cluster tree to determine the total number of prompt clusters for every event
		// then we scan over the tree again, applying the cuts, and using the above result to veto multi-cluster events

		std::map<int, bool> passCluster;
		std::map<int, int> promptClusterCount;  // {eventNumber: count}

		Long64_t nClusterEntries = tCluster->GetEntries();

		// Pass 1 --> calculate the total number of prompt clusters per event
		std::cout << "PASS 1: Counting prompt clusters for " << nClusterEntries << " cluster entries..." << std::endl;
    	for (Long64_t i = 0; i < nClusterEntries; ++i) {
        	tCluster->GetEntry(i);
        
        	// Prompt window check: bunchTime < bunch_time_cutoff + time_to_prompt_end (prompt_window_end)
        	if (bunchTime < prompt_window_end) {
            	promptClusterCount[clusterEN]++;
        	}
    	}
		std::cout << "Finished Pass 1. Found counts for " << promptClusterCount.size() << " unique events." << std::endl;
		std::cout << "PASS 2: Applying selection cuts for clusters..." << std::endl;

		// Pass 2 --> Apply the preselection cuts
	
		// Loop over the Cluster Tree again to apply selection cuts
    	for (Long64_t i = 0; i < nClusterEntries; ++i) {
        	tCluster->GetEntry(i);
        
			// If a cluster passes, mark the event number as 'true' in the map.
			// If it doesn't pass, we do nothing. The event will either remain 
            // absent from the map (and thus default to false later) or already 
            // be true from a previous passing cluster.

			// Retrieve the calculated prompt cluster count for this event
        	int count = promptClusterCount[clusterEN];

			// -----------------------------
			// Cluster pre-selection (this will also include more advanced cluster selection cuts)

			// 1. Only one prompt cluster
        	if (count > 1) {
            	continue; // Skip the entire event's clusters if more than one prompt cluster was found
        	}

			// these next two are likely taken care of by cut #1

			// 2. Must be the primary cluster (earliest)
        	if (clusterNumber != 0) {
            	continue;
        	}

			// 3. Must be a prompt cluster
			if (bunchTime > prompt_window_end) {
            	continue;
        	}

			// 4. Min 10 hits
			if (clusterHits < 10) {
            	continue;
        	}

			// 5. Must be within the beam spill
			if (bunchTime > bunch_time_cutoff) {
            	continue;
        	}

			// 6. Bunch cutting (3.10ns sigma)
			if (!bunch_cutting(bunchTime)) {
            	continue;
        	}

			// 7. CCB cut (<= 0.6)
			if (clusterCB > 0.6) {
				continue;
			}

			// 8. Energy reco cut
			double currentRecoE = reco_energy(clusterPE);
        	if (currentRecoE < 5.0 || currentRecoE > 12.0) {
            	continue;
        	}

			
			// 9. Charge Vector Parabolic cut
			if (hitX && hitY) {
            	// cw[2] (Z) and cw[0] (X) are passed to the parabola function
            	std::array<double, 3> cw = pairwise_relative_direction(*hitX, *hitY, *hitZ, *hitPE, *hitT);

            	// is_inside_rotated_parabola(Z_val, X_val, ...)
            	if (is_inside_rotated_parabola(cw[2], cw[0]) == false) {
                	continue;
            	}
        	} else {
             	// Skip if hit vectors weren't properly loaded
             	std::cout << "Charge vector error: Hits weren't properly loaded!" << std::endl;
             	continue;
        	}

			// 10. Reconstructed vertex in FV
			// Pass recoVtxX, recoVtxY, recoVtxZ (in meters) and set 'reco' flag to true.
        	if (!FidVol(recoVtX, recoVtY, recoVtZ, true)) {
            	continue;
        	}

			
			// If all checks pass, mark the event number as 'true' in the map.
        	passCluster[clusterEN] = true;

        } 
			
		std::cout << "Found " << passCluster.size() << " events with passing clusters." << std::endl;
    	//}


		// *********************************
		// Created/added branches
		//

		bool NCQE_MC_Signal = false;   // true signal

		int NCQE_EventCategory = -1;
		// 0 = nuNCQE
		// 1 = External
		// 2 = Out-of-FV
		// 3 = CC
		// 4 = NCother (includes MEC, NC RES, NC DIS, NC on Hydrogen (which may be labeled as QEL), NC EL)
		// 5 = nubarNCQE

		bool NCQE_Selected = false;   // reco
		
		// total weights
		std::vector<double> All_weight;
		std::vector<double> flux_All;

		// branch names
		TBranch *NCQEsig    = tTrig->Branch("NCQE_MC_Signal",&NCQE_MC_Signal, "NCQE_MC_Signal/O");
		TBranch *NCQEevtcat = tTrig->Branch("NCQE_EventCategory",&NCQE_EventCategory, "NCQE_EventCategory/I");
		TBranch *NCQEreco   = tTrig->Branch("NCQE_Selected", &NCQE_Selected, "NCQEreco/O");
		TBranch *WAll       = tTrig->Branch("weight_All_UBGenie",&All_weight);
        TBranch *WfAll      = tTrig->Branch("weight_flux_all",&flux_All);


		// *********************************
		// Loop over trigger tree and fill branches
		//

		Long64_t nTrigEntries = tTrig->GetEntries();
    	std::cout << "Filling new branch for " << nTrigEntries << " trigger entries..." << std::endl;

    	// Loop over the Trigger Tree
    	for (Long64_t i = 0; i < nTrigEntries; ++i) {
        	tTrig->GetEntry(i); // Reads the eventNumber (evNTrig) and other existing data
        
        	// Reset the branch variable for the new entry
        	NCQE_Selected = false; 
        	NCQE_MC_Signal = false;
			NCQE_EventCategory = -1;

			bool hasMRDactivity = false;
			bool hasFMVactivity = false;
			bool extVol = false;
			bool trueFV = false;

			// MRD activity flag
            if ( (MRDhitT && !MRDhitT->empty()) || (tankMRDCoinc == 1) ) {
				hasMRDactivity = true;
        	}
			
			// Veto activity flag
			if (vetoHit == 1 or noveto == 0) {
            	hasFMVactivity = true;
        	}

			// ***************************************
			// Reco flag
			//
        	// Look up the event number in the map
        	auto it = passCluster.find(evNTrig);
        	if (it != passCluster.end() && it->second) {
            	// The event was found in the map and the value is true (at least one cluster passed) - now check the MRD / veto condition
				if (!hasMRDactivity && !hasFMVactivity) {
					NCQE_Selected = true;
				}
        	} 
        	// Otherwise, it remains false (either not in the map, or had only non-passing clusters)


			// ***************************************
			// True flag and event categorization
			//

			// did the event originate outside the water volume?
			extVol = External(nuvtxx, nuvtxy, nuvtxz);

			// did the event occur within the FV?
			trueFV = FidVol(nuvtxx, nuvtxy, nuvtxz, false);

			// Check background categoriesin an order that ensures mutual exclusivity

			// External
			if (extVol) {
                NCQE_EventCategory = 1;
            }

			// Out-of-FV
			else if (!extVol && !trueFV) {
                NCQE_EventCategory = 2;
            }

			// At this point, all remaining events are External == 0 && FV == 1.
            // These are events inside the Fiducial Volume.

			// SIGNAL (nuNCQE)
			else if (isNC && isQEL && ((nuPDG == 14) || (nuPDG == 12)) && !isMEC && (targetZ == 8)) {
				NCQE_MC_Signal = true;
				NCQE_EventCategory = 0;
			}

			// nubarNCQE
			else if (isNC && isQEL && ((nuPDG == -14) || (nuPDG == -12)) && !isMEC && (targetZ == 8)) {
                NCQE_EventCategory = 5;
            }

			// CC
			else if (!isNC) {
                NCQE_EventCategory = 3;
            }

			// NCother
            // External == 0 && FV == 1 && trueNC == 1 && NOT(QEL && TargetZ=8 && !MEC && NuPDG>0)
            // This is the default for all other NC events inside the FV.
            // Since all previous cases (External, Out-of-FV, nuNCQE, nubarNCQE, CC)
            // have been filtered out, any remaining event must be NCother.
            // (The explicit conditions for nuNCQE and nubarNCQE cover the 'NOT' part 
            // of the NCother definition).
            else if (isNC) { // Equivalent to trueNC == 1 and all other conditions failed
                NCQE_EventCategory = 4;
            }

        	// Fill the new branch entry
        	NCQEreco->Fill();
			NCQEsig->Fill();
			NCQEevtcat->Fill();		


			//XSec weight appendage (what MicroBooNE uses)
        	All_weight.insert(All_weight.end(), All0_weight->begin(), All0_weight->end());
			All_weight.insert(All_weight.end(), All1_weight->begin(), All1_weight->end());
			All_weight.insert(All_weight.end(), All2_weight->begin(), All2_weight->end());
			All_weight.insert(All_weight.end(), All3_weight->begin(), All3_weight->end());
			All_weight.insert(All_weight.end(), All4_weight->begin(), All4_weight->end());

			// lump all flux uncertainties together
			for(int j = 0; j < 1000; j++){
            	flux_All.push_back(
				flux_horncurrent->at(j) *
				flux_expskin->at(j) *
				flux_piplus->at(j) *
				flux_piminus->at(j) *
				flux_kplus->at(j) *
				flux_kminus->at(j) *
				flux_kzero->at(j) *
				flux_pionine->at(j) *
				flux_pionqe->at(j) *
				flux_piontot->at(j) *
				flux_nucine->at(j) *
				flux_nucqe->at(j) *
				flux_nuctot->at(j)
				);
        	}

			WAll->Fill();
        	WfAll->Fill();

        	All_weight.clear();
        	flux_All.clear();


    	}
		

		// *********************************
		// Cleanup
		//

        tTrig->Write("",TObject::kOverwrite);
        f->Close();
		delete f;

	}   // end of file loop
}

// done
