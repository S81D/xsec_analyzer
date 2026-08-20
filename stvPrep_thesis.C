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
        double radius   = 70.;  //cm
        double y_min    = -100.; //cm
        double y_max    = 100.;  //cm
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
    const double bunch_sigma = 3.59;    // determined based on MC bunch width (update to include BRF fitting)
    const double avg_diff = 0.5081;     // weighted avg between NCQE and CCinc bunch positions

	// bunch variation (+/- 0.5ns)
	//const double avg_diff = -0.5081;

	// CCinc centers extracted from bunch fitting (update if needed)
    static const double CC_centers[] = {
    				30.55620074944292, 50.23596632813184, 68.03027546306257, 87.65149609251257, 106.38480577249247,
                    125.28116343989417, 144.32309666797107, 163.1397586239191, 181.97778394194054, 201.06274310271073,
                    220.16343853423476, 238.30288809166953, 258.38514881200354, 276.8146398611303, 296.0884356855441,
                    314.35016323142816, 333.72419652049575, 352.9437944077438, 371.260196759933, 390.08407553949974,
                    409.44913114396746, 427.94715055258985, 447.770917943728, 465.9339836950206, 484.7186644716273,
                    504.4287698808072, 523.0814143742606, 542.3141201678939, 561.2254420066273, 579.8111323601904,
                    597.9945900965694, 618.0700927506651, 636.5621477504429, 655.7874418315001, 674.4764151744993,
                    693.6905256979796, 712.4819567023824, 731.0362230741698, 750.2744038347329, 769.2102243656747,
                    788.2453849590231, 807.0011262101622, 825.9067358626738, 845.3047059731147, 864.0653260882543,
                    883.0323070032683, 901.6429155761502, 920.6808374077715, 939.6712240835648, 958.6990235603214,
                    977.72856819438, 996.2787338241341, 1015.08703056768, 1034.3375283354696, 1053.9759092673044,
                    1072.052984198524, 1091.129322973113, 1110.3389344568238, 1128.965140612168, 1148.3985550156349,
                    1166.982403872106, 1186.220604492118, 1204.4016122615021, 1223.770899436999, 1242.691800785879,
                    1261.9422162459668, 1280.3202276072477, 1299.702530449622, 1318.214516959161, 1337.6134384701843,
                    1356.5308849686603, 1375.0603202255497, 1393.91335999822, 1413.0147004160135, 1431.9961602521892,
                    1451.5391995383018, 1470.035773429933, 1488.7824718348304, 1507.6311809765648, 1526.7701548484445,
                   	1545.2663160965799
	};

    const int N = sizeof(CC_centers) / sizeof(double);

	// first 1/3 of the spill
	//const int max_bunch_index = 25;  // matches centers[25] in your Python
    //if (cluster_time > (CC_centers[max_bunch_index] - avg_diff + bunch_sigma)) {
    //    return false;
    //}

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

	//const double DIRT_SCALE = 0.18;   // based on sideband

	// ********************************
	// load libraries (modify this path to any compiled toolanalysis directory)
	gSystem->Load("/exp/annie/app/users/doran/NC_CC_Nov_6_2025/lib/libDataModel.so");
	gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");


	// loop for executing over N files
	for (int i = 0; i < 4000; ++i) {       // change 1000 to N to match N input files
		std::string file =
			"/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/stv/MC_" + 
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

			// 6. Bunch cutting (3.59ns sigma)
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

		//double dirt_muon = 1.0;             // apply weighting to external events		

		// total weights
		std::vector<double> All_weight;
		std::vector<double> flux_All;

		// branch names
		TBranch *NCQEsig    = tTrig->Branch("NCQE_MC_Signal",&NCQE_MC_Signal, "NCQE_MC_Signal/O");
		TBranch *NCQEevtcat = tTrig->Branch("NCQE_EventCategory",&NCQE_EventCategory, "NCQE_EventCategory/I");
		TBranch *NCQEreco   = tTrig->Branch("NCQE_Selected", &NCQE_Selected, "NCQEreco/O");
		TBranch *WAll       = tTrig->Branch("weight_All_UBGenie",&All_weight);
        TBranch *WfAll      = tTrig->Branch("weight_flux_all",&flux_All);
		//TBranch *DirtMu     = tTrig->Branch("DirtMu",&dirt_muon);
	

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
			//dirt_muon = -1;

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
				//dirt_muon = DIRT_SCALE;
            }

			// Out-of-FV
			else if (!extVol && !trueFV) {
                NCQE_EventCategory = 2;
				//dirt_muon = 1.0;
            }

			// At this point, all remaining events are External == 0 && FV == 1.
            // These are events inside the Fiducial Volume.

			// SIGNAL (nuNCQE)
			else if (isNC && isQEL && ((nuPDG == 14) || (nuPDG == 12)) && !isMEC && (targetZ == 8)) {
				NCQE_MC_Signal = true;
				NCQE_EventCategory = 0;
				//dirt_muon = 1.0;
			}

			// nubarNCQE
			else if (isNC && isQEL && ((nuPDG == -14) || (nuPDG == -12)) && !isMEC && (targetZ == 8)) {
                NCQE_EventCategory = 5;
				//dirt_muon = 1.0;
            }

			// CC
			else if (!isNC) {
                NCQE_EventCategory = 3;
				//dirt_muon = 1.0;
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
				//dirt_muon = 1.0;
            }

        	// Fill the new branch entry
        	NCQEreco->Fill();
			NCQEsig->Fill();
			NCQEevtcat->Fill();		
			//DirtMu->Fill();

			//XSec weight appendage (what MicroBooNE uses)
        	All_weight.insert(All_weight.end(), All0_weight->begin(), All0_weight->end());

			// lump all flux uncertainties together (loop over number of multisim)
			for(int j = 0; j < 100; j++){
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
