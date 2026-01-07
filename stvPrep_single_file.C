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

bool FidVol(double x, double y, double z){
        double radius   = 100.;  //cm
        double y_min    = -100.; //cm
        double y_max    = 100.;  //cm
        double z_center = 168.1; //cm
        double y_offset = 14.46; //cm
        if(y+y_offset > y_min && y+y_offset < y_max && radius > std::sqrt((z - z_center)*(z - z_center) + x*x)){
                return true;
        }
        else {
                return false;
        }
}

void stvPrep(){
        //Open file and trees
        TFile *f = new TFile("/exp/annie/data/users/jminock/standard_tank_ntuples/PhaseIITree_400k_ntuple.root","update");
        gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDataModel.so");
//      gSystem->Load("/exp/annie/app/users/jminock/ToolAnalysis/lib/libDict.so");
        gInterpreter->GenerateDictionary("map<string,vector<double>>", "map;string;vector");
        TTree *tTrig = (TTree*)f->Get("phaseIITriggerTree");

        //insert variables here
        int rnTrig;
        int rnMRD;
        int evNTrig;
        int evNMRD;
        int trigNum;
        int trigword, hasTank, hasMRD, tankMRDCoinc, noveto;
        int nrings;
        double nuE, nuvtxx, nuvtxy, nuvtxz, nupx, nupy, nupz;
        double fslpx, fslpy, fslpz, mcfslE;
        double fslvtxx, fslvtxy, fslvtxz;
        int isCC, isQEL, isRES, isDIS, isCOH, isMEC, mcfslpdg;
        double mcvtxx, mcvtxy, mcvtxz, mcdirx, mcdiry, mcdirz, mcangle, mcmuonE, mctanktracklength, mcmrdtracklength;
        int hasPi0, hasPiP, hasPiM, hasPiPC, hasPiMC, hasKP, hasKM, hasKPC, hasKMC;
        int numMRDTracks;
        int mcentersmrd, mcexitsmrd, mcpenetratesmrd;

        bool mcfv;
	double mcp, mcct;
	bool mc_no_mesons;
	bool recofv, recoMRD;
	double recop, recopc;

        int simpleflag, simplefv;
        double simpleenergy, simplecostheta, simplept, simplemrdenergy, simplemrdtrack, simpletanktrack;
        double simplevtxx, simplevtxy, simplevtxz, simplestopvtxx, simplestopvtxy, simplestopvtxz;
        double simplemrdstartx, simplemrdstarty, simplemrdstartz, simplemrdstopx, simplemrdstopy, simplemrdstopz;
        double rcsr, rcmr;
        double recovtxx, recovtxy, recovtxz;
        int julieflag, juliefv, nummrdlayers;
        double julieenergy, juliecostheta, juliept, juliemrdtrack;
        double julievtxx, julievtxy, julievtxz, juliestopvtxx, juliestopvtxy, juliestopvtxz;
        double juliemrdstartx, juliemrdstarty, juliemrdstartz, juliemrdstopx, juliemrdstopy, juliemrdstopz;
        vector<double>* MRDTrackAngle = new vector<double>();
        vector<double>* MRDTrackAngleError = new vector<double>();
        vector<double>* MRDPenetrationDepth = new vector<double>();
        vector<double>* MRDTrackLength = new vector<double>();
        vector<double>* MRDEntryPointRadius = new vector<double>();
        vector<double>* MRDEnergyLoss = new vector<double>();
        vector<double>* MRDEnergyLossError = new vector<double>();
        vector<double>* MRDTrackStartX = new vector<double>();
        vector<double>* MRDTrackStartY = new vector<double>();
        vector<double>* MRDTrackStartZ = new vector<double>();
        vector<double>* MRDTrackStopX = new vector<double>();
        vector<double>* MRDTrackStopY = new vector<double>();
        vector<double>* MRDTrackStopZ = new vector<double>();
        vector<bool>* MRDSide = new vector<bool>();
        vector<bool>* MRDStop = new vector<bool>();
        vector<bool>* MRDThrough = new vector<bool>();
//      vector<double>* MRDTrackLengthTrig = new vector<double>();
//      vector<double>* MRDTrackLengthMRD = new vector<double>();
        vector<double>* All_weight = new vector<double>();
	vector<double>* All0_weight = new vector<double>();
	vector<double>* All1_weight = new vector<double>();
	vector<double>* All2_weight = new vector<double>();
	vector<double>* All3_weight = new vector<double>();
	vector<double>* All4_weight = new vector<double>();
	vector<double>* All5_weight = new vector<double>();
	vector<double>* flux_All = new vector<double>();
	vector<double>* flux_horncurrent = new vector<double>();
	vector<double>* flux_expskin = new vector<double>();
	vector<double>* flux_piplus = new vector<double>();
	vector<double>* flux_piminus = new vector<double>();
	vector<double>* flux_kplus = new vector<double>();
	vector<double>* flux_kminus = new vector<double>();
	vector<double>* flux_kzero = new vector<double>();
	vector<double>* flux_pionine = new vector<double>();
	vector<double>* flux_pionqe = new vector<double>();
	vector<double>* flux_piontot = new vector<double>();
	vector<double>* flux_nucine = new vector<double>();
	vector<double>* flux_nucqe = new vector<double>();
	vector<double>* flux_nuctot = new vector<double>();
        map<string,vector<double>>* xsecweights = new map<string,vector<double>>();
        map<string,vector<double>>* fluxweights = new map<string,vector<double>>();
        vector<double>* TCV_weight = new vector<double>();
        //TBranch *bwgts = 0;
        //Set branch addresses
//      tTrig->SetBranchAddress("weight_All_UBGenie",&All_weight);
//      tTrig->SetBranchAddress("weight_TunedCentralValue_UBGenie",&TCV_weight);
//      TBranch *TCV = tTrig->Branch("weight_TunedCentralValue_UBGenie",&TCV_weight);

        tTrig->SetBranchAddress("runNumber",&rnTrig);
        tTrig->SetBranchAddress("eventNumber",&evNTrig);
        tTrig->SetBranchAddress("trigword",&trigword);
        tTrig->SetBranchAddress("HasTank",&hasTank);
        tTrig->SetBranchAddress("HasMRD",&hasMRD);
        tTrig->SetBranchAddress("TankMRDCoinc",&tankMRDCoinc);
        tTrig->SetBranchAddress("NoVeto",&noveto);
        tTrig->SetBranchAddress("numMRDTracks",&numMRDTracks);
        tTrig->SetBranchAddress("MRDTrackAngle",&MRDTrackAngle);
        tTrig->SetBranchAddress("MRDTrackAngleError",&MRDTrackAngleError);
        tTrig->SetBranchAddress("MRDPenetrationDepth",&MRDPenetrationDepth);
        tTrig->SetBranchAddress("MRDTrackLength",&MRDTrackLength);
        tTrig->SetBranchAddress("MRDEntryPointRadius",&MRDEntryPointRadius);
        tTrig->SetBranchAddress("MRDEnergyLoss",&MRDEnergyLoss);
        tTrig->SetBranchAddress("MRDEnergyLossError",&MRDEnergyLossError);
        tTrig->SetBranchAddress("MRDTrackStartX",&MRDTrackStartX);
        tTrig->SetBranchAddress("MRDTrackStartY",&MRDTrackStartY);
        tTrig->SetBranchAddress("MRDTrackStartZ",&MRDTrackStartZ);
        tTrig->SetBranchAddress("MRDTrackStopX",&MRDTrackStopX);
        tTrig->SetBranchAddress("MRDTrackStopY",&MRDTrackStopY);
        tTrig->SetBranchAddress("MRDTrackStopZ",&MRDTrackStopZ);
        tTrig->SetBranchAddress("MRDSide",&MRDSide);
        tTrig->SetBranchAddress("MRDStop",&MRDStop);
        tTrig->SetBranchAddress("MRDThrough",&MRDThrough);

        tTrig->SetBranchAddress("trueCC",&isCC);
        tTrig->SetBranchAddress("trueQEL",&isQEL);
        tTrig->SetBranchAddress("trueRES",&isRES);
        tTrig->SetBranchAddress("trueDIS",&isDIS);
        tTrig->SetBranchAddress("trueCOH",&isCOH);
        tTrig->SetBranchAddress("trueMEC",&isMEC);
        tTrig->SetBranchAddress("trueNuIntxVtx_X",&nuvtxx);
        tTrig->SetBranchAddress("trueNuIntxVtx_Y",&nuvtxy);
        tTrig->SetBranchAddress("trueNuIntxVtx_Z",&nuvtxz);
        tTrig->SetBranchAddress("trueNeutrinoEnergy",&nuE);
        tTrig->SetBranchAddress("trueNeutrinoMomentum_X",&nupx);
        tTrig->SetBranchAddress("trueNeutrinoMomentum_Y",&nupy);
        tTrig->SetBranchAddress("trueNeutrinoMomentum_Z",&nupz);
        tTrig->SetBranchAddress("truePi0",&hasPi0);
        tTrig->SetBranchAddress("truePiPlus",&hasPiP);
        tTrig->SetBranchAddress("truePiMinus",&hasPiM);
        tTrig->SetBranchAddress("truePiPlusCher",&hasPiPC);
        tTrig->SetBranchAddress("truePiMinusCher",&hasPiMC);
        tTrig->SetBranchAddress("trueKPlus",&hasKP);
        tTrig->SetBranchAddress("trueKMinus",&hasKM);
        tTrig->SetBranchAddress("trueKPlusCher",&hasKPC);
        tTrig->SetBranchAddress("trueKMinusCher",&hasKMC);
        tTrig->SetBranchAddress("trueFSLMomentum_X",&fslpx);
        tTrig->SetBranchAddress("trueFSLMomentum_Y",&fslpy);
        tTrig->SetBranchAddress("trueFSLMomentum_Z",&fslpz);
        tTrig->SetBranchAddress("trueFSLVtx_X",&fslvtxx);
        tTrig->SetBranchAddress("trueFSLVtx_Y",&fslvtxy);
        tTrig->SetBranchAddress("trueFSLVtx_Z",&fslvtxz);
        tTrig->SetBranchAddress("trueFSLEnergy",&mcfslE);
        tTrig->SetBranchAddress("trueFSLPdg",&mcfslpdg);
        tTrig->SetBranchAddress("triggerNumber",&trigNum);
//      tTrig->SetBranchAddress("XSecWeights",&xsecweights);
//      tTrig->SetBranchAddress("FluxWeights",&fluxweights);

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


        //Simple Reco
        tTrig->SetBranchAddress("simpleRecoFlag",&simpleflag);
        tTrig->SetBranchAddress("simpleRecoEnergy",&simpleenergy);
        tTrig->SetBranchAddress("simpleRecoVtxX",&simplevtxx);
        tTrig->SetBranchAddress("simpleRecoVtxY",&simplevtxy);
        tTrig->SetBranchAddress("simpleRecoVtxZ",&simplevtxz);
        tTrig->SetBranchAddress("simpleRecoStopVtxX",&simplestopvtxx);
        tTrig->SetBranchAddress("simpleRecoStopVtxY",&simplestopvtxy);
        tTrig->SetBranchAddress("simpleRecoStopVtxZ",&simplestopvtxz);
        tTrig->SetBranchAddress("simpleRecoCosTheta",&simplecostheta);
        tTrig->SetBranchAddress("simpleRecoPt",&simplept);
        tTrig->SetBranchAddress("simpleRecoFV",&simplefv);
        tTrig->SetBranchAddress("simpleRecoMrdEnergyLoss",&simplemrdenergy);
        tTrig->SetBranchAddress("simpleRecoTrackLengthInMRD",&simplemrdtrack);
        tTrig->SetBranchAddress("simpleRecoTrackLengthInTank",&simpletanktrack);
        tTrig->SetBranchAddress("simpleRecoMRDStartX",&simplemrdstartx);
        tTrig->SetBranchAddress("simpleRecoMRDStartY",&simplemrdstarty);
        tTrig->SetBranchAddress("simpleRecoMRDStartZ",&simplemrdstartz);
        tTrig->SetBranchAddress("simpleRecoMRDStopX",&simplemrdstopx);
        tTrig->SetBranchAddress("simpleRecoMRDStopY",&simplemrdstopy);
        tTrig->SetBranchAddress("simpleRecoMRDStopZ",&simplemrdstopz);

        //Ring Counting Reco
        tTrig->SetBranchAddress("RCSRPred",&rcsr);
        tTrig->SetBranchAddress("RCMRPred",&rcmr);

        //true Muon info from WCSim
        tTrig->SetBranchAddress("trueMuonEnergy",&mcmuonE);
        tTrig->SetBranchAddress("trueVtxX",&mcvtxx);
        tTrig->SetBranchAddress("trueVtxY",&mcvtxy);
        tTrig->SetBranchAddress("trueVtxZ",&mcvtxz);
        tTrig->SetBranchAddress("trueDirX",&mcdirx);
        tTrig->SetBranchAddress("trueDirY",&mcdiry);
        tTrig->SetBranchAddress("trueDirZ",&mcdirz);
        tTrig->SetBranchAddress("trueTrackLengthInWater",&mctanktracklength);
        tTrig->SetBranchAddress("trueTrackLengthInMRD",&mcmrdtracklength);
        tTrig->SetBranchAddress("trueAngle",&mcangle);
//      tTrig->SetBranchAddress("trueEntersMRD",&mcentersmrd);
//      tTrig->SetBranchAddress("trueExitsMRD",&mcexitsmrd);
//      tTrig->SetBranchAddress("truePenetratesMRD",&mcpenetratesmrd);
        TBranch *TFV    = tTrig->Branch("trueFV",&mcfv);
	TBranch *RFV    = tTrig->Branch("recoFV",&recofv);
	TBranch *TMuP   = tTrig->Branch("trueMuonMomentum",&mcp);
	TBranch *TCT    = tTrig->Branch("trueCosTheta",&mcct);
	TBranch *RMuP   = tTrig->Branch("simpleRecoMomentum",&recop);
	TBranch *RMuPC  = tTrig->Branch("simpleRecoMomentumCor",&recopc);
	TBranch *TnoPi  = tTrig->Branch("true_no_mesons",&mc_no_mesons);
	TBranch *RinMRD = tTrig->Branch("reco_contained_in_MRD",&recoMRD);
	TBranch *WAll = tTrig->Branch("weight_All_UBGenie",&All_weight);
	TBranch *WfAll = tTrig->Branch("weight_flux_all",&flux_All);

//      tMRD->SetBranchAddress("runNumber",&rnMRD);
//      tMRD->SetBranchAddress("eventNumber",&evNMRD);
//      tMRD->SetBranchAddress("MRDTrackLength",&MRDTrackLengthMRD);


        double muon_m = 105.66;
        Long64_t nentriesTrig = tTrig->GetEntries();
        std::cout << "TriggerTree: " << nentriesTrig << std::endl;
        //fill histograms
        for (Long64_t i = 0; i < nentriesTrig; i++) {
//      for (Long64_t i = 0; i < 2; i++) {
                tTrig->GetEntry(i);
                if(i%1000 == 0) std::cout << i << std::endl;
                bool inFV = FidVol(nuvtxx,nuvtxy,nuvtxz);
                bool simplerecoFV = FidVol(simplevtxx*100.,simplevtxy*100.,simplevtxz*100.);
                double fslp = std::sqrt(fslpx*fslpx + fslpy*fslpy + fslpz*fslpz);
                bool hasPi = ((hasPi0) || (hasPiP) || (hasPiM));
                bool hasVisPi = ((hasPi0) || (hasPiPC) || (hasPiMC));
                bool hasNonMuon = ((hasPi) || (hasKP) || (hasKM));
                bool hasVisNonMuon = ((hasVisPi) || (hasKPC) || (hasKMC));

                mcfv = inFV;
		recofv = simplerecoFV;
		mcp = std::sqrt(mcmuonE*mcmuonE - muon_m*muon_m);
		mcct = std::cos(mcangle*M_PI/180.);
		recop  = std::sqrt(simpleenergy*simpleenergy - muon_m*muon_m);
		recopc = std::sqrt(simpleenergy*simpleenergy - muon_m*muon_m)*1.25 - 137.;
		mc_no_mesons = !(hasVisNonMuon);
		recoMRD = numMRDTracks == 1 ? MRDStop->at(0) : false;

		//check all the weights have the same number of universes
/*		if(All0_weight->size() != 100 || All1_weight->size() != 100 || All2_weight->size() != 100 || All3_weight->size() != 100 || All4_weight->size() != 100) {
			std::cerr << "[ERROR] wrong universe size " << i << ": " << All0_weight->size() << " " << All1_weight->size() << " " << All2_weight->size() << " " << All3_weight->size() << " " << All4_weight->size() << std::endl;
			return false;
		}
		if(flux_horncurrent->size() != 1000 || flux_expskin->size() != 1000 || flux_piplus->size() != 1000 || flux_piminus->size() != 1000 || flux_kplus->size() != 1000 || flux_kminus->size() != 1000 || flux_kzero->size() != 1000 || flux_pionine->size() != 1000 || flux_pionqe->size() != 1000 || flux_piontot->size() != 1000 || flux_nucine->size() != 1000 || flux_nucqe->size() != 1000 || flux_nuctot->size() != 1000 ) {
			std::cerr << "[ERROR] wrong universe size " << i << ": " << flux_horncurrent->size() << " " << flux_expskin->size() << " " << flux_piplus->size() << " " << flux_piminus->size() << " " << flux_kplus->size() << " " << flux_kminus->size() << " " << flux_kzero->size() << " " << flux_pionine->size() << " " << flux_pionqe->size() << " " << flux_piontot->size() << " " << flux_nucine->size() << " " << flux_nucqe->size() << " " << flux_nuctot->size() << std::endl;
			return false;
		}
*/
		//XSec weight appendage
		All_weight->insert(All_weight->end(), All0_weight->begin(), All0_weight->end());
		All_weight->insert(All_weight->end(), All1_weight->begin(), All1_weight->end());
		All_weight->insert(All_weight->end(), All2_weight->begin(), All2_weight->end());
		All_weight->insert(All_weight->end(), All3_weight->begin(), All3_weight->end());
		All_weight->insert(All_weight->end(), All4_weight->begin(), All4_weight->end());

//		for(int j = 0; j < All0_weight->size(); j++){
//			All_weight->push_back(All0_weight->at(j));
//			if(All_weight->at(j) != All0_weight->at(j)) std::cout << "[ERROR] weight mismatch: " << i << " " << j << " " << All_weight->at(j) << " " << All0_weight->at(j) << std::endl;
//			std::cout << All0_weight->at(j) << std::endl;
//		}

		for(int j = 0; j < 1000; j++){
			flux_All->push_back(flux_horncurrent->at(j)*flux_expskin->at(j)*flux_piplus->at(j)*flux_piminus->at(j)*flux_kplus->at(j)*flux_kminus->at(j)*flux_kzero->at(j)*flux_pionine->at(j)*flux_pionqe->at(j)*flux_piontot->at(j)*flux_nucine->at(j)*flux_nucqe->at(j)*flux_nuctot->at(j));
		}
//All and flux_all

                TFV->Fill();
		RFV->Fill();
		TMuP->Fill();
		TCT->Fill();
		RMuP->Fill();
		RMuPC->Fill();
		TnoPi->Fill();
		RinMRD->Fill();
		WAll->Fill();
		WfAll->Fill();

		All_weight->clear();
		flux_All->clear();
	}


        tTrig->Write("",TObject::kOverwrite);
//      tTrig->ResetBranchAddresses();
        delete f;
}
