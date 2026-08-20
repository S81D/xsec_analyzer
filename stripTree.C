#include "TFile.h"
#include "TTree.h"
#include <iostream>

void stripTree() {

    const char* infile =
		"/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/GENIE-Empirical-4.417e20-ntuple_dirt_downsampled.root";

    const char* outfile =
        "/exp/annie/data/users/doran/GENIE_reweight_MC_ToolChain/GENIE-Empirical-4.417e20-ntuple_dirt_downsampled_stripped.root";

    // Open input file
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error opening input file" << std::endl;
        return;
    }

    // Get trigger tree
    TTree *T = (TTree*)fin->Get("phaseIITriggerTree");
    if (!T) {
        std::cerr << "Error: phaseIITriggerTree not found" << std::endl;
        fin->Close();
        return;
    }

    // Disable ALL branches
    T->SetBranchStatus("*", 0);

    // ----------------------------
    // Enable ONLY desired branches
    // ----------------------------

    // GENIE / cross-section weights
    //T->SetBranchStatus("weight_All0_UBGenie", 1);
    
	//T->SetBranchStatus("weight_All_UBGenie", 1);
    T->SetBranchStatus("weight_TunedCentralValue_UBGenie", 1);
    
	//T->SetBranchStatus("weight_AxFFCCQEshape_UBGenie", 1);
    //T->SetBranchStatus("weight_DecayAngMEC_UBGenie", 1);
    //T->SetBranchStatus("weight_NormCCCOH_UBGenie", 1);
    //T->SetBranchStatus("weight_NormNCCOH_UBGenie", 1);
    //T->SetBranchStatus("weight_RPA_CCQE_UBGenie", 1);
    
	//T->SetBranchStatus("weight_RootinoFix_UBGenie", 1);
    
	//T->SetBranchStatus("weight_ThetaDelta2NRad_UBGenie", 1);
    //T->SetBranchStatus("weight_Theta_Delta2Npi_UBGenie", 1);
    //T->SetBranchStatus("weight_VecFFCCQEshape_UBGenie", 1);
    //T->SetBranchStatus("weight_XSecShape_CCMEC_UBGenie", 1);

    // Flux systematics
    //T->SetBranchStatus("weight_flux_all", 1);
    
	//T->SetBranchStatus("weight_horncurrent_FluxUnisim", 1);
    //T->SetBranchStatus("weight_expskin_FluxUnisim", 1);
    //T->SetBranchStatus("weight_pioninexsec_FluxUnisim", 1);
    //T->SetBranchStatus("weight_piontotxsec_FluxUnisim", 1);
    //T->SetBranchStatus("weight_pionqexsec_FluxUnisim", 1);
    //T->SetBranchStatus("weight_nucleoninexsec_FluxUnisim", 1);
    //T->SetBranchStatus("weight_nucleontotxsec_FluxUnisim", 1);
    //T->SetBranchStatus("weight_nucleonqexsec_FluxUnisim", 1);

    // Hadron production weights
    //T->SetBranchStatus("weight_piplus_PrimaryHadronSWCentralSplineVariation", 1);
    //T->SetBranchStatus("weight_piminus_PrimaryHadronSWCentralSplineVariation", 1);
    //T->SetBranchStatus("weight_kminus_PrimaryHadronNormalization", 1);
    //T->SetBranchStatus("weight_kplus_PrimaryHadronFeynmanScaling", 1);
    //T->SetBranchStatus("weight_kzero_PrimaryHadronSanfordWang", 1);

    // NCQE flags
    T->SetBranchStatus("NCQE_MC_Signal", 1);
    T->SetBranchStatus("NCQE_EventCategory", 1);
    T->SetBranchStatus("NCQE_Selected", 1);

    // ----------------------------
    // Clone & write tree
    // ----------------------------

    TFile *fout = TFile::Open(outfile, "RECREATE");
    fout->cd();

    TTree *Tnew = T->CloneTree();
    Tnew->Write();

    fout->Close();
    fin->Close();

    std::cout << "Successfully wrote phaseIITriggerTree with selected branches to:\n  "
              << outfile << std::endl;
}

