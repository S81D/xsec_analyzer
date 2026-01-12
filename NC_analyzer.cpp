/*Author: S. Doran <doran@iastate.edu>
 *
 * Usage: ./NC_analyzer
 *
 */

#include <iostream>
#include <iomanip>
#include <string> 
#include <vector>
#include <map>
#include <cmath>
#include "TFile.h"
#include "TMatrixD.h"
#include "TH1D.h"
#include "MCC9SystematicsCalculator.hh"
#include "includes/AnnieGeometryTools.hh"

void XS_extractor() {

    // 1 for the money
    std::cout << "\n______________________________________________" << std::endl;
    std::cout << "\nExecuting script..." << std::endl;

    // ....................................
    // Initialization
    //

    // Systematics calculator object
    std::cout << "\nIntializing Systematics Calculator...\n" << std::endl;
    auto* mcc9 = new MCC9SystematicsCalculator(
        "output.root",
        "systcalc.conf" );
    const auto &syst = *mcc9;
    std::cout << "\nSystematics Calculator initialized" << std::endl;


    // 2 for the show
    // ....................................
    // Execution
    //
    
    // binning
    int num_true_bins = mcc9->true_bins_.size();
    int num_reco_bins = mcc9->reco_bins_.size();
    std::cout<<"\nNum true bins = "<< num_true_bins << std::endl;
    std::cout << "Num reco bins = " << num_reco_bins << std::endl;


    // Grab conversion factor for going from event rate --> XS
    std::cout << "\nPreparing conversion factor...\n" << std::endl;

    // XS equation: (N-B) / [eff * N_targets * Phi]
    // conversion factor = N_targets * Phi [cm^2 / nucleon], then divide by 10^-38
    double total_pot = mcc9->total_bnb_data_pot_;   // from Data / FakeData file
    double integ_flux = integrated_numu_flux_in_FV( total_pot );
    double num_Ar = num_O_targets_in_FV();
    double conv_factor = (num_Ar * integ_flux)/1e38;

    // covariance matrix and event rates
    std::cout << "\nGrabbing covariance matrix...\n" << std::endl;

    auto true_signal = syst.get_cv_true_signal();    // POT-scaled GENIE CV true signal (NC truth events)
    auto reco_signal = syst.get_cv_reco_signal();    // POT-scaled GENIE (weighted) total reco

    auto meas = syst.get_measured_events();
    const auto& data_signal = meas.reco_signal_;     // data signal estimator (how many reconstructed signal events we have)
    const auto& data_covmat = meas.cov_matrix_;      // full covariance matrix (stat + syst) on observed reco events

    // Evaluate the total and partial covariance matrices in reco space
    auto* matrix_map_ptr = syst.get_covariances().release();
    auto& matrix_map = *matrix_map_ptr;

    // assumes single bin XS
    double N_true_MC   = true_signal->operator()(0,0);    // GENIE truth NCQE events in FV
    double N_reco_MC   = reco_signal->operator()(0,0);    // GENIE‐predicted reconstructed NCQE events
    double N_data      = meas.reco_signal_->operator()(0,0);  // Estimator of reconstructed NCQE in data
    double Var_data    = (*meas.cov_matrix_)(0,0);        // full covariance
    double sigma_data  = std::sqrt(Var_data);             // total uncertainty


    // correct N_data by the GENIE efficiency (takes you from reconstructed --> true); this is a simple scalar correction as its a 1D measurement
    double eff = N_reco_MC / N_true_MC;

    // XS = N / (eff * conv_factor)

    // data cross section
    double xs_data = N_data / (eff * conv_factor);
    double xs_data_err = sigma_data / (eff * conv_factor);

    // GENIE truth XS (with stat uncertainty)
    double sigma_mc_truth = std::sqrt(
        matrix_map.at("MCstats").cov_matrix_->GetBinContent(1,1)
    );
    double xs_genie_err = sigma_mc_truth / conv_factor;
    double xs_genie = N_true_MC / conv_factor;


    // Thanks Gemini for such a lovely print out table :)

    // 3 to get ready
    std::cout << "\n\n\n" << std::string(60, '=') << std::endl;
    std::cout << std::setw(40) << std::left << " NEUTRAL CURRENT CROSS SECTION EXTRACTION" << std::endl;
    std::cout << std::string(60, '=') << std::endl;

    std::cout << "\n[1] INPUT STATISTICS & EFFICIENCY" << std::endl;
    std::cout << std::setw(35) << std::left << "  - Total Data POT:" << std::scientific << std::setprecision(2) << total_pot << std::endl;
    std::cout << std::setw(35) << std::left << "  - Integrated Flux (Phi):" << integ_flux << " cm^-2" << std::endl;
    std::cout << std::setw(35) << std::left << "  - Targets in FV (N_t):" << num_Ar << std::endl;
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << std::setw(35) << std::left << "  - GENIE True Signal (N_true): " << N_true_MC << std::endl;
    std::cout << std::setw(35) << std::left << "  - GENIE Reco Signal (N_reco): " << N_reco_MC << std::endl;
    std::cout << std::setw(35) << std::left << "  - GENIE Efficiency: " << (eff * 100.0) << " %" << std::endl;
    std::cout << std::setw(35) << std::left << "  - Data Signal Estimator (N_data): " << N_data << " events" << std::endl;

    std::cout << "\n[2] UNCERTAINTY BREAKDOWN (on N_data)" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    std::cout << std::setw(25) << std::left << "  Source" << std::setw(12) << "Sigma" << std::setw(10) << "Frac Err" << std::endl;
    std::cout << std::string(50, '-') << std::endl;
    
    for (const auto& [name, cov] : matrix_map) {
        double sigma = std::sqrt(cov.cov_matrix_->GetBinContent(1,1));
        double frac = (N_data > 0) ? (sigma / N_data) * 100.0 : 0.0;
        
        std::cout << "  " << std::setw(23) << std::left << name 
                  << std::setw(12) << std::setprecision(3) << sigma 
                  << std::setprecision(1) << frac << " %" << std::endl;
    }
    std::cout << std::string(50, '-') << std::endl;
    std::cout << "  " << std::setw(23) << "TOTAL UNCERTAINTY" 
              << std::setw(12) << sigma_data 
              << std::setprecision(1) << (sigma_data/N_data)*100.0 << " %" << std::endl;

    std::cout << "\n[3] CROSS SECTION EQUATION\n" << std::endl;
    std::cout << "      N_data" << std::endl;
    std::cout << "XS = ------------------------" << std::endl;
    std::cout << "      eps * N_t * Phi" << std::endl;
    
    // Switch to scientific just for this substitution to handle large Nt and Phi
    std::cout << "\n      " << std::scientific << std::setprecision(3) << N_data << std::endl;
    std::cout << "XS = ----------------------------------------------------" << std::endl;
    std::cout << "      (" << eff << ") * (" << num_Ar << ") * (" << integ_flux << ")" << std::endl;
    
    // Switch back to fixed for the Final Results section
    std::cout << std::fixed;

    std::cout << "\n" << std::string(60, '=') << std::endl;
    std::cout << " FINAL RESULTS [10^-38 cm^2 / Oxygen]" << std::endl;
    std::cout << std::string(60, '-') << std::endl;
    
    // Convert to 10^-38 units for display if desired, or keep raw
    std::cout << "  DATA  XS: " << std::fixed << std::setprecision(3) << xs_data << " +/- " << xs_data_err << std::endl;
    std::cout << "  GENIE XS: " << xs_genie << " +/- " << xs_genie_err << " (MC stat)" << std::endl;
    
    double pull = (xs_data - xs_genie) / xs_data_err;
    std::cout << "  Data/MC Ratio: " << (xs_data / xs_genie) << std::endl;
    std::cout << "  Agreement:     " << std::abs(pull) << " sigma" << std::endl;
    std::cout << std::string(60, '=') << "\n" << std::endl;

}


// and 4 to go
int main(int argc, char* argv[]) {
   XS_extractor();
   return 0;
}