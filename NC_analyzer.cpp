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
#include "WienerSVDUnfolder.hh"

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

    std::cout << "\nIntializing File Properties Manager...\n" << std::endl;
	  auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "file_properties.txt" );
    std::cout << "\nFile Properties Manager initialized" << std::endl;

	auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "file_properties.txt" );


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
    double total_pot = mcc9->total_bnb_data_pot_;   // from Data (FakeData) file
    double total_pot = mcc9->total_bnb_data_pot_;   // from Data / FakeData file
    double integ_flux = integrated_numu_flux_in_FV( total_pot );
    double num_Ar = num_O_targets_in_FV();
    double conv_factor = (num_Ar * integ_flux)/1e38;

	std::cout << "Total POT = " << total_pot << std::endl;
	std::cout << "Integrated flux = " << integ_flux << std::endl;
	std::cout << "N targets = " << num_Ar << std::endl;
	std::cout << "Conversion factor = " << conv_factor << std::endl;

    // covariance matrix and event rates
    std::cout << "\nGrabbing covariance matrix...\n" << std::endl;

    // Evaluate the total and partial covariance matrices in reco space
    auto matrix_map = syst.get_covariances();

	std::cout << "\nStarting the unfolding -----------------" << std::endl;

    // Perform background subtraction and unfolding to get a measurement of event
    // counts in (regularized) true space
    std::unique_ptr< Unfolder > unfolder (new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kFirstDeriv ) );
    auto result = unfolder->unfold( *mcc9 );

    // unfolding diagnostics
    TMatrixD* signal = result.unfolded_signal_.get();
    TMatrixD* S = result.response_matrix_.get();
    TMatrixD* cov = result.cov_matrix_.get();
    std::cout << "\n\n  === Unfolding Output ===\n\n";

    if (signal) {
      std::cout << "  [Unfolded signal: true bin counts]\n";
      for (int i = 0; i < signal->GetNrows(); ++i) {
        std::cout << "    True bin " << i << ": " << (*signal)(i, 0) << " events\n";
      }
    } else {
      std::cout << "  [ERROR] unfolded_signal_ is null\n";
    }

    if (S) {
      std::cout << "  [Smearceptance matrix S (Reco x True)]\n";      // efficiency
      for (int i = 0; i < S->GetNrows(); ++i) {
        for (int j = 0; j < S->GetNcols(); ++j) {
          std::cout << "    S(" << i << "," << j << ") = " << (*S)(i,j) << "\n";
        }
      }
    } else {
      std::cout << "  [ERROR] response_matrix_ (smearceptance) is null\n";
    }

    if (cov) {
      std::cout << "  [Unfolded covariance: sqrt(diagonal elements)]\n";
      for (int i = 0; i < cov->GetNrows(); ++i) {
        std::cout << "    Bin " << i << ": ±" << std::sqrt((*cov)(i,i)) << "\n";
      }
    } else {
      std::cout << "  [ERROR] cov_matrix_ is null\n";
    }

    // end diagnostics


    std::cout << "\nUnfolding completed -----------------" << std::endl;
    std::cout << "\nPost-processing covariance matrices..\n" << std::endl;

    // The Error Propagation Matrix (A_c) tells us how a shift in Reco-space
    // translates to a shift in the Unfolded Truth-space.
    const TMatrixD& err_prop = *result.err_prop_matrix_;
    TMatrixD err_prop_tr( TMatrixD::kTransposed, err_prop );

    // This map will store your unfolded Truth-space systematic matrices
    std::map< std::string, std::unique_ptr<TMatrixD> > unfolded_cov_matrix_map;

    for ( const auto& matrix_pair : *matrix_map ) {
        const std::string& matrix_key = matrix_pair.first;
        // This is the Reco-space matrix (e.g., "flux", "genie")
        auto temp_cov_mat = matrix_pair.second.get_matrix();

        // Matrix multiplication: Cov_truth = A_c * Cov_reco * A_c^T
        TMatrixD temp_mat( *temp_cov_mat, TMatrixD::EMatrixCreatorsOp2::kMult, err_prop_tr );
        unfolded_cov_matrix_map[ matrix_key ] = std::make_unique< TMatrixD >(
            err_prop, TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat );

        // Print the breakdown
        double diag_val = (*unfolded_cov_matrix_map[matrix_key])(0,0);
        double fractional_err = std::sqrt(diag_val) / (*signal)(0,0);
        
        std::cout << "  Systematic [" << std::setw(15) << matrix_key << "]: ±" 
                  << std::sqrt(diag_val) << " events (" 
                  << fractional_err * 100.0 << "%)" << std::endl;
    }

    // Final Cross Section Calculation
    double val_xsec = (*signal)(0,0) / conv_factor; 
    double err_xsec = std::sqrt((*cov)(0,0)) / conv_factor;

    std::cout << "\n\n================================================" << std::endl;
    std::cout << "FINAL MEASUREMENT:" << std::endl;
    std::cout << "Cross Section: " << val_xsec << " x 10^-38 cm^2" << std::endl;
    std::cout << "Total Error:   ±" << err_xsec << " (" << (err_xsec/val_xsec)*100.0 << "%)" << std::endl;
    std::cout << "================================================" << std::endl;


    // 3. VALIDATION STANZA (Truth Comparison)
    std::cout << "\n\n--- Validation Stanza ---\n" << std::endl;

    // This gets the 1D projection used for the unfolding input
    TH1D* reco_data_hist = mcc9->cv_universe().hist_reco_.get(); 
    std::cout << "Raw selected data events (reco bins): " << reco_data_hist->Integral() << std::endl;

    // Get the GENIE CV Truth prediction (The "MicroBooNETune" equivalent)
    const Universe& cv_univ = mcc9->cv_universe();
    std::cout << "GENIE truth prediction (true bins):\n" << std::endl;
    for (int b = 0; b < num_true_bins; ++b) {
        // hist_true_ contains the actual number of signal events predicted in the FV
        double genie_pred = cv_univ.hist_true_->GetBinContent( b + 1 );
        double genie_err  = cv_univ.hist_true_->GetBinError(b + 1);
        std::cout << "  True bin " << b << ": " << genie_pred << " ± " << genie_err << " events" << std::endl;

        if (b == 0) {  // if true signal, we want the predicted XS
          double model_xs     = genie_pred / conv_factor;
          double model_xs_err = genie_err / conv_factor;
          // Calculate the "True" Cross Section from the model to see if they match
          std::cout << "-------------------------------------------------------" << std::endl;
          std::cout << "  >>> GENIE MODEL PREDICTION <<<" << std::endl;
          std::cout << "  Events (CV): " << genie_pred << " +/- " << genie_err << " (MC Stats)" << std::endl;
          std::cout << "  Model XS:    " << model_xs << " +/- " << model_xs_err << " [10^-38 cm^2]" << std::endl;
          std::cout << "-------------------------------------------------------\n" << std::endl;
        }
    }

    // end function
    std::cout << "\n" <<std::endl;
	std::cout << "\nStarting the unfolding -----------------" << std::endl;
	// Perform background subtraction and unfolding to get a measurement of event
  // counts in (regularized) true space
	//std::unique_ptr< Unfolder > unfolder (new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kFirstDeriv ) );
    //auto result = unfolder->unfold( *mcc9 );

	std::unique_ptr< Unfolder > unfolder (new WienerSVDUnfolder( true, WienerSVDUnfolder::RegularizationMatrixType::kIdentity ) );
    auto result = unfolder->unfold( *mcc9 );

    //double evts = result.unfolded_signal_->operator()( 0, 0 );
    //double error = std::sqrt( std::max(0., result.cov_matrix_->operator()( 0, 0 )) );


  // Begin Burke Edits: unfolding diagnostics
  TMatrixD* signal = result.unfolded_signal_.get();
  TMatrixD* S = result.response_matrix_.get();
  TMatrixD* cov = result.cov_matrix_.get();
  std::cout << "\n  === Unfolding Output ===\n";

  if (signal) {
    std::cout << "  [Unfolded signal: true bin counts]\n";
    for (int i = 0; i < signal->GetNrows(); ++i) {
      std::cout << "    True bin " << i << ": " << (*signal)(i, 0) << " events\n";
    }
  } else {
    std::cout << "  [ERROR] unfolded_signal_ is null\n";
  }

  if (S) {
    std::cout << "  [Smearceptance matrix S (Reco x True)]\n";
    for (int i = 0; i < S->GetNrows(); ++i) {
      for (int j = 0; j < S->GetNcols(); ++j) {
        std::cout << "    S(" << i << "," << j << ") = " << (*S)(i,j) << "\n";
      }
    }
  } else {
    std::cout << "  [ERROR] response_matrix_ (smearceptance) is null\n";
  }

  if (cov) {
    std::cout << "  [Unfolded covariance: sqrt(diagonal elements)]\n";
    for (int i = 0; i < cov->GetNrows(); ++i) {
      std::cout << "    Bin " << i << ": ±" << std::sqrt((*cov)(i,i)) << "\n";
    }
  } else {
    std::cout << "  [ERROR] cov_matrix_ is null\n";
  }

  // End Burke Edits
  //CrossSectionResult xsec( result );

  std::cout << "Unfolding completed -----------------" << std::endl;
  std::cout << "\nPost-processing covariance matrices.." << std::endl;

  // The Error Propagation Matrix (A_c) tells us how a shift in Reco-space
  // translates to a shift in the Unfolded Truth-space.
  const TMatrixD& err_prop = *result.err_prop_matrix_;
  TMatrixD err_prop_tr( TMatrixD::kTransposed, err_prop );

  // This map will store your unfolded Truth-space systematic matrices
  std::map< std::string, std::unique_ptr<TMatrixD> > unfolded_cov_matrix_map;

  for ( const auto& matrix_pair : *matrix_map ) {
      const std::string& matrix_key = matrix_pair.first;
      // This is the Reco-space matrix (e.g., "flux", "genie")
      auto temp_cov_mat = matrix_pair.second.get_matrix();

      // Matrix multiplication: Cov_truth = A_c * Cov_reco * A_c^T
      TMatrixD temp_mat( *temp_cov_mat, TMatrixD::EMatrixCreatorsOp2::kMult, err_prop_tr );
      unfolded_cov_matrix_map[ matrix_key ] = std::make_unique< TMatrixD >(
          err_prop, TMatrixD::EMatrixCreatorsOp2::kMult, temp_mat );

      // Print the breakdown for your presentation!
      double diag_val = (*unfolded_cov_matrix_map[matrix_key])(0,0);
      double fractional_err = std::sqrt(diag_val) / (*signal)(0,0);
      
      std::cout << "  Systematic [" << std::setw(15) << matrix_key << "]: ±" 
                << std::sqrt(diag_val) << " events (" 
                << fractional_err * 100.0 << "%)" << std::endl;
  }

  // Final Cross Section Calculation (Correcting the / vs * issue)
  double val_xsec = (*signal)(0,0) / conv_factor; 
  double err_xsec = std::sqrt((*cov)(0,0)) / conv_factor;

  std::cout << "\n================================================" << std::endl;
  std::cout << "FINAL MEASUREMENT (True Bin 0):" << std::endl;
  std::cout << "Cross Section: " << val_xsec << " x 10^-38 cm^2" << std::endl;
  std::cout << "Total Error:   ±" << err_xsec << " (" << (err_xsec/val_xsec)*100.0 << "%)" << std::endl;
  std::cout << "================================================" << std::endl;

	// 3. VALIDATION STANZA (Truth Comparison)
  std::cout << "\n--- Validation Stanza ---" << std::endl;
  
  // Get the "Raw" reco data counts (N_obs)
  //TH1D* reco_data_hist = mcc9->data_hists2d_.at( NFT::kOnBNB ).get();
  //std::cout << "Raw selected data events (reco bins): " << reco_data_hist->Integral() << std::endl;

	// This gets the 1D projection used for the unfolding input
	TH1D* reco_data_hist = mcc9->cv_universe().hist_reco_.get(); 
	std::cout << "Raw selected data events (reco bins): " << reco_data_hist->Integral() << std::endl;

  
	// Get the GENIE CV Truth prediction (The "MicroBooNETune" equivalent)
  const Universe& cv_univ = mcc9->cv_universe();
  std::cout << "GENIE truth prediction (true bins):\n" << std::endl;
  for (int b = 0; b < num_true_bins; ++b) {
      // hist_true_ contains the actual number of signal events predicted in the FV
      double genie_pred = cv_univ.hist_true_->GetBinContent( b + 1 );
      double genie_err  = cv_univ.hist_true_->GetBinError(b + 1);
      std::cout << "  True bin " << b << ": " << genie_pred << " ± " << genie_err << " events" << std::endl;

      if (b == 0) {  // if true signal, we want the predicted XS
        double model_xs     = genie_pred / conv_factor;
        double model_xs_err = genie_err / conv_factor;
        // Calculate the "True" Cross Section from the model to see if they match
        std::cout << "-------------------------------------------------------" << std::endl;
        std::cout << "  >>> GENIE MODEL PREDICTION <<<" << std::endl;
        std::cout << "  Events (CV): " << genie_pred << " +/- " << genie_err << " (MC Stats)" << std::endl;
        std::cout << "  Model XS:    " << model_xs << " +/- " << model_xs_err << " [10^-38 cm^2]" << std::endl;
        std::cout << "-------------------------------------------------------\n" << std::endl;
      }
  }
}


// and 4 to go
int main(int argc, char* argv[]) {
   XS_extractor();
   return 0;
}

// done
