/*Author: S. Doran <doran@iastate.edu>
 *
 * Usage: ./NC_analyzer files.txt
 *
 */

// Modified from chi_square_cc0pi_christian.cpp for the NCQE (single bin) total flux integrated cross section

#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <set>
#include "TChain.h"
#include "TFile.h"
#include "TParameter.h"
#include "TTree.h"
#include "TVector3.h"
#include "EventCategory.hh"
#include "TreeUtils.hh"
#include <fstream>
#include <sstream>
#include "TCanvas.h"
#include "TH2D.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TApplication.h"
#include "TPaveLabel.h"
#include <cassert>
#include <set>
#include <vector>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLegend.h>
#include <TNtuple.h>
#include <TPad.h>
#include "TColor.h"
#include "TInterpreter.h"
#include <algorithm>
#include "FilePropertiesManager.hh"
#include "MCC9SystematicsCalculator.hh"
#include "TMatrixT.h"
#include "WienerSVDUnfolder.hh"
#include <iomanip>
#include "TMatrixD.h"
#include "TLatex.h"
#include "RooStats/RooStatsUtils.h"
#include "TDecompChol.h"
#include "TDecompSVD.h"
#include "includes/PlotUtils.hh"
#include "includes/AnnieGeometryTools.hh"


// function will:
// - extract event rates (CV universe, data, ext data)
// - extract total stat + syst uncertainty on event rate
// - extract fractional uncertainties
// - compute cross section and full uncertainties


// GENIE CV (truth) --> GENIE CV signal events
// GENIE CV prediction --> for the reco

// ExtBNB --> Off-beam background (Off-beam (ExtBNB))

// OnBNB
// - Fake data (GENIE CV)
// - Fake data (NuWro)
// - Data

void Test_plot() {

    std::cout << "\nExecuting Test_plot()..." << std::endl;

    // ....................................
    // Initialization
    //

    gROOT->SetBatch(false);
    gStyle->SetOptStat(0);

    // PDF output file names
    char text_title_pdf[2024];
    char RootName[1024];
	  sprintf(text_title_pdf, "NC_univmake_plots.pdf(");
    std::string text_title_pdf_string(text_title_pdf);

    // load file_properties.txt
    auto& fpm = FilePropertiesManager::Instance();
    fpm.load_file_properties( "file_properties.txt" );
    std::cout << "\nfile_properties.txt loaded" << std::endl;

    // Systematics calculator object
    std::cout << "\nIntializing Systematics Calculator...\n" << std::endl;
    auto* mcc9 = new MCC9SystematicsCalculator(
        "output.root",
        "systcalc.conf" );
    const auto &syst = *mcc9;
    std::cout << "\nSystematics Calculator initialized" << std::endl;   

    // ....................................
    // Execution
    //
    
    // binning
    std::cout<<"\nNum true bins = "<< mcc9->true_bins_.size() << std::endl;
    std::cout << "Num reco bins = " << mcc9->reco_bins_.size() << std::endl;

    std::cout << "\nPreparing conversion factor...\n" << std::endl;

    // conversion factor
    double total_pot = mcc9->total_bnb_data_pot_;
    std::cout<<"Total POT: "<< total_pot << std::endl;
    double integ_flux = integrated_numu_flux_in_FV( total_pot );
    std::cout << "Integrated neutrino flux: " << integ_flux << std::endl;
    double num_O = num_O_targets_in_FV();
    std::cout << "N_targets: " << num_O << std::endl;

    // XS equation: (N-B) / [eff * N_targets * Phi]
    // conversion factor = N_targets * Phi [cm^2 / nucleon], then divide by 10^-38
    double conv_factor = (num_O * integ_flux)/1e38;
    std::cout << "conv_factor: " << conv_factor << std::endl;


    // from MicroBooNE (leave uncommented but good for validating we're doing it right)

    //int num_bins = syst.get_num_signal_true_bins();

    // Get CV universe
    //const Universe* cv_univ = &syst.cv_universe();
    //std::cout << "\nCV universe counts per bin:" << std::endl;
    //for (int b = 0; b < num_bins; ++b) {
    //    double true_evts = cv_univ->hist_true_->GetBinContent(b+1);
    //    std::cout << "  " << true_evts << "\n";
    //}

    // Fake Data
    //if (syst.fake_data_universe()) {
    //    const Universe* fake_univ = syst.fake_data_universe().get();
    //    for (int b = 0; b < num_bins; ++b) {
    //        double fake_evts = fake_univ->hist_true_->GetBinContent(b+1);
    //        std::cout << "  " << fake_evts << std::endl;
    //    }
    //}


    TH1D* reco_bnb_hist = syst.data_hists_.at( NFT::kOnBNB ).get();      // reconstructed data (real or fake)
    std::cout << "\nData loaded" << std::endl;

    TH1D* reco_ext_hist = syst.data_hists_.at( NFT::kExtBNB ).get();     // reconstructed bkg data
    std::cout << "Off-beam background data loaded" << std::endl;

    std::cout << "-------------------------------" << std::endl;
    std::cout << "Reco (BNB) counts: " << reco_bnb_hist->Integral() << std::endl;
    std::cout << "Reco (EXT) counts: " << reco_ext_hist->Integral() << std::endl;

    std::cout << "\nSubtracting background..." << std::endl;
    syst.data_hists_.at( NFT::kOnBNB ).get()->Add(syst.data_hists_.at( NFT::kExtBNB ).get(),-1);   // new BNB data histogram = BNB - off-beam bkg
    TH1D* after_subtraction = syst.data_hists_.at(NFT::kOnBNB).get();
    std::cout << "Reco (BNB-EXT) counts: " << after_subtraction->Integral() << std::endl;


    TH1D* genie_cv_truth = mcc9->cv_universe().hist_true_.get();        // GENIE CV (truth)
    int num_true_bins = genie_cv_truth->GetNbinsX();
    std::cout << "\nGENIE truth CV loaded" << std::endl;

    TH1D* genie_cv_reco = mcc9->cv_universe().hist_reco_.get();         // GENIE CV (reco)
    int num_reco_bins = genie_cv_reco->GetNbinsX();
    std::cout << "GENIE reco CV loaded\n" << std::endl;


    std::cout << "\nGrabbing covariance matrix...\n" << std::endl;

    auto true_signal = syst.get_cv_true_signal();     // true signal (CV)
    auto reco_signal = syst.get_cv_reco_signal();     // reco signal (CV)
    auto meas = syst.get_measured_events();           // reconstructed events from data (would be good to double check this is background subtracted!)

    const auto& data_signal = meas.reco_signal_;     // reconstructed event count
    const auto& data_covmat = meas.cov_matrix_;      // full covariance matrix (stat + syst)

    // access individual covariance components to get the fractional uncertainties
    auto* matrix_map_ptr = syst.get_covariances().release();
    auto& matrix_map = *matrix_map_ptr;

    
    // Print out event rates
    std::cout << "  [signal + reco bin counts]\n";
    for ( int t = 0; t < num_true_bins; ++t ) { 
        double evts = 0.;
        double error = 0.;
        if ( t < num_true_bins - 1) {
            ///////////////////////////////////////////////////
            std::cout<<"\nCV_true_signal: "<<true_signal->operator()( t, 0 )<<std::endl;
            std::cout<<"CV_reco_signal: "<<reco_signal->operator()( t, 0 )<<std::endl;
            std::cout<<"Background subtracted data [reco_signal]: "<<data_signal->operator()( t, 0 )<<std::endl;
            /////////////////////////////////////////////////         
        }
    }
    // total error (stat + syst)
    std::cout << "\n[Unfolded covariance: sqrt(diagonal elements)]\n";
    for (int i = 0; i < data_covmat->GetNrows(); ++i) {
            std::cout << "  Bin " << i << ": ±" << std::sqrt((*data_covmat)(i,i)) << "\n";   // total error on the reconstructed event rate
    }


    std::cout << "\nFractional uncertainty breakdown:" << std::endl;
    std::cout << "-------------------------------" << std::endl;
    
    double N = data_signal->operator()(0,0);          // reconstructed event count
    double N_err = std::sqrt((*data_covmat)(0,0));    // total error on reconstructed event count
    double N_CV_reco = reco_signal->operator()(0,0);  // CV reco
    double N_CV_true = true_signal->operator()(0,0);  // CV truth


    // fractional uncertainties (for the reconstructed event rate)
    // include variables like this to pull out individual uncertainties
    double sigma_mcstat = 0.0;  // MC stat errors
    for (const auto& pair : matrix_map) {

        const std::string& name = pair.first;
        const auto& cov = pair.second;

        TH2D* covmat = cov.cov_matrix_.get();
        double sigma = std::sqrt( covmat->GetBinContent(1,1) );
        std::cout << name << " sigma = " << sigma << std::endl;
        double frac  = (N > 0.) ? sigma / N : 0.;

        std::cout << name << " frac err = "
                << frac * 100. << "%\n";

        // MC statistical uncertainty
        if (pair.first == "MCstats") {
            TH2D* covmat = pair.second.cov_matrix_.get();
            sigma_mcstat = std::sqrt( covmat->GetBinContent(1,1) );
        }
    }


    // now we can calculate a simple scalar efficiency for our XS calculation
    // No unfolding required for single bin, total cross section (unfolding does give us the same signal counts though which is consistent)
    double eff =  N_CV_reco / N_CV_true;   // use reco CV / truth CV (ordinarily this is taken from the smearing)
    std::cout << "\nefficiency (CV reco / CV true) = (" << N_CV_reco << "/" << N_CV_true << ") = " << eff << std::endl;


    // XS equation: (N-B) / [eff * N_targets * Phi]
    // conversion factor = N_targets * Phi [cm^2 / nucleon], then divide by 10^-38

    // calculate cross section from data (already background subtracted)
    double val_xsec = N / (eff * conv_factor);
    double err_xsec = N_err / (eff * conv_factor);

    // calculate CV cross section for closure / validation
    // "Does the measurement machinery reproduce the known truth within its own uncertainties?"
    double val_xsec_CV = N_CV_reco / (eff * conv_factor);
    double err_xsec_CV_MCstat = sigma_mcstat / (eff * conv_factor);    // MC stat only

    std::cout << "\nXS calculated -----------------" << std::endl;
    std::cout << "Data cross section = " << val_xsec << " ± " << err_xsec << " cm^2 / oxygen\n" << std::endl;
    std::cout << "GENIE CV cross section = " << val_xsec_CV << " ± " << err_xsec_CV_MCstat << " (MC stat only) cm^2 / oxygen\n" << std::endl;


    //
    // ******************************************
    // Plot
    //
    std::cout << "\n\nPlotting -----------------\n" << std::endl;

    // ExtBNB vs onBNB
    // reconstructed event rate plot (reco CV vs reco data)
    // cross section plot

    // define bin edges for plotting histograms
    //std::vector<double> nbins = {0, 1};   // bin edges

    // PDF canvas
    TCanvas* c4 = new TCanvas("c4");
	  c4->Print(text_title_pdf); 


    // .........................
    // ExtBNB vs onBNB
    Draw_HIST(
        reco_bnb_hist,
        "Fake data (GENIE CV)",
        reco_ext_hist,
        "Off-Beam background",
        "",
        "",
        "Events",
        false,  // norm area
        true,   // set grid
        false,  // bin width norm
        -99,
        c4,
        text_title_pdf_string);


    // .........................
    // ExtBNB vs onBNB
    Draw_HIST(
        reco_bnb_hist,
        "Fake data (GENIE CV)",
        genie_cv_reco,
        "GENIE CV Prection",
        "",
        "",
        "Events",
        false,  // norm area
        true,   // set grid
        false,  // bin width norm
        -99,
        c4,
        text_title_pdf_string);


    // .........................
    // Cross section

    TH1D* h_xsec_data = new TH1D(
        "h_xsec_data",
        "; ;#sigma_{#nu NCQE} [10^{-38} cm^{2} / oxygen]",
        1, 0.5, 1.5
    );

    TH1D* h_xsec_cv = new TH1D(
        "h_xsec_cv",
        "; ;#sigma_{#nu NCQE} [10^{-38} cm^{2} / oxygen]",
        1, 0.5, 1.5
    );

    // Data
    h_xsec_data->SetBinContent(1, val_xsec);
    h_xsec_data->SetBinError(1, err_xsec);
    h_xsec_data->SetMarkerStyle(20);
    h_xsec_data->SetMarkerSize(1.4);
    h_xsec_data->SetLineColor(kBlack);
    h_xsec_data->SetLineWidth(2);
    h_xsec_data->SetStats(0);
    h_xsec_data->GetXaxis()->SetTickLength(0);
    h_xsec_data->GetXaxis()->SetLabelSize(0);


    // GENIE CV band (MC stat only)
    double x_min = 0.5;
    double x_max = 1.5;
    double y_min = val_xsec_CV - err_xsec_CV_MCstat;
    double y_max = val_xsec_CV + err_xsec_CV_MCstat;

    TLine* genie_cv_line = new TLine(
        x_min,
        val_xsec_CV,
        x_max,
        val_xsec_CV
    );
    genie_cv_line->SetLineColor(kAzure);
    genie_cv_line->SetLineWidth(2);
    genie_cv_line->SetLineStyle(2);  // dashed (optional)

    // +/- MC stat error band
    TLine* genie_cv_line_up = new TLine(
        x_min, val_xsec_CV + err_xsec_CV_MCstat,
        x_max, val_xsec_CV + err_xsec_CV_MCstat
    );

    TLine* genie_cv_line_dn = new TLine(
        x_min, val_xsec_CV - err_xsec_CV_MCstat,
        x_max, val_xsec_CV - err_xsec_CV_MCstat
    );

    genie_cv_line_up->SetLineColor(kAzure);
    genie_cv_line_dn->SetLineColor(kAzure);

    genie_cv_line_up->SetLineWidth(2);
    genie_cv_line_dn->SetLineWidth(2);

    genie_cv_line_up->SetLineStyle(2);
    genie_cv_line_dn->SetLineStyle(2);

    TCanvas* c_xs = new TCanvas("c_xs", "NCQE Cross Section", 600, 500);
    h_xsec_data->GetYaxis()->SetTitleOffset(1.3);
    h_xsec_data->GetXaxis()->SetLabelSize(0);  // hide x labels (single bin)

    double y_min_plot = 0.5;
    double y_max_plot = 7.0;
    h_xsec_data->SetMinimum(y_min_plot);
    h_xsec_data->SetMaximum(y_max_plot);

    h_xsec_data->Draw("E1");
    //h_xsec_cv->Draw("E1 SAME");
    genie_cv_line_up->Draw("SAME");
    genie_cv_line_dn->Draw("SAME");
    h_xsec_data->Draw("E1 SAME");  // redraw point on top
    genie_cv_line->Draw("SAME");

    TLegend* leg = new TLegend(0.55, 0.70, 0.85, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(h_xsec_data, "Fake BNB Data (stat. #oplus syst.)", "lep");
    //leg->AddEntry(h_xsec_cv,   "GENIE CV (MC stat only)",     "lep");
    leg->AddEntry(genie_cv_line, "GENIE CV Truth (MC stat. only)",     "l");

    leg->Draw();
    c_xs->SaveAs(text_title_pdf);



    std::cout<<"\nSaving PDF...\n"<<std::endl; 

  	sprintf(text_title_pdf, "NC_univmake_plots.pdf)");
  	c4->Print(text_title_pdf);
    std::cout<<"\n"; 
}


// main function
int main(int argc, char* argv[]) {
   Test_plot();
   return 0;
}
