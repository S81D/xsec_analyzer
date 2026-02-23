// Header of geometry tool 
//// Hist tools 
/// Christian Nguyen 
#pragma once

#include <algorithm>
#include <cmath>

// ROOT includes
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "THStack.h"
#include "TLegend.h"
#include "TObjArray.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TGraph2D.h>

// STV analysis includes

//#include "SliceBinning.hh"
//#include "SliceHistogram.hh"
#include "TLatex.h"
#include "TLine.h"
#include "TH2Poly.h"
#include "UBTH2Poly.h" 
#include "GridCanvas.hh"
 
//#include "ConfigMakerUtils.hh"
//#include "../EventCategory.hh"
#include "NamedCategory.hh"
#include "HistUtils.h"
#include "TVector3.h"


TGraph* createBoundaryGraph(const std::string& plotType);


//extern double FV_X_MIN =   21.5;
//extern double FV_X_MAX =  234.85;
//
//extern double FV_Y_MIN = -95.0;
//extern double FV_Y_MAX =  95.0;
//
//extern double FV_Z_MIN =   21.5;
//extern double FV_Z_MAX =  966.8;
//
////ANNIE FV Boundaries
//extern double FV_RAD = 100.0;
//
//extern double A_FV_Y_MIN = -100.0;
//extern double A_FV_Y_MAX =  100.0;
//
//extern double A_Z_CTR = 168.1; //Z center of tank


extern double FV_X_MIN;
extern double FV_X_MAX;

extern double FV_Y_MIN;
extern double FV_Y_MAX;

extern double FV_Z_MIN;
extern double FV_Z_MAX;

//ANNIE FV Boundaries
extern double FV_RAD;

extern double A_FV_Y_MIN;
extern double A_FV_Y_MAX;

extern double A_Z_CTR; //Z center of tank





// Use a template here so that this function can take float or double values as
// input
template <typename Number> bool point_inside_FV( Number x, Number y, Number z )
{
  bool x_inside_FV = ( FV_X_MIN < x ) && ( x < FV_X_MAX );
  bool y_inside_FV = ( FV_Y_MIN < y ) && ( y < FV_Y_MAX );
  bool z_inside_FV = ( FV_Z_MIN < z ) && ( z < FV_Z_MAX );
  return ( x_inside_FV && y_inside_FV && z_inside_FV );
}

inline bool point_inside_FV( const TVector3& pos ) {
  return point_inside_FV( pos.X(), pos.Y(), pos.Z() );
}

template <typename Number> bool point_inside_AFV( Number x, Number y, Number z )
{
  if(y > A_FV_Y_MIN && y < A_FV_Y_MAX && FV_RAD > std::sqrt((z - A_Z_CTR)*(z - A_Z_CTR) + x*x) && (z - A_Z_CTR < FV_RAD)){
    return true;
  }
  else {
    return false;
  }
}

inline bool point_inside_AFV( const TVector3& pos ) {
  return point_inside_AFV( pos.X(), pos.Y(), pos.Z() );
}

// Returns the number of Ar nuclei inside the fiducial volume
inline double num_Ar_targets_in_FV() {
  double volume = ( FV_X_MAX - FV_X_MIN ) * ( FV_Y_MAX - FV_Y_MIN )
    * ( FV_Z_MAX - FV_Z_MIN ); // cm^3
  constexpr double m_mol_Ar = 39.948; // g/mol
  constexpr double N_Avogadro = 6.02214076e23; // mol^(-1)
  constexpr double mass_density_LAr = 1.3836; // g/cm^3

  double num_Ar = volume * mass_density_LAr * N_Avogadro / m_mol_Ar;
  return num_Ar;
}

// Returns the number of O nuclei inside the fiducial volume
inline double num_O_targets_in_FV() {
  double volume = FV_RAD * FV_RAD * M_PI * ( A_FV_Y_MAX - A_FV_Y_MIN ); // cm^3
  constexpr double m_mol_H2O = 18.015; // g/mol
  constexpr double N_Avogadro = 6.02214076e23; // mol^(-1)
  constexpr double mass_density_H2O = 1.; // g/cm^3

  // if you are assuming pure water
  //double num_O = volume * mass_density_H2O * N_Avogadro / m_mol_H2O;

  // if considering Gd-loading (its literally a 0.1% difference --> should work out to double num_O = 3.147e29;)
  double gd_fraction = 0.001;  // mass fraction of Gd (0.001 = 0.1% for Phase II)
  double mass_H2O = volume * mass_density_H2O * (1.0 - gd_fraction);  // g (mass of water within the FV considering Gd loading)
  double num_O = mass_H20 * N_Avogadro / m_mol_H2O;  // number of O16 targets given Gd loading

  //double num_O = 3.147e29; // accounting for Gd loading (< 1% difference) 
 
  return num_O;
}

// Returns the total BNB neutrino flux (nu / cm^2) in the fiducial
// volume as a function of a given beam exposure (measured in
// protons-on-target)
// NOTE: This is currently approximated using the flux in the *active volume*.
// NOTE: this has been adapted from MicroBooNE + ANNIE CC to give you TOTAL nu instead of just numu
inline double integrated_numu_flux_in_FV( double pot ) {
  //
  // for NC, we must include both numu and nue --> the value below is the combined flux
  constexpr double numu_per_cm2_per_POT_in_AV = 2.08722e-8;  // 2.26256e-8 is what James uses (just nu mu, in his condensed FV)
                                                             // derived values for me (r < 1m, |y| < 1.5m:          
                                                             // Integrated νμ flux: 2.028605343287099e-08 per POT per cm^2
                                                             // Integrated νe flux: 1.120153603922698e-10 per POT per cm^2
                                                             // updated FV (r < 0.7m, |y| < 1.0m:          
                                                             // Integrated νμ flux: 2.075647937553975e-08 per POT per cm^2
                                                             // Integrated νe flux: 1.1571293632036711e-10 per POT per cm^2
  double flux = pot * numu_per_cm2_per_POT_in_AV; // numu / cm^2
  return flux;
}
