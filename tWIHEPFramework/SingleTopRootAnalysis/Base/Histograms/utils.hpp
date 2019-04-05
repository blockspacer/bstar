/******************************************************************************
 * utils.hpp                                                                  *
 *                                                                            *
 * Helper functions which could be useful                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Functions                                                                  *
 *    Double_t etaPhysics2Detector()    -- Detector eta -> Physics eta        *
 *    Double_t etaDetector2Physics()    -- Detector eta -> Physics eta        *
 *    Double_t Neutrino_Pz()            -- Get neutrino pz from W mass        *
 *    Double_t Neutrino_Pz_Top()        -- Get neut pz from W with top=178Gev *
 *    Double_t Neutrino_Pz_Base()       -- Assume lepton & nuetrino massless  *
 *    Int_t FindBestJet()               -- Find jet giving top = 178 GeV      *
 *    Int_t FindBestJet()               -- Find jet giving top = 178 GeV      *
 *                                                                            *
 * History                                                                    *
 *       9 Sep 2004 - Helper functions                                        *
 *      17 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/
#ifndef utils_h
#define utils_h

#include <TLorentzVector.h>
#include <math.h>
#include <vector>
#include "SingleTopRootAnalysis/Particles/Recon/Jet.hpp"

// Convert Physics eta into Detector eta
Double_t etaPhysics2Detector(const Double_t& etain, const Double_t& zvertex, const Double_t& R = 72.0);

// Convert Detector eta into Physics eta
Double_t etaDetector2Physics(const Double_t& etain, const Double_t& zvertex, const Double_t& R = 72.0);

// Find the neutrino pz
Double_t Neutrino_Pz( Double_t nu_px, Double_t nu_py, 
		      Double_t l_px, Double_t l_py, Double_t l_pz, Double_t l_e, 
		      Double_t& scf, Double_t nu_pz_true = 0.0);

// Find the neutrino pz with top mass constraint
Double_t Neutrino_Pz_Top( Double_t nu_px, Double_t nu_py, 
			  Double_t l_px, Double_t l_py, Double_t l_pz, Double_t l_e, 
			  Double_t b_px, Double_t b_py, Double_t b_pz, Double_t b_e, 
			  Double_t nu_pz_true = 0.0);

// Find the neutrino pz with neutrino and lepton mass set to 0
void Neutrino_Pz_Base( Double_t nu_px, Double_t nu_py,                            // neutrino 
		       Double_t l_px, Double_t l_py, Double_t l_pz, Double_t l_e, // lepton
		       Double_t& scf,Double_t &solution1, Double_t &solution2 );   // return variables


// Find the "best jet", which is the one that gives a top mass closest to 175GeV.
Int_t findBestJet(Int_t njets, const std::vector<Jet> jetVector, const TLorentzVector W, Double_t targetmass = 175.0);

// Find the "best jet" which is the one that gives a top mass closest to 175GeV and the ideal W.
// That is, combine the top mass constraint and the W mass constraint
Int_t findBestJetW(Int_t njets, TLorentzVector* jets, TLorentzVector ele, TLorentzVector nu, Double_t targetmass = 175.0);
#endif
