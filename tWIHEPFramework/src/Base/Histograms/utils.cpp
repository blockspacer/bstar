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

#include "TMath.h"
#include "SingleTopRootAnalysis/Base/Histograms/utils.hpp"
#include <string>
#include <iostream>

using namespace std;

/******************************************************************************
 * Double_t etaDetector2Physics()                                             *
 *                                                                            *
 * Convert Detector eta into Physics eta                                      *
 * R is the distance from beamline in cm(~72 is CPS)                          *
 * This simply calls the method below with -zvertex                           *      
 *                                                                            *
 * Input:  None                                                               *
 * Output: Physics Eta                                                        *
 ******************************************************************************/
Double_t etaDetector2Physics(const Double_t& etain,const Double_t& zvertex,const Double_t& R)
{
  return(etaPhysics2Detector(etain,-zvertex,R));
} //etaDetector2Physics()

/******************************************************************************
 * Double_t etaPhysics2Detector()                                             *
 *                                                                            *
 * Convert Physics eta into Detector eta                                      *
 * R is the distance from beamline in cm(~72 is CPS)                          *
 *                                                                            *
 * Input:  None                                                               *
 * Output: Physics Eta                                                        *
 ******************************************************************************/
Double_t etaPhysics2Detector(const Double_t& physEta, const Double_t& zvertex,const Double_t& R) {

  Double_t detEta     = 0.0;
  Double_t Theta      = 0.0;
  Double_t ThetaPrime = 0.0;
  Double_t d          = 0.0;

  //Do case where physEta is positive
  if (physEta >= 0.0) {
    //physEta is positive, zvertex is positive
    if (zvertex >= 0.0) {
      Theta = 0.0;
      ThetaPrime = 0.0;
      Theta = 2.0 * TMath::ATan(exp(-1.0 * physEta));
      ThetaPrime = TMath::ATan(R / (zvertex + R / TMath::Tan(Theta)));
      detEta = -1.0 * TMath::Log(TMath::Tan(ThetaPrime/2.0));
    } //if
    //physEta is positive, zvertex is negative
    else {
      Theta = 0.;
      ThetaPrime = 0.;
      Theta = 2.0 * TMath::ATan(exp(-1.0 * physEta));
      d = R / TMath::Tan(Theta);
      if (d > (-1.0 * zvertex)) {
        ThetaPrime = TMath::ATan(R/(zvertex+d));
        detEta = -1.*TMath::Log(TMath::Tan(ThetaPrime/2.));
      } //if
      //case where positive physEta becomes negative detector Eta
      else {
        ThetaPrime = TMath::ATan(R / (-1.0 * (zvertex+d)));
        detEta = TMath::Log(TMath::Tan(ThetaPrime / 2.0));
      } //else
    } //else
  } //if

  //Do case where physEta is negative
  else {
    //physEta is negative, zvertex is positive
    if (zvertex >= 0) {
      Theta = 0.0;
      ThetaPrime = 0.0;
      Theta = 2.0 * TMath::ATan(exp(physEta));
      d = R / TMath::Tan(Theta);
      if (d > zvertex) {
        ThetaPrime = TMath::ATan(R / (d-zvertex));
        detEta = TMath::Log(TMath::Tan(ThetaPrime/2.));
      } //if
      // Case where negative physEta becomes positive detector Eta
      else {
        ThetaPrime = TMath::ATan(R / (zvertex-d));
        detEta = -1.0 * TMath::Log(TMath::Tan(ThetaPrime/2.0));
      } //else
    } //if

    //physEta is negative, zvertex is negative
    else {
      Theta = 0.0;
      ThetaPrime = 0.0;
      Theta = 2.0 * TMath::ATan(exp(physEta));
      ThetaPrime = TMath::ATan(R / (-1.0 * zvertex + R/TMath::Tan(Theta)));
      detEta = TMath::Log(TMath::Tan(ThetaPrime / 2.0));
    } //else
  } //else
  return detEta;
} //etaPhysics2Detector

/******************************************************************************
 * Double_t Neutrino_Pz()                                                     *
 *                                                                            *
 * Get neutrino Pz from the W mass constraint.                                *
 * Pick the solution with the smaller pz                                      *
 *                                                                            *
 * This function is from the D0 top_tree package.                             *
 * Modification: only do scaling if scf is passed in (if it is not -99999.)   *
 * Note: if input scf is set to a positive num then a scale factor is applied *
 *  This means we only apply and return the scale factor if scf is >=0.       *
 *                                                                            *
 * Input:  Lepton and neutrino info                                           *
 * Output: neutrino pz                                                        *
 ******************************************************************************/
Double_t Neutrino_Pz( Double_t nu_px, Double_t nu_py, 
		    Double_t l_px, Double_t l_py, Double_t l_pz, Double_t l_e, 
		    Double_t& scf,Double_t nu_pz_true)
{
  Double_t solution1;
  Double_t solution2;

  // Find the two possible solutions:
  Neutrino_Pz_Base( nu_px, nu_py, l_px, l_py, l_pz, l_e, // input arguments 
		    scf, solution1, solution2);          // return variables

  // Now choose the lower one of the two
  //if(nu_pz_true>100)
  //cout<<"possible solutions: "<<solution1<<" and "<<solution2;
  //cout<<" and truth is "<<nu_pz_true<<" and pz_l "<<l_pz<<endl;
  // choose solution with smallest |l_pz| a la Run I
  // alternatives don't work !!
  Double_t nu_pz = fabs (solution1) < fabs (solution2) ? solution1 : solution2 ;

    return nu_pz;
  
} //Neutrino_Pz
 
/******************************************************************************
 * Double_t Neutrino_Pz_Top()                                                 *
 *                                                                            *
 * Get neutrino Pz from the W mass constraint.                                *
 * Pick the solution giving a top mass closest to 178 Gev                     *
 *                                                                            *
 * Input:  Lepton and neutrino info                                           *
 * Output: neutrino pz                                                        *
 ******************************************************************************/
Double_t Neutrino_Pz_Top( Double_t nu_px, Double_t nu_py, 
			Double_t l_px, Double_t l_py, Double_t l_pz, Double_t l_e, 
			Double_t b_px, Double_t b_py, Double_t b_pz, Double_t b_e, 
			Double_t nu_pz_true)
{
  Double_t solution1;
  Double_t solution2;
  Double_t scf = -99999.0;

  // Find the two possible solutions:
  Neutrino_Pz_Base( nu_px, nu_py, l_px, l_py, l_pz, l_e,  // input arguments 
		    scf, solution1, solution2);           // return variables

  // Now choose a solution
  //cout<<"  possible solutions: "<<solution1<<" and "<<solution2;
  //cout<<" and truth is "<<nu_pz_true<<endl;
  
  // compute the two solutions to the top mass:

  // compute two W invariant masses:
  Double_t nu_e1 = sqrt(pow(nu_px,2)+pow(nu_py,2)+pow(solution1,2));
  //Double_t w_m1=sqrt( pow(nu_e1+l_e,2)-
  //	    (pow(nu_px+l_px,2)+pow(nu_py+l_py,2)+pow(solution1+l_pz,2)));
  Double_t nu_e2 = sqrt(pow(nu_px,2)+pow(nu_py,2)+pow(solution2,2));
  //Double_t w_m2=sqrt( pow(nu_e2+l_e,2)-
  //	    (pow(nu_px+l_px,2)+pow(nu_py+l_py,2)+pow(solution2+l_pz,2)));

  Double_t bl_px   = l_px  + b_px;
  Double_t bl_py   = l_py  + b_py;
  Double_t bl_pz   = l_pz  + b_pz;
  Double_t bl_e    = l_e   + b_e;
  Double_t top_e1  = bl_e  + nu_e1;
  Double_t top_e2  = bl_e  + nu_e2;
  Double_t top_px  = bl_px + nu_px;
  Double_t top_py  = bl_py + nu_py;
  Double_t top_pz1 = bl_pz + solution1;
  Double_t top_m1  = sqrt(pow(top_e1,2) - pow(top_px,2) - pow(top_py,2) - pow(top_pz1,2));
  Double_t top_pz2 = bl_pz + solution2;
  Double_t top_m2  = sqrt(pow(top_e2,2) - pow(top_px,2) - pow(top_py,2) - pow(top_pz2,2));
  //cout<<" reco Top mass 1 is "<<top_m1<<", e1= "<<top_e1<<", pz1= "<<top_pz1<<endl;
  //cout<<" reco Top mass 2 is "<<top_m2<<", e2= "<<top_e2<<", pz2= "<<top_pz2<<endl;
  
  Double_t nu_pz = fabs(top_m1 - 178.0) < fabs(top_m2 - 178.0) ? solution1 : solution2;
  
  return nu_pz;
  
} //Neutrino_Pz_Top
 
/******************************************************************************
 * void Neutrino_Pz_Base()                                                    *
 *                                                                            *
 * Get neutrino Pz from the W mass constraint.                                *
 * Assume the lepton and neutrino are massless                                *
 *                                                                            *
 * Input:  Lepton and neutrino info                                           *
 * Output: neutrino pz                                                        *
 ******************************************************************************/
void Neutrino_Pz_Base( Double_t nu_px, Double_t nu_py, 
		       Double_t l_px, Double_t l_py, Double_t l_pz, Double_t l_e,
		       Double_t& scf,Double_t &solution1, Double_t &solution2)
{
  Bool_t usescf = kTRUE;
  if(scf < 0.0) {
    usescf = kFALSE;
    //cout<<"not using scale factor!"<<endl;
  } //if

  scf = 1.0;
  Double_t nu_pz = 0.;
  Double_t l_pt  = TMath::Sqrt(pow(l_px,2) + pow(l_py,2));
  
  Double_t Mw    = 80.42;
  Double_t Mt;
  Double_t nu_pt = TMath::Sqrt(pow(nu_px,2)+pow(nu_py,2));
  Mt = TMath::Sqrt(pow(l_pt+nu_pt,2) - pow(l_px+nu_px,2) - pow(l_py+nu_py,2));
  
  //  cout<<"   W mass is "<<Mw<<", mt is "<<Mt<<endl;
  Double_t A;
  //rs: only change the scale factor if the user wants to.
  if (Mt < Mw || !usescf) A = pow(Mw,2) / 2.;
  // assume Mt=Mw, rescale MET accordingly
  else       {        
    cout<<"  scaling met"<<endl;
    A = pow(Mt,2) / 2.;
    Double_t k = nu_pt * l_pt - nu_px * l_px - nu_py * l_py;
    k = (k == 0 ? 0.00001 : k);
    scf = 0.5 * pow(Mw,2) / k;
    nu_px *= scf;
    nu_py *= scf;
    nu_pt = sqrt(pow(nu_px,2) + pow(nu_py,2));
  } //else
  
  // If the electron goes along the beam direction (Px==Py==0),
  // then the algorithm below doesn't work, so do the calculation here.
  Double_t S_denom = pow(l_pz,2) - pow(l_e,2);
  if(S_denom == 0.0) {
    //cout<<"lepton along beam direction, Mt is "<<Mt<<endl;
    //cout<<"lpx "<<l_px<<", l_py "<<l_py<<", l_pz "<<l_pz<<", n_px "<<nu_px<<", n_py "<<nu_py<<endl;
    // in this case there is one exact solution
    Double_t mwplz = pow(Mw,2) / fabs(l_pz);
    nu_pz= (pow(nu_px,2) + pow(nu_py,2) - pow(mwplz/2.,2)) / mwplz;
    solution1 = nu_pz;
    solution2 = nu_pz;
    return;
  } //if
  
  Double_t B = nu_px * l_px + nu_py * l_py;
  Double_t C;
  Double_t C1 = 1.0E0 + pow(nu_pt,2) * (pow(l_pz,2)-pow(l_e,2)) / pow(A+B,2);
  if(C1 <= 0.0) {
    //cout<<"Warning in Neutrino_Pz! negative C1 "<<C1<<", Setting it to 0."<<endl;
    // in this case C will be zero and there is only one solution, 
    // so simply return the one possible solution
    solution1 = (-(A+B) * l_pz) / (pow(l_pz,2) - pow(l_e,2));
    solution2 = solution1;
    return;
  } //if
  else {
    C = TMath::Sqrt(C1);
    solution1 = (-(A+B)*l_pz + (A+B)*l_e*C) / S_denom;
    solution2 = (-(A+B)*l_pz - (A+B)*l_e*C) / S_denom;
    return;
  } //else
  
} //Neutrino_Pz_Base

//
// find the "best jet", which is the one that gives a top mass
// closest to 175GeV.
// return the index in the array that points to the best jet.
/******************************************************************************
 * void findBestJet()                                                         *
 *                                                                            *
 * Find best jet                                                              *
 * Best jet gives a top mass closest to 175 GeV                               *
 *                                                                            *
 * Input:  Jet and W info, target top mass                                    *
 * Output: neutrino pz                                                        *
 ******************************************************************************/
Int_t findBestJet(Int_t njets, const std::vector<Jet> jetVector, const TLorentzVector W, Double_t targetMass)
{
  // Loop over the jets 
  // For each iteration reconstruct the invariant mass of the jet-W system.
  TLorentzVector top;
  Double_t mt, diff, mindiff = 999999999.0;
  Double_t mtarg = pow(targetMass,2);

  Int_t foundJet = 0;
  //cout<<"in findBestJet there are njets "<<njets<<endl;
  for(Int_t i=0;i<njets;i++) {
    top = W + jetVector[i];
    // for some reason M() returns a negative mass every other time,
    // so we have to make sure we take the absolute value
    mt = top.Mag2();
    diff = fabs(mt-mtarg);
    //cout<<"  jet "<<i<<", mass "<<top.M()<<", diff "<<diff<<endl;
    if(diff < mindiff) {
      mindiff = diff;
      foundJet = i;
      //cout<<"  best mass for jet "<<i<<" is "<<top.M()<<endl;
    } //if
  } // end of loop over jets

  return foundJet;
} //findBestJet
