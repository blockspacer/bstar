#ifndef MET_Cleaning_Utils_h
#define MET_Cleaning_Utils_h

#include <vector>
#include "TLorentzVector.h"
#include "SingleTopRootAnalysis/Base/Dictionary/Top_MET.hpp"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// class MET_Cleaning_Utils
// Proves a standard interface to modify the missing ET
// -enables modification fo the missing ET for default corrections to the selected objects and for studies of systematic uncertainties. 
// -implements the MET only uncertainties for the CellOut term and the SoftJet term.   
//
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class MET_Cleaning_Utils {

 public:
  MET_Cleaning_Utils();
  virtual ~MET_Cleaning_Utils();

  // Common funtions
  float DeltaPhi(float phi1, float phi2);
  float DeltaEta(float eta1, float eta2);
  float dr(float phi1, float eta1,float phi2, float eta2);

  //Cleaning Tools
  bool isJetSafe(std::vector<TLorentzVector> &jet, std::vector<float> *BCH_CORR_JET,
		 std::vector<float> *BCH_CORR_CELL, bool isdata, float shift);
  bool isJetSafe(std::vector<TLorentzVector> &jet, std::vector<float> *BCH_CORR_JET,
                 std::vector<float> *BCH_CORR_CELL, bool isdata, float shift, int RunNumber);
  bool HasLArFailure(int RunNumber);

  bool hasUnmatchTrackJet(std::vector<TLorentzVector> &jet,float JCone,std::vector<TLorentzVector> &trackjet,
			  float TJCone,float minpt);
  bool isJetAffected(std::vector<TLorentzVector> &jet, float iso, float JCone, float minpt );
  bool isJetAffected(TLorentzVector &jet, float iso, float JCone, float minpt );
  bool iswithinDeadFEB(float eta, float phi, float objCone, float objIso,bool dosquare=false);


  // Validation Tools
  float dphi_MET_MPT(float &met_phi,std::vector<float> *trk_pt,std::vector<float> *trk_eta,
		     std::vector<float> *trk_d0_wrtPV,std::vector<float> *trk_z0_wrtPV,
		     std::vector<float> *trk_theta_wrtPV,std::vector<float> *trk_phi_wrtPV,
		     std::vector<int> *trk_nPixHits,std::vector<int> *trk_nSCTHits);
  TLorentzVector MPT_Loose(std::vector<float> *trk_pt,std::vector<float> *trk_eta,std::vector<float> *trk_d0_wrtPV,
			   std::vector<float> *trk_z0_wrtPV,std::vector<float> *trk_theta_wrtPV,
			   std::vector<float> *trk_phi_wrtPV,std::vector<int> *trk_nPixHits,
			   std::vector<int> *trk_nSCTHits );
  float dphi_Jet_MET(std::vector<TLorentzVector> &jet, float met_phi);
  float JetTrackJetRatio(std::vector<TLorentzVector> &jet,TLorentzVector &trackjet);

  float JetTrackJetphi(TLorentzVector &jet,std::vector<TLorentzVector> &trackjet,
		       std::vector<TLorentzVector> &jettruth,float smear, bool isMC);
  float dphi_MET_MPT(float met_phi,float met_truthphi, float smear, bool isMC,std::vector<float> *trk_pt,std::vector<float> *trk_eta,
                     std::vector<float> *trk_d0_wrtPV,std::vector<float> *trk_z0_wrtPV,
                     std::vector<float> *trk_theta_wrtPV,std::vector<float> *trk_phi_wrtPV,
                     std::vector<int> *trk_nPixHits,std::vector<int> *trk_nSCTHits);
  int Jetmatching(TLorentzVector &jet,std::vector<TLorentzVector> &jetColl2,float Cone);

  
 private:
  std::vector<float> eta_min, eta_max, phi_min, phi_max;
  float GeV;
  float jet_THRESHOLD;
  float jet_THRESHOLD_1sigma;
  float jet_iso;
  Top_MET top_met;
};
#endif
