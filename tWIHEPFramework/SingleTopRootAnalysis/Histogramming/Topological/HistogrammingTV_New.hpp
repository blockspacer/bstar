/******************************************************************************
 * HistogrammingTV_New.hpp                                                  *
 *                                                                            *
 * Fills a tree with specific variables.                                      *
 * Should be called after preselection cuts in bin/analysis                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTV_New()              -- Parameterized Constructor       *
 *    ~HistogrammingTV_New()             -- Destructor                      *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTV_New"    *
 *                                                                            *
 * Private Data Members of this class                                         *
 *        TFile *_skimFile                                                    *
 *        TTree *_skimTree;                                                   *
          Lots of variables, histograms                                       *
 *                                                                            *
 * History                                                                    *
 *  arguments: - pointer to input particles class                             *
 *             - write: whether to write out the NN skim or not               *
 *             - name: name of the NNvars skim file                           *
 *             - doSpin: whether to include top spin correlation vars in NN   *                                                       
 *                                                                            *
 * History                                                                    *
 *      6 Nov 2009 - P. Ryan - added angular variables and cleaned up         *
 *****************************************************************************/

#ifndef HistogrammingTV_New_h
#define HistogrammingTV_New_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/ScaleWjets.hpp"

#include "SingleTopRootAnalysis/Base/Dictionary/MuidCBScaleEffFactors.hpp"

#include "SingleTopRootAnalysis/Base/Dictionary/MuonSF.hpp"

//#include "config/scalefactors/muon_SF_pLHC.h"
//#include "config/scalefactors/muon_SF_EPS.h"
//#include "config/scalefactors/electron_SF_pLHC_sgtop.h"
#include "config/scalefactors/electron_SF_EPS_sgtop.h"
#include "config/scalefactors/WJetSF_EPS.h"

#include <string>
#include <vector>
#include <sstream>

#include <map>
#include <cmath>

//#include "CLHEP/Vector/LorentzVector.h"

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVectorD.h"


using namespace std;


class HistogrammingTV_New : public HistoCut 
{

public:

  HistogrammingTV_New(EventContainer *obj,TFile *skimFilez);
  
  // destructor
  ~HistogrammingTV_New();

  // methods:
  void BookHistogram();
  bool Apply();
  
  inline std::string GetCutName (void) { return "HistogrammingTV_New"; };

 private:

  TFile *skimFile;
  TTree *skimTree;
  std::vector<Double_t> bweightvec;
  ScaleWjets* scpt10;
  ScaleWjets* scpt20;
  ScaleWjets* scopt2;
  ScaleWjets* scopt3;
  TLorentzVector *jet;
  MuidCBScaleEffFactors MSF;
  MuonSF MUSFTrigID; 
 
  Float_t actualIntPerXing;
  Float_t averageIntPerXing;

  Float_t jetweight_opt2;
  Float_t jetweight_opt3;
  Float_t jetweight_pt10;
  Float_t jetweight_pt20;
  // skim tree variables
  Float_t mulbn;
  Float_t P_AllJetsSum;
 Float_t distns;
  Float_t distbunch;
  Float_t bcid;
  Int_t Pvtxall_n;

  Int_t safejetevent;
  Int_t safejeteventup;
  Int_t safejeteventdown;

  //------------------------------------------------------
  // Event Weight
  //------------------------------------------------------
  Float_t EventWeight;
  Float_t TotalWeight;

  Float_t MuonSFTrigger;
  Float_t MuonSFID;
  Float_t MuonSFRecoIDTrigger;
  Float_t ElectronSFTrigger;
  Float_t ElectronSFRecoID;
  Float_t ElectronSFRecoIDTrigger;
  Float_t WJetSF;

  Float_t  MuonSFTriggerErr;
  Float_t  MuonSFTriggerErrdown;
  Float_t  MuonSFIDErr;
  
  Float_t  ElectronSFTriggerErr;
  Float_t  ElectronSFRecoIDErr;
  
  Float_t WJetSFErrTtbarDown;
  Float_t WJetSFErrTtbarUp;
  Float_t WJetSFErrSgTopDown;
  Float_t WJetSFErrSgTopUp;
  Float_t WJetSFErrZjetsDown;
  Float_t WJetSFErrZjetsUp;
  Float_t WJetSFErrDibDown;
  Float_t WJetSFErrDibUp;
  Float_t WJetSFErrQCDDown;
  Float_t WJetSFErrQCDUp;
  Float_t WJetSFErrJESDown;
  Float_t WJetSFErrJESUp;
  Float_t WJetSFErrBtagSFDown;
  Float_t WJetSFErrBtagSFUp;
  Float_t WJetSFErrLqtagSFDown;
  Float_t WJetSFErrLqtagSFUp;
  
  Float_t WJetSFErrSysterrPercent;
  Float_t WJetSFErrStaterrPercent;
  Float_t WJetSFErrTotalerrPercent;


  Int_t HFOR;
  Int_t Pvtx_n;
  Int_t pdgid;
 
  Int_t NWDecayLeptons;
  Int_t NWDecayLightQuarks;

  Float_t PileupWeight; //tagging weight per event
  Float_t TaggingWeight; //tagging weight per event
  Float_t TaggingWeight_Bup;//tagging weight per event applying the b SF up shift
  Float_t TaggingWeight_Bdown;//tagging weight per event applying the b SF down shift
  Float_t TaggingWeight_Lqup;//tagging weight per event applying the Lq SF up shift
  Float_t TaggingWeight_Lqdown;//tagging weight per event applying the Lq SF down shift

  Float_t ElectronDetEta;
  Float_t QTimesEta;

  Double_t Cos_Lepton_Jet2_RJet1Top;

  Float_t JetProbWeight1;
  Float_t JetProbWeight2;
  Float_t JetProbWeight3;
  Float_t JetProbWeight4;

  Float_t SV0Weight1;
  Float_t SV0Weight2;
  Float_t SV0Weight3;
  Float_t SV0Weight4;


  Float_t JetFitterCombNNWeight1;
  Float_t JetFitterCombNNWeight2;
  Float_t JetFitterCombNNWeight3;
  Float_t JetFitterCombNNWeight4;


  Float_t IP3DSV1Weight1;
  Float_t IP3DSV1Weight2;
  Float_t IP3DSV1Weight3;
  Float_t IP3DSV1Weight4;

  Int_t JetTruthLabel;

  Float_t         Lepton_cluster_pt;
   Float_t         Lepton_cluster_eta;
   Float_t         Lepton_cluster_phi;
   Float_t         Lepton_cluster_E;
   Float_t         Lepton_cluster_Et;
   Float_t         Lepton_track_pt;
   Float_t         Lepton_track_eta;

 Float_t MuonSFX;
  Float_t MuonSFXError;

  Float_t UJet1_Jet_rapidity;
  Float_t BJet1_Jet_rapidity;
  Float_t Lepton_rapidity;
  Float_t TopBJet1_rapidity;


  //------------------------------------------------------
  // Tree Variables
  //------------------------------------------------------
  //informational variables first
   Int_t           BJet1_Index;
   Int_t           BJet2_Index;
   Int_t           BestJet_Index;
   Int_t           LastJet_Index;
   Int_t           NotBestJet1_Index;
   Int_t           NotBestJet2_Index;
   Int_t           UJet1_Index;
   Int_t           UJet2_Index;
   Int_t           TriangularCut1;
   Int_t           TriangularCut2;
   Int_t           TriangularCut3;

  Int_t           Jet1_Index;
  Int_t           Jet2_Index;
  Int_t           Jet3_Index;

  Int_t           Jet1eta_Index;
  Int_t           Jet2eta_Index;
  Int_t           Jet3eta_Index;

   Int_t           isBestJet;
   Int_t           isNeutrino;
  
   Int_t          RunNumber;
   Int_t          EventNumber;
   Double_t        weight;

   Int_t           Lepton_nBLHits;
   Int_t           Lepton_nPixHits;
   Int_t           Lepton_nSCTHits;
   Int_t           Lepton_nTRTHits;
   Int_t           Lepton_nTRTOutliers;
   Int_t           Lepton_nTRTHighTOutliers;
   Float_t         Lepton_track_d0;
   Float_t         Lepton_track_z0;
   Float_t         Lepton_track_err_d0;
   Float_t         Lepton_track_err_z0;
   Float_t         Lepton_track_d0_exPV;
   Float_t         Lepton_track_z0_exPV;
   Float_t         Lepton_track_err_d0_exPV;
   Float_t         Lepton_track_err_z0_exPV;

   Int_t           Lepton_pdgId;
   Float_t         Lepton_pt;
   Float_t         Lepton_eta;
   Float_t         Lepton_phi;
   Float_t         Lepton_E;
   Float_t         Lepton_Et;
   Float_t         Lepton_charge;
   Float_t         Lepton_etcone20;
   Float_t         Lepton_etcone30;
   Float_t         Lepton_etcone40;
   Float_t         Lepton_ptcone20;
   Float_t         Lepton_ptcone30;
   Float_t         Lepton_ptcone40;
   Float_t         MissingEt_etx;
   Float_t         MissingEt_ety;
   Float_t         MissingEt_et;
   Float_t         MissingEt_sumet;
   Float_t         MissingEt_phi;
   Int_t           Jet_n;
   Int_t           BJet_n;
   Int_t           BJetSV0_n;
 
   Double_t        Aplanarity_AllJetsLepton;
   Double_t        Centrality_AllJetsLepton;
   Double_t        Centrality_Jet1Jet2;
   Double_t        Centrality_LastJetLepton;
   Double_t        CosDeltaPhi_BJet1BJet2;
   Double_t        CosDeltaPhi_BJet1Lepton;
   Double_t        CosDeltaPhi_BJet1UJet1;
   Double_t        CosDeltaPhi_BJet1UJet2;
   Double_t        CosDeltaPhi_BJet2Lepton;
   Double_t        CosDeltaPhi_LeptonMissingEt;
   Double_t        CosDeltaPhi_LeptonNeutrino;
   Double_t        CosDeltaPhi_UJet1Lepton;
   Double_t        CosDeltaPhi_UJet1UJet2;
   Double_t        CosDeltaPhi_UJet2Lepton;
   Double_t        Cos_BJet1RAllJets_AllJets;
   Double_t        Cos_Jet1RAllJets_AllJets;
   Double_t        Cos_Jet1RCM_Zaxis;
   Double_t        Cos_LeptonRW_WRBestJetTop;
   Double_t        Cos_Lepton_UJet1_RBJet1Top;
   Double_t        Cos_Lepton_Zaxis_RBJet1Top;
   Double_t        Cos_Lepton_Zaxis_RBestJetTop;
   Double_t        Cos_NotBestJet1RAllJets_AllJets;
   Double_t        Cos_UJet1RAllJets_AllJets;
   Double_t        DeltaEta_BJet1BJet2;
   Double_t        DeltaEta_BJet1Lepton;
   Double_t        DeltaEta_BJet1UJet1;
   Double_t        DeltaEta_BJet1UJet2;
   Double_t        DeltaEta_BJet2Lepton;
   Double_t        DeltaEta_Jet1Lepton;
   Double_t        DeltaEta_Jet2Lepton;
   Double_t        DeltaEta_Jet3Lepton;
   Double_t        DeltaEta_Jet1Jet2;
   Double_t        DeltaEta_Jet2Jet3;
   Double_t        DeltaEta_LeptonNeutrino;
   Double_t        DeltaEta_UJet1Lepton;
   Double_t        DeltaEta_UJet1UJet2;
   Double_t        DeltaEta_UJet2Lepton;
   Double_t        DeltaPhi_BJet1BJet2;
   Double_t        DeltaPhi_BJet1Lepton;
   Double_t        DeltaPhi_BJet1UJet1;
   Double_t        DeltaPhi_BJet1UJet2;
   Double_t        DeltaPhi_BJet2Lepton;
   Double_t        DeltaPhi_Jet1Jet2;
   Double_t        DeltaPhi_Jet1Lepton;
   Double_t        DeltaPhi_Jet1MissingEt;
   Double_t        DeltaPhi_Jet1Neutrino;
   Double_t        DeltaPhi_Jet2Lepton;
   Double_t        DeltaPhi_Jet2MissingEt;
   Double_t        DeltaPhi_Jet2Neutrino;
   Double_t        DeltaPhi_Jet3Lepton;
   Double_t        DeltaPhi_Jet3MissingEt;
   Double_t        DeltaPhi_Jet3Neutrino;
   Double_t        DeltaPhi_LeptonMissingEt;
   Double_t        DeltaPhi_LeptonNeutrino;
   Double_t        DeltaPhi_UJet1Lepton;
   Double_t        DeltaPhi_UJet1UJet2;
   Double_t        DeltaPhi_UJet2Lepton;
   Double_t        DeltaPt_BJet1UJet1;
   Double_t        DeltaPt_BJet1UJet2;
   Double_t        DeltaPt_Jet1Jet2;
   Double_t        DeltaRMin_LeptonAllJets;
   Double_t        DeltaR_BJet1BJet2;
   Double_t        DeltaR_BJet1Lepton;
   Double_t        DeltaR_BJet1UJet1;
   Double_t        DeltaR_BJet1UJet2;
   Double_t        DeltaR_BJet2Lepton;
   Double_t        DeltaR_Jet1Jet2;
   Double_t        DeltaR_Jet1Lepton;
   Double_t        DeltaR_Jet2Lepton;
   Double_t        DeltaR_Jet3Lepton;
   Double_t        DeltaR_UJet1Lepton;
   Double_t        DeltaR_UJet1UJet2;
   Double_t        DeltaR_UJet2Lepton;
   Double_t        H_AllJets;
   Double_t        H_AllJetsLeptonNeutrino;
   Double_t        H_AllJetsMinusBJet1;
   Double_t        H_AllJetsMinusBestJet;
   Double_t        H_Jet1Jet2;
   Double_t        H_Jet1Jet2LeptonMissingEt;
   Double_t        H_Jet1Jet2LeptonNeutrino;
   Double_t        H_Jet1Jet2LeptonNeutrinoMinusBJet1;
   Double_t        Ht_AllJets;
   Double_t        Ht_AllJetsLepton;
   Double_t        Ht_AllJetsLeptonMissingEt;
   Double_t        Ht_AllJetsMinusBJet1;
   Double_t        Ht_AllJetsMinusBestJet;
   Double_t        Ht_Jet1Jet2;
   Double_t        Ht_Jet1Jet2LeptonMissingEt;
   Double_t        Ht_Jet1Jet2LeptonNeutrino;
   Double_t        Ht_Jet1Jet2LeptonNeutrinoMinusBJet1;
   Double_t        Ht_LeptonMissingEt;
   Double_t        Ht_LeptonNeutrino;
   Double_t        InvariantMass_AllJets;
   Double_t        InvariantMass_AllJetsLeptonNeutrino;
   Double_t        InvariantMass_AllJetsMinusBJet1;
   Double_t        InvariantMass_AllJetsMinusBestJet;
   Double_t        InvariantMass_BestJetLeptonNeutrino;
   Double_t        Jet1Jet2LeptonNeutrino_InvariantMass;
   Double_t        Jet1Jet2LeptonNeutrino_Pt;
   Double_t        Jet1Jet2_InvariantMass;
   Double_t        Jet1Jet2_TransverseMass;
   Double_t        Jet1Jet2_p;
   Double_t        Jet1Jet2_pt;
   Double_t        Mass_BJet1BJet2;
   Double_t        MaxJetEta;
   Double_t        MinJetEta;
   Double_t        Neutrino_pz;
   Double_t        P_AllJets;
   Double_t        Pt_AllJetsMinusBJet1;
   Double_t        Pt_AllJetsMinusBestJet;
   Double_t        Pt_BJet1BJet2;
   Double_t        Sphericity_AllJetsLepton;
   Double_t        TopBJet1_eta;
   Double_t        TopBJet1_mass;
   Double_t        TopBJet1_phi;
   Double_t        TopBJet1_pt;
   Double_t        TopBJet2_eta;
   Double_t        TopBJet2_mass;
   Double_t        TopBJet2_phi;
   Double_t        TopBJet2_pt;
   Double_t        TopBestJet_eta;
   Double_t        TopBestJet_mass;
   Double_t        TopBestJet_phi;
   Double_t        TopBestJet_pt;
   Double_t        TopJet1_eta;
   Double_t        TopJet1_mass;
   Double_t        TopJet1_phi;
   Double_t        TopJet1_pt;
   Double_t        TopJet2_eta;
   Double_t        TopJet2_mass;
   Double_t        TopJet2_phi;
   Double_t        TopJet2_pt;
   Double_t        TopJet3_eta;
   Double_t        TopJet3_mass;
   Double_t        TopJet3_phi;
   Double_t        TopJet3_pt;
   Double_t        W_TransverseMass;
   Double_t        W_charge;
   Double_t        W_eta;
   Double_t        W_mass;
   Double_t        W_pt;
   Double_t        WprimefinalBJet1TopBJet2_mass;
   Double_t        WprimefinalBJet1TopBJet2_eta; 
   Double_t        WprimefinalBJet1TopBJet2_phi; 
   Double_t        WprimefinalBJet1TopBJet2_pt;  
   Double_t        WprimefinalBJet2TopBJet1_mass;
   Double_t        WprimefinalBJet2TopBJet1_eta; 
   Double_t        WprimefinalBJet2TopBJet1_phi; 
   Double_t        WprimefinalBJet2TopBJet1_pt;  
   Double_t        WprimefinalJet1TopJet2_mass;
   Double_t        WprimefinalJet1TopJet2_eta; 
   Double_t        WprimefinalJet1TopJet2_phi; 
   Double_t        WprimefinalJet1TopJet2_pt;  
   Double_t        WprimefinalJet2TopJet1_mass;
   Double_t        WprimefinalJet2TopJet1_eta; 
   Double_t        WprimefinalJet2TopJet1_phi; 
   Double_t        WprimefinalJet2TopJet1_pt;  

  //Variables from vectors
  Int_t          BJet1_BJet_flavor_track_n;
  Int_t          BJet1_BJet_flavor_svp_ntrk;
  Int_t          BJet1_BJet_flavor_svp_ntrkv;
  Int_t          BJet1_BJet_flavor_svp_ntrkj;

  Int_t          BJet2_BJet_flavor_track_n;
  Int_t          BJet2_BJet_flavor_svp_ntrk;
  Int_t          BJet2_BJet_flavor_svp_ntrkv;
  Int_t          BJet2_BJet_flavor_svp_ntrkj;


  Float_t          BJet1_nTrk;
  Float_t          BJet1_sumPtTrk;
  Float_t          BJet1_flavor_sv0p_mass;
  Float_t          BJet1_flavor_sv0p_efrc;
  Float_t          BJet1_flavor_sv0p_x;
  Float_t          BJet1_flavor_sv0p_y;
  Float_t          BJet1_flavor_sv0p_z;
  Float_t          BJet1_flavor_sv0p_chi2;

 Float_t          UJet1_nTrk;
  Float_t          UJet1_sumPtTrk;
  Float_t          UJet1_flavor_sv0p_mass;
  Float_t          UJet1_flavor_sv0p_efrc;
  Float_t          UJet1_flavor_sv0p_x;
  Float_t          UJet1_flavor_sv0p_y;
  Float_t          UJet1_flavor_sv0p_z;
  Float_t          UJet1_flavor_sv0p_chi2;

 Float_t          BJet2_nTrk;
  Float_t          BJet2_sumPtTrk;
  Float_t          BJet2_flavor_sv0p_mass;
  Float_t          BJet2_flavor_sv0p_efrc;
  Float_t          BJet2_flavor_sv0p_x;
  Float_t          BJet2_flavor_sv0p_y;
  Float_t          BJet2_flavor_sv0p_z;
  Float_t          BJet2_flavor_sv0p_chi2;

 Float_t          UJet2_nTrk;
  Float_t          UJet2_sumPtTrk;
  Float_t          UJet2_flavor_sv0p_mass;
  Float_t          UJet2_flavor_sv0p_efrc;
  Float_t          UJet2_flavor_sv0p_x;
  Float_t          UJet2_flavor_sv0p_y;
  Float_t          UJet2_flavor_sv0p_z;
  Float_t          UJet2_flavor_sv0p_chi2;
  
   
  Float_t          Jet1_Jet_pt;
   Float_t          Jet1_Jet_eta;
   Float_t          Jet1_Jet_phi;
   Float_t          Jet1_Jet_E;
   Float_t          Jet1_Jet_Et;
   Float_t          Jet1_Jet_charge;

  Float_t          Jet2_Jet_pt;
   Float_t          Jet2_Jet_eta;
   Float_t          Jet2_Jet_phi;
   Float_t          Jet2_Jet_E;
   Float_t          Jet2_Jet_Et;
   Float_t          Jet2_Jet_charge;

  Float_t          Jet3_Jet_pt;
   Float_t          Jet3_Jet_eta;
   Float_t          Jet3_Jet_phi;
   Float_t          Jet3_Jet_E;
   Float_t          Jet3_Jet_Et;
   Float_t          Jet3_Jet_charge;

  Float_t          Jet4_Jet_pt;
   Float_t          Jet4_Jet_eta;
   Float_t          Jet4_Jet_phi;
   Float_t          Jet4_Jet_E;
   Float_t          Jet4_Jet_Et;
   Float_t          Jet4_Jet_charge;

  Float_t          Jet1eta_Jet_pt;
   Float_t          Jet1eta_Jet_eta;
   Float_t          Jet1eta_Jet_phi;
   Float_t          Jet1eta_Jet_E;
   Float_t          Jet1eta_Jet_Et;
   Float_t          Jet1eta_Jet_charge;

  Float_t          Jet2eta_Jet_pt;
   Float_t          Jet2eta_Jet_eta;
   Float_t          Jet2eta_Jet_phi;
   Float_t          Jet2eta_Jet_E;
   Float_t          Jet2eta_Jet_Et;
   Float_t          Jet2eta_Jet_charge;

  Float_t          Jet3eta_Jet_pt;
   Float_t          Jet3eta_Jet_eta;
   Float_t          Jet3eta_Jet_phi;
   Float_t          Jet3eta_Jet_E;
   Float_t          Jet3eta_Jet_Et;
   Float_t          Jet3eta_Jet_charge;

   Float_t          BJet1_Jet_pt;
   Float_t          BJet1_Jet_eta;
   Float_t          BJet1_Jet_phi;
   Float_t          BJet1_Jet_E;
   Float_t          BJet1_Jet_Et;
   Float_t          BJet1_Jet_charge;
    
   Float_t          BJet2_Jet_pt;
   Float_t          BJet2_Jet_eta;
   Float_t          BJet2_Jet_phi;
   Float_t          BJet2_Jet_E;
   Float_t          BJet2_Jet_Et;
   Float_t          BJet2_Jet_charge;

   Float_t          BestJet_Jet_pt;
   Float_t          BestJet_Jet_eta;
   Float_t          BestJet_Jet_phi;
   Float_t          BestJet_Jet_E;
   Float_t          BestJet_Jet_Et;
   Float_t          BestJet_Jet_charge;

   Float_t          LastJet_Jet_pt;
  Float_t          LastJet_Jet_p;
   Float_t          LastJet_Jet_eta;
   Float_t          LastJet_Jet_phi;
   Float_t          LastJet_Jet_E;
   Float_t          LastJet_Jet_Et;
   Float_t          LastJet_Jet_charge;

   Float_t          NotBestJet1_Jet_pt;
   Float_t          NotBestJet1_Jet_eta;
   Float_t          NotBestJet1_Jet_phi;
   Float_t          NotBestJet1_Jet_E;
   Float_t          NotBestJet1_Jet_Et;
   Float_t          NotBestJet1_Jet_charge;

   Float_t          NotBestJet2_Jet_pt;
   Float_t          NotBestJet2_Jet_eta;
   Float_t          NotBestJet2_Jet_phi;
   Float_t          NotBestJet2_Jet_E;
   Float_t          NotBestJet2_Jet_Et;
   Float_t          NotBestJet2_Jet_charge;

   Float_t          UJet1_Jet_pt;
   Float_t          UJet1_Jet_eta;
   Float_t          UJet1_Jet_phi;
   Float_t          UJet1_Jet_E;
   Float_t          UJet1_Jet_Et;
   Float_t          UJet1_Jet_charge;


   Float_t          UJet2_Jet_pt;
   Float_t          UJet2_Jet_eta;
   Float_t          UJet2_Jet_phi;
   Float_t          UJet2_Jet_E;
   Float_t          UJet2_Jet_Et;
   Float_t          UJet2_Jet_charge;
  /**
  //------------------------------------------------------
  // System Mass Histograms
  //------------------------------------------------------
  Float_t InvariantMass_AllJets;
  Float_t InvariantMass_Jet1Jet2LeptonMET;
  Float_t InvariantMass_AllJetsLeptonMET;
  Float_t InvariantMass_LeptonMETBestJet;
  Float_t InvariantMass_AllJets_MinusBestJet;
  Float_t InvariantMass_AllJets_MinusBTaggedJet;
  Float_t InvariantMass_Jet1Jet2;
  Float_t TransverseMass_Jet1Jet2;
  Float_t Shat;
  //Float_t Neutrino_Mass;
  Float_t Neutrino_Pz;
    
  //------------------------------------------------------
  // H Histograms
  //------------------------------------------------------
  Float_t H_Jet1Jet2LeptonMET;
  Float_t H_Jet1Jet2LeptonMET_MinusBJet;
  Float_t H_AllJetsLeptonMET;
  Float_t H_AllJets;
  Float_t H_AllJets_MinusBestJet;
  Float_t H_AllJets_MinusBTaggedJet;
  Float_t H_Jet1Jet2;
    
  //------------------------------------------------------
  // Ht Histograms
  //------------------------------------------------------
  Float_t HT; //I added this
  Float_t HT_Jet1Jet2LeptonMET;
  Float_t HT_LeptonMET;
  Float_t HT_Jet1Jet2LeptonMET_MinusBJet;
  Float_t HT_AllJetsLeptonMET;
  Float_t HT_AllJetsLepton;
  Float_t HT_AllJets;
  Float_t HT_LastJet;
  Float_t P_LastJet;
  Float_t HT_AllJets_MinusBestJet;
  Float_t HT_AllJets_MinusBTaggedJet;
  Float_t HT_Jet1Jet2;

  //------------------------------------------------------
  // Jet Histograms
  //------------------------------------------------------
  Float_t Jet1Pt;
  Float_t Jet2Pt;
  Float_t Jet3Pt;
  Float_t Jet4Pt;
  Float_t Jet1Eta;
  Float_t Jet2Eta;
  Float_t Jet3Eta;
  Float_t Jet4Eta;
  Float_t Jet1Phi;
  Float_t Jet2Phi;
  Float_t Jet3Phi;
  Float_t Jet4Phi;
  Float_t Jet1Pt_NotBest;
  Float_t Jet2Pt_NotBest;
  Float_t Pt_Jet1Jet2;
  Float_t Pt_AllJets_MinusBestJet;
  Float_t Pt_AllJets_MinusBTaggedJet;
  Float_t DeltaPt_Jet1Jet2;
  int NJets;
  //int NGoodJets;
  int NTaggedJets;
  int NUntaggedJets;
  Float_t BestJetPt;
  Float_t DeltaRJet1Jet2;
  Float_t DeltaPhiLeadingJetLepton;
  Float_t DeltaPhiLeadingJetMissingEt;
  Float_t DeltaPhiLeptonMissingEt;
  Float_t DeltaRLeadingJetLepton;
  Float_t Jet12_Centrality;
  Float_t AllJetsLepton_Centrality;
  Float_t LastJetLepton_Centrality;
  Float_t MaxJetEta;
  Float_t MinJetEta;

  Float_t JetProbWeight1;
  Float_t JetProbWeight2;
  Float_t JetProbWeight3;
  Float_t JetProbWeight4;

  Float_t SV0Weight1;
  Float_t SV0Weight2;
  Float_t SV0Weight3;
  Float_t SV0Weight4;

  Int_t JetTruthLabel;
  //------------------------------------------------------
  // Missing Et Histograms
  //------------------------------------------------------
  // Float_t METPt;
  Float_t MET;
    
  //------------------------------------------------------
  // W Histograms
  //------------------------------------------------------
  Float_t WTransverseMass;
  Float_t WPt;
  int WCharge;
  
  //------------------------------------------------------
  // Top Histograms
  //------------------------------------------------------
  Float_t BestTopMass;
  Float_t BTaggedTopMass;
  Float_t LeadingJetTopMass;
    
  //------------------------------------------------------
  // Lepton Histograms
  //------------------------------------------------------
  Float_t LeptonPt;
  Float_t LeptonEta;
  //  Float_t QTimesEtaLepton;
  Float_t LeptonChargePDGID;
  //Float_t EtCone20;
  //------- 
  Int_t LeptonIsTight;
  Float_t LeptonEtCone20;
  Float_t LeptonEtCone30;
  Float_t LeptonEtCone40;
  Float_t LeptonPtCone30;
  Int_t  ElectronQuality;
  Int_t  IsRobustMedium;
  Int_t IsemTrackBlayer;
  Int_t IsemTrackMatchEoverP;
  //------------------------------------------------------
  // UnTagged Jet Histograms
  //------------------------------------------------------
  Float_t LeadingUntaggedJetPt;
  Float_t SecondUntaggedJetPt;
  Float_t LeadingUntaggedJetEta;
  // Float_t QTimesEta;
  Float_t DeltaPt_tag_untag;
  Float_t DeltaRLeadingUntaggedJetLepton;
  //------------------------------------------------------
  // b-Tagged Jet Histograms
  //------------------------------------------------------
  Float_t LeadingBTaggedJetPt;
  Float_t SecondBTaggedJetPt;
  Float_t LeadingBTaggedJetEta;
  Float_t SecondBTaggedJetEta;
  Float_t DeltaRLeadingBTaggedJetLepton;
  Float_t DeltaRSecondBTaggedJetLepton;
  Float_t DeltaRBTaggedJet1BTaggedJet2;
  //------------------------------------------------------
  // Angular Variable Histograms
  //------------------------------------------------------
  Float_t Cos_UntaggedJetLepton_BTaggedTop;
  Float_t Cos_Jet1AllJets_AllJets;
  Float_t Cos_BTaggedJetAllJets_AllJets;
  Float_t Cos_UntaggedJetAllJets_AllJets;
  Float_t Cos_NotBestJetAllJets_AllJets;
  Float_t Cos_LeptonQZ_BestTop;
  Float_t Cos_wHelicityThetaStar;
  Float_t Cos_Jet1ThetaStar;
  **/

 
};


#endif


