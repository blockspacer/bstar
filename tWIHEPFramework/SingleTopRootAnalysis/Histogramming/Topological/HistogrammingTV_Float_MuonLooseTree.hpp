/******************************************************************************
 * HistogrammingTV_Float_MuonLooseTree.hpp                                                  *
 *                                                                            *
 * Fills a tree with specific variables.                                      *
 * Should be called after preselection cuts in bin/analysis                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTV_Float_MuonLooseTree()              -- Parameterized Constructor       *
 *    ~HistogrammingTV_Float_MuonLooseTree()             -- Destructor                      *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTV_Float_MuonLooseTree"    *
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

#ifndef HistogrammingTV_Float_MuonLooseTree_h
#define HistogrammingTV_Float_MuonLooseTree_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;


class HistogrammingTV_Float_MuonLooseTree : public HistoCut 
{

public:

  HistogrammingTV_Float_MuonLooseTree(EventContainer *obj,TFile *skimFilez);
  
  // destructor
  ~HistogrammingTV_Float_MuonLooseTree();

  // methods:
  void BookHistogram();
  bool Apply();
  
  inline std::string GetCutName (void) { return "HistogrammingTV_Float_MuonLooseTree"; };

private:
  TFile *_skimFile;
  TTree *_skimTree;
  std::vector<Double_t> _bweightvec;

  // skim tree variables
 
  //------------------------------------------------------
  // Event Weight
  //------------------------------------------------------
  Float_t _EventWeight;
  Int_t _HFOR;
   Int_t Pvtx_n;
 Int_t pdgid;
 
  Int_t _NWDecayLeptons;
  Int_t _NWDecayLightQuarks;

 Float_t _PileupWeight; //tagging weight per event
  Float_t _TaggingWeight; //tagging weight per event
  Float_t _TaggingWeight_Bup;//tagging weight per event applying the b SF up shift
  Float_t _TaggingWeight_Bdown;//tagging weight per event applying the b SF down shift
  Float_t _TaggingWeight_Lqup;//tagging weight per event applying the Lq SF up shift
  Float_t _TaggingWeight_Lqdown;//tagging weight per event applying the Lq SF down shift





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
   Float_t          BJet1_Jet_emfraction;
   Float_t          BJet1_Jet_jvtxf;
    
   Float_t          BJet2_Jet_pt;
   Float_t          BJet2_Jet_eta;
   Float_t          BJet2_Jet_phi;
   Float_t          BJet2_Jet_E;
   Float_t          BJet2_Jet_Et;
   Float_t          BJet2_Jet_charge;
   Float_t          BJet2_Jet_emfraction;
   Float_t          BJet2_Jet_jvtxf;

   Float_t          BestJet_Jet_pt;
   Float_t          BestJet_Jet_eta;
   Float_t          BestJet_Jet_phi;
   Float_t          BestJet_Jet_E;
   Float_t          BestJet_Jet_Et;
   Float_t          BestJet_Jet_charge;
   Float_t          BestJet_Jet_emfraction;
   Float_t          BestJet_Jet_jvtxf;

   Float_t          LastJet_Jet_pt;
   Float_t          LastJet_Jet_eta;
   Float_t          LastJet_Jet_phi;
   Float_t          LastJet_Jet_E;
   Float_t          LastJet_Jet_Et;
   Float_t          LastJet_Jet_charge;
   Float_t          LastJet_Jet_emfraction;
   Float_t          LastJet_Jet_jvtxf;

   Float_t          NotBestJet1_Jet_pt;
   Float_t          NotBestJet1_Jet_eta;
   Float_t          NotBestJet1_Jet_phi;
   Float_t          NotBestJet1_Jet_E;
   Float_t          NotBestJet1_Jet_Et;
   Float_t          NotBestJet1_Jet_charge;
   Float_t          NotBestJet1_Jet_emfraction;
   Float_t          NotBestJet1_Jet_jvtxf;

   Float_t          NotBestJet2_Jet_pt;
   Float_t          NotBestJet2_Jet_eta;
   Float_t          NotBestJet2_Jet_phi;
   Float_t          NotBestJet2_Jet_E;
   Float_t          NotBestJet2_Jet_Et;
   Float_t          NotBestJet2_Jet_charge;
   Float_t          NotBestJet2_Jet_emfraction;
   Float_t          NotBestJet2_Jet_jvtxf;

   Float_t          UJet1_Jet_pt;
   Float_t          UJet1_Jet_eta;
   Float_t          UJet1_Jet_phi;
   Float_t          UJet1_Jet_E;
   Float_t          UJet1_Jet_Et;
   Float_t          UJet1_Jet_charge;
   Float_t          UJet1_Jet_emfraction;
   Float_t          UJet1_Jet_jvtxf;


   Float_t          UJet2_Jet_pt;
   Float_t          UJet2_Jet_eta;
   Float_t          UJet2_Jet_phi;
   Float_t          UJet2_Jet_E;
   Float_t          UJet2_Jet_Et;
   Float_t          UJet2_Jet_charge;
   Float_t          UJet2_Jet_emfraction;
   Float_t          UJet2_Jet_jvtxf;


  /**
   std::vector<float>   *Jet_pt;
   std::vector<float>   *Jet_eta;
   std::vector<float>   *Jet_phi;
   std::vector<float>   *Jet_E;
   std::vector<float>   *Jet_Et;
   std::vector<float>   *Jet_charge;
   std::vector<float>   *Jet_emfraction;
   std::vector<float>   *Jet_jvtxf;

   std::vector<int>     *BJet_flavor_track_n;
   std::vector<int>     *BJet_flavor_svp_ntrk;
   std::vector<int>     *BJet_flavor_svp_ntrkv;
   std::vector<int>     *BJet_flavor_svp_ntrkj;
   
   std::vector<float>   *BJet_pt;
   std::vector<float>   *BJet_eta;
   std::vector<float>   *BJet_phi;
   std::vector<float>   *BJet_E;
   std::vector<float>   *BJet_Et;
   std::vector<float>   *BJet_charge;
   std::vector<float>   *BJet_emfraction;
   std::vector<float>   *BJet_jvtxf;
 
  //theses are variables not included for now:
  std::vector<bool>    *Jet_isGood;
   std::vector<bool>    *Jet_isUgly;
   std::vector<bool>    *Jet_isBad;
   std::vector<int>     *Jet_flavor_track_n;
   std::vector<int>     *Jet_flavor_svp_n2t;
   std::vector<int>     *Jet_flavor_svp_ntrk;
   std::vector<int>     *Jet_flavor_svp_ntrkv;
   std::vector<int>     *Jet_flavor_svp_ntrkj;
   std::vector<int>     *Jet_flavor_svp_ndof;
   std::vector<int>     *Jet_flavor_sv0p_n2t;
   std::vector<int>     *Jet_flavor_sv0p_ntrk;
   std::vector<int>     *Jet_flavor_sv0p_ntrkv;
   std::vector<int>     *Jet_flavor_sv0p_ntrkj;
   std::vector<int>     *Jet_flavor_sv0p_ndof;
   std::vector<int>     *Jet_flavor_truth_trueflav;
   std::vector<int>     *Jet_flavor_truth_BHadronpdg;
   std::vector<float>   *Jet_jvtxf;
   std::vector<float>   *Jet_flavor_svp_mass;
   std::vector<float>   *Jet_flavor_svp_efrc;
   std::vector<float>   *Jet_flavor_svp_x;
   std::vector<float>   *Jet_flavor_svp_y;
   std::vector<float>   *Jet_flavor_svp_z;
   std::vector<float>   *Jet_flavor_svp_chi2;
   std::vector<float>   *Jet_flavor_svp_err_x;
   std::vector<float>   *Jet_flavor_svp_err_y;
   std::vector<float>   *Jet_flavor_svp_err_z;
   std::vector<float>   *Jet_flavor_svp_cov_x;
   std::vector<float>   *Jet_flavor_svp_cov_y;
   std::vector<float>   *Jet_flavor_svp_cov_z;
   std::vector<float>   *Jet_flavor_svp_cov_xy;
   std::vector<float>   *Jet_flavor_svp_cov_xz;
   std::vector<float>   *Jet_flavor_svp_cov_yz;
   std::vector<float>   *Jet_flavor_sv0p_mass;
   std::vector<float>   *Jet_flavor_sv0p_efrc;
   std::vector<float>   *Jet_flavor_sv0p_x;
   std::vector<float>   *Jet_flavor_sv0p_y;
   std::vector<float>   *Jet_flavor_sv0p_z;
   std::vector<float>   *Jet_flavor_sv0p_chi2;
   std::vector<float>   *Jet_flavor_sv0p_err_x;
   std::vector<float>   *Jet_flavor_sv0p_err_y;
   std::vector<float>   *Jet_flavor_sv0p_err_z;
   std::vector<float>   *Jet_flavor_sv0p_cov_x;
   std::vector<float>   *Jet_flavor_sv0p_cov_y;
   std::vector<float>   *Jet_flavor_sv0p_cov_z;
   std::vector<float>   *Jet_flavor_sv0p_cov_xy;
   std::vector<float>   *Jet_flavor_sv0p_cov_xz;
   std::vector<float>   *Jet_flavor_sv0p_cov_yz;
   std::vector<double>  *Jet_flavor_weight;
   std::vector<double>  *Jet_flavor_weight_TrackCounting2D;
   std::vector<double>  *Jet_flavor_weight_JetProb;
   std::vector<double>  *Jet_flavor_weight_IP1D;
   std::vector<double>  *Jet_flavor_weight_IP2D;
   std::vector<double>  *Jet_flavor_weight_IP3D;
   std::vector<double>  *Jet_flavor_weight_SV0;
   std::vector<double>  *Jet_flavor_weight_SV1;
   std::vector<double>  *Jet_flavor_weight_SV2;
   std::vector<double>  *Jet_flavor_weight_JetFitterTag;
   std::vector<double>  *Jet_flavor_weight_JetFitterCOMB;
   std::vector<double>  *Jet_flavor_weight_JetFitterCOMBNN;
   std::vector<double>  *Jet_flavor_weight_SoftMuonTag;
   std::vector<double>  *Jet_flavor_weight_SoftElectronTag;
   std::vector<double>  *Jet_flavor_truth_dRMinToB;
   std::vector<double>  *Jet_flavor_truth_dRMinToC;
   std::vector<double>  *Jet_flavor_truth_dRMinToT;

  std::vector<bool>    *BJet_isGood;
  std::vector<bool>    *BJet_isUgly;
   std::vector<bool>    *BJet_isBad;
  
  std::vector<int>     *BJet_flavor_svp_n2t;
   std::vector<int>     *BJet_flavor_svp_ntrk;
   std::vector<int>     *BJet_flavor_svp_ntrkv;
   std::vector<int>     *BJet_flavor_svp_ntrkj;
   std::vector<int>     *BJet_flavor_svp_ndof;
   std::vector<int>     *BJet_flavor_sv0p_n2t;
   std::vector<int>     *BJet_flavor_sv0p_ntrk;
   std::vector<int>     *BJet_flavor_sv0p_ntrkv;
   std::vector<int>     *BJet_flavor_sv0p_ntrkj;
   std::vector<int>     *BJet_flavor_sv0p_ndof;
   std::vector<int>     *BJet_flavor_truth_trueflav;
   std::vector<int>     *BJet_flavor_truth_BHadronpdg;
  std::vector<float>   *BJet_flavor_svp_mass;
   std::vector<float>   *BJet_flavor_svp_efrc;
   std::vector<float>   *BJet_flavor_svp_x;
   std::vector<float>   *BJet_flavor_svp_y;
   std::vector<float>   *BJet_flavor_svp_z;
   std::vector<float>   *BJet_flavor_svp_chi2;
   std::vector<float>   *BJet_flavor_svp_err_x;
   std::vector<float>   *BJet_flavor_svp_err_y;
   std::vector<float>   *BJet_flavor_svp_err_z;
   std::vector<float>   *BJet_flavor_svp_cov_x;
   std::vector<float>   *BJet_flavor_svp_cov_y;
   std::vector<float>   *BJet_flavor_svp_cov_z;
   std::vector<float>   *BJet_flavor_svp_cov_xy;
   std::vector<float>   *BJet_flavor_svp_cov_xz;
   std::vector<float>   *BJet_flavor_svp_cov_yz;
   std::vector<float>   *BJet_flavor_sv0p_mass;
   std::vector<float>   *BJet_flavor_sv0p_efrc;
   std::vector<float>   *BJet_flavor_sv0p_x;
   std::vector<float>   *BJet_flavor_sv0p_y;
   std::vector<float>   *BJet_flavor_sv0p_z;
   std::vector<float>   *BJet_flavor_sv0p_chi2;
   std::vector<float>   *BJet_flavor_sv0p_err_x;
   std::vector<float>   *BJet_flavor_sv0p_err_y;
   std::vector<float>   *BJet_flavor_sv0p_err_z;
   std::vector<float>   *BJet_flavor_sv0p_cov_x;
   std::vector<float>   *BJet_flavor_sv0p_cov_y;
   std::vector<float>   *BJet_flavor_sv0p_cov_z;
   std::vector<float>   *BJet_flavor_sv0p_cov_xy;
   std::vector<float>   *BJet_flavor_sv0p_cov_xz;
   std::vector<float>   *BJet_flavor_sv0p_cov_yz;
   std::vector<double>  *BJet_flavor_weight;
   std::vector<double>  *BJet_flavor_weight_TrackCounting2D;
   std::vector<double>  *BJet_flavor_weight_JetProb;
   std::vector<double>  *BJet_flavor_weight_IP1D;
   std::vector<double>  *BJet_flavor_weight_IP2D;
   std::vector<double>  *BJet_flavor_weight_IP3D;
   std::vector<double>  *BJet_flavor_weight_SV0;
   std::vector<double>  *BJet_flavor_weight_SV1;
   std::vector<double>  *BJet_flavor_weight_SV2;
   std::vector<double>  *BJet_flavor_weight_JetFitterTag;
   std::vector<double>  *BJet_flavor_weight_JetFitterCOMB;
   std::vector<double>  *BJet_flavor_weight_JetFitterCOMBNN;
   std::vector<double>  *BJet_flavor_weight_SoftMuonTag;
   std::vector<double>  *BJet_flavor_weight_SoftElectronTag;
   std::vector<double>  *BJet_flavor_truth_dRMinToB;
   std::vector<double>  *BJet_flavor_truth_dRMinToC;
   std::vector<double>  *BJet_flavor_truth_dRMinToT;
  **/
  //------------------------------------------------------
  // System Mass Histograms
  //------------------------------------------------------
  Float_t _InvariantMass_AllJets;
  Float_t _InvariantMass_Jet1Jet2LeptonMET;
  Float_t _InvariantMass_AllJetsLeptonMET;
  Float_t _InvariantMass_LeptonMETBestJet;
  Float_t _InvariantMass_AllJets_MinusBestJet;
  Float_t _InvariantMass_AllJets_MinusBTaggedJet;
  Float_t _InvariantMass_Jet1Jet2;
  Float_t _TransverseMass_Jet1Jet2;
  Float_t _Shat;
  //Float_t _Neutrino_Mass;
  Float_t _Neutrino_Pz;
    
  //------------------------------------------------------
  // H Histograms
  //------------------------------------------------------
  Float_t _H_Jet1Jet2LeptonMET;
  Float_t _H_Jet1Jet2LeptonMET_MinusBJet;
  Float_t _H_AllJetsLeptonMET;
  Float_t _H_AllJets;
  Float_t _H_AllJets_MinusBestJet;
  Float_t _H_AllJets_MinusBTaggedJet;
  Float_t _H_Jet1Jet2;
    
  //------------------------------------------------------
  // Ht Histograms
  //------------------------------------------------------
  Float_t _HT; //I added this
  Float_t _HT_Jet1Jet2LeptonMET;
  Float_t _HT_LeptonMET;
  Float_t _HT_Jet1Jet2LeptonMET_MinusBJet;
  Float_t _HT_AllJetsLeptonMET;
  Float_t _HT_AllJetsLepton;
  Float_t _HT_AllJets;
  Float_t _HT_LastJet;
  Float_t _P_LastJet;
  Float_t _HT_AllJets_MinusBestJet;
  Float_t _HT_AllJets_MinusBTaggedJet;
  Float_t _HT_Jet1Jet2;

  //------------------------------------------------------
  // Jet Histograms
  //------------------------------------------------------
  Float_t _Jet1Pt;
  Float_t _Jet2Pt;
  Float_t _Jet3Pt;
  Float_t _Jet4Pt;
  Float_t _Jet1Eta;
  Float_t _Jet2Eta;
  Float_t _Jet3Eta;
  Float_t _Jet4Eta;
  Float_t _Jet1Phi;
  Float_t _Jet2Phi;
  Float_t _Jet3Phi;
  Float_t _Jet4Phi;
  Float_t _Jet1Pt_NotBest;
  Float_t _Jet2Pt_NotBest;
  Float_t _Pt_Jet1Jet2;
  Float_t _Pt_AllJets_MinusBestJet;
  Float_t _Pt_AllJets_MinusBTaggedJet;
  Float_t _DeltaPt_Jet1Jet2;
  int _NJets;
  //int _NGoodJets;
  int _NTaggedJets;
  int _NUntaggedJets;
  Float_t _BestJetPt;
  Float_t _DeltaRJet1Jet2;
  Float_t _DeltaPhiLeadingJetLepton;
  Float_t _DeltaPhiLeadingJetMissingEt;
  Float_t _DeltaPhiLeptonMissingEt;
  Float_t _DeltaRLeadingJetLepton;
  Float_t _Jet12_Centrality;
  Float_t _AllJetsLepton_Centrality;
  Float_t _LastJetLepton_Centrality;
  Float_t _MaxJetEta;
  Float_t _MinJetEta;

  Float_t _JetProbWeight1;
  Float_t _JetProbWeight2;
  Float_t _JetProbWeight3;
  Float_t _JetProbWeight4;

  Float_t _SV0Weight1;
  Float_t _SV0Weight2;
  Float_t _SV0Weight3;
  Float_t _SV0Weight4;

  Int_t _JetTruthLabel;
  //------------------------------------------------------
  // Missing Et Histograms
  //------------------------------------------------------
  // Float_t _METPt;
  Float_t _MET;
    
  //------------------------------------------------------
  // W Histograms
  //------------------------------------------------------
  Float_t _WTransverseMass;
  Float_t _WPt;
  int _WCharge;
  
  //------------------------------------------------------
  // Top Histograms
  //------------------------------------------------------
  Float_t _BestTopMass;
  Float_t _BTaggedTopMass;
  Float_t _LeadingJetTopMass;
    
  //------------------------------------------------------
  // Lepton Histograms
  //------------------------------------------------------
  Float_t _LeptonPt;
  Float_t _LeptonEta;
  //  Float_t _QTimesEtaLepton;
  Float_t _LeptonChargePDGID;
  //Float_t _EtCone20;
  //------- 
  Int_t _LeptonIsTight;
  Float_t _LeptonEtCone20;
  Float_t _LeptonEtCone30;
  Float_t _LeptonEtCone40;
  Float_t _LeptonPtCone30;
  Int_t  _ElectronQuality;
  Int_t  _IsRobustMedium;
  Int_t _IsemTrackBlayer;
  Int_t _IsemTrackMatchEoverP;
  //------------------------------------------------------
  // UnTagged Jet Histograms
  //------------------------------------------------------
  Float_t _LeadingUntaggedJetPt;
  Float_t _SecondUntaggedJetPt;
  Float_t _LeadingUntaggedJetEta;
  // Float_t _QTimesEta;
  Float_t _DeltaPt_tag_untag;
  Float_t _DeltaRLeadingUntaggedJetLepton;
  //------------------------------------------------------
  // b-Tagged Jet Histograms
  //------------------------------------------------------
  Float_t _LeadingBTaggedJetPt;
  Float_t _SecondBTaggedJetPt;
  Float_t _LeadingBTaggedJetEta;
  Float_t _SecondBTaggedJetEta;
  Float_t _DeltaRLeadingBTaggedJetLepton;
  Float_t _DeltaRSecondBTaggedJetLepton;
  Float_t _DeltaRBTaggedJet1BTaggedJet2;
  //------------------------------------------------------
  // Angular Variable Histograms
  //------------------------------------------------------
  Float_t _Cos_UntaggedJetLepton_BTaggedTop;
  Float_t _Cos_Jet1AllJets_AllJets;
  Float_t _Cos_BTaggedJetAllJets_AllJets;
  Float_t _Cos_UntaggedJetAllJets_AllJets;
  Float_t _Cos_NotBestJetAllJets_AllJets;
  Float_t _Cos_LeptonQZ_BestTop;
  Float_t _Cos_wHelicityThetaStar;
  Float_t _Cos_Jet1ThetaStar;


 
};


#endif


