//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 23 15:27:36 2011 by ROOT version 5.28/00a
// from TTree ElectronPreTagTree/ElectronPreTagTree
// found on file: /msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root
//////////////////////////////////////////////////////////

#ifndef ElectronPreTagTree_h
#define ElectronPreTagTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ElectronPreTagTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          Lepton_loose;
   UInt_t          Lepton_medium;
   UInt_t          Lepton_mediumIso;
   UInt_t          Lepton_tight;
   UInt_t          Lepton_tightIso;
   UInt_t          Lepton_trackBlayer;
   UInt_t          Lepton_trackMatchEoverP;
   UInt_t          Lepton_trackTRThits;
   UInt_t          Lepton_trackTRTratio;
   UInt_t          Lepton_trackA0Tight;
   UInt_t          Lepton_trackMatchEtaTight;
   UInt_t          Lepton_trackMatchPhi;
   UInt_t          Lepton_conversionMatch;
   Int_t           Lepton_Barcode;
   Int_t           Lepton_pdgId;
   Float_t         Lepton_charge;
   Float_t         Lepton_Etcone45;
   Float_t         Lepton_Etcone20;
   Float_t         Lepton_Etcone30;
   Float_t         Lepton_Etcone40;
   Float_t         Lepton_ptcone20;
   Float_t         Lepton_ptcone30;
   Float_t         Lepton_ptcone40;
   Float_t         Lepton_cluster_pt;
   Float_t         Lepton_cluster_eta;
   Float_t         Lepton_cluster_phi;
   Float_t         Lepton_cluster_E;
   Float_t         Lepton_cluster_Et;
   Float_t         Lepton_cluster_track_Et;
   Float_t         Lepton_track_d0;
   Float_t         Lepton_track_z0;
   Float_t         Lepton_track_phi;
   Float_t         Lepton_track_theta;
   Float_t         Lepton_track_qoverp;
   Float_t         Lepton_track_pt;
   Float_t         Lepton_track_eta;
   Float_t         MissingEt_etx;
   Float_t         MissingEt_ety;
   Float_t         MissingEt_et;
   Float_t         MissingEt_sumet;
   Float_t         MissingEt_phi;
   Int_t           Jet_n;
   std::vector<int>     *Jet_Barcode;
   std::vector<int>     *Jet_flavor_truth_trueflav;
   std::vector<int>     *Jet_flavor_truth_BHadronpdg;
   std::vector<float>   *Jet_pt;
   std::vector<float>   *Jet_eta;
   std::vector<float>   *Jet_phi;
   std::vector<float>   *Jet_E;
   std::vector<float>   *Jet_Et;
   std::vector<float>   *Jet_emscale_pt;
   std::vector<float>   *Jet_emscale_eta;
   std::vector<float>   *Jet_emscale_phi;
   std::vector<float>   *Jet_emscale_E;
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
   Int_t           BJet_n;
   std::vector<int>     *BJet_Barcode;
   std::vector<int>     *BJet_flavor_truth_trueflav;
   std::vector<int>     *BJet_flavor_truth_BHadronpdg;
   std::vector<float>   *BJet_pt;
   std::vector<float>   *BJet_eta;
   std::vector<float>   *BJet_phi;
   std::vector<float>   *BJet_E;
   std::vector<float>   *BJet_Et;
   std::vector<float>   *BJet_emscale_pt;
   std::vector<float>   *BJet_emscale_eta;
   std::vector<float>   *BJet_emscale_phi;
   std::vector<float>   *BJet_emscale_E;
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
   std::vector<int>     *BJet_index;
   UInt_t          RunNumber;
   UInt_t          EventNumber;
   UInt_t          lbn;
   Double_t        weight;
   Bool_t          Lepton_inTrigger;
   Bool_t          isElectronTriggerPassed;
   Int_t           BJet1_Index;
   Int_t           BJet2_Index;
   Int_t           BestJet_Index;
   Int_t           LastJet_Index;
   Int_t           NotBestJet1_Index;
   Int_t           NotBestJet2_Index;
   Int_t           TriangularCut1;
   Int_t           TriangularCut2;
   Int_t           TriangularCut3;
   Int_t           UJet1_Index;
   Int_t           UJet2_Index;
   Int_t           isBestJet;
   Int_t           isNeutrino;
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
   Double_t        Neutrino_mass;
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
   Double_t        W_phi;
   Double_t        W_pt;
   Bool_t          isLeptonTriggerPassed;
   Bool_t          isLeptonTriggerMatched;
   Bool_t          isNoEleMuOverlap;

   // List of branches
   TBranch        *b_Lepton_loose;   //!
   TBranch        *b_Lepton_medium;   //!
   TBranch        *b_Lepton_mediumIso;   //!
   TBranch        *b_Lepton_tight;   //!
   TBranch        *b_Lepton_tightIso;   //!
   TBranch        *b_Lepton_trackBlayer;   //!
   TBranch        *b_Lepton_trackMatchEoverP;   //!
   TBranch        *b_Lepton_trackTRThits;   //!
   TBranch        *b_Lepton_trackTRTratio;   //!
   TBranch        *b_Lepton_trackA0Tight;   //!
   TBranch        *b_Lepton_trackMatchEtaTight;   //!
   TBranch        *b_Lepton_trackMatchPhi;   //!
   TBranch        *b_Lepton_conversionMatch;   //!
   TBranch        *b_Lepton_Barcode;   //!
   TBranch        *b_Lepton_pdgId;   //!
   TBranch        *b_Lepton_charge;   //!
   TBranch        *b_Lepton_Etcone45;   //!
   TBranch        *b_Lepton_Etcone20;   //!
   TBranch        *b_Lepton_Etcone30;   //!
   TBranch        *b_Lepton_Etcone40;   //!
   TBranch        *b_Lepton_ptcone20;   //!
   TBranch        *b_Lepton_ptcone30;   //!
   TBranch        *b_Lepton_ptcone40;   //!
   TBranch        *b_Lepton_cluster_pt;   //!
   TBranch        *b_Lepton_cluster_eta;   //!
   TBranch        *b_Lepton_cluster_phi;   //!
   TBranch        *b_Lepton_cluster_E;   //!
   TBranch        *b_Lepton_cluster_Et;   //!
   TBranch        *b_Lepton_cluster_track_Et;   //!
   TBranch        *b_Lepton_track_d0;   //!
   TBranch        *b_Lepton_track_z0;   //!
   TBranch        *b_Lepton_track_phi;   //!
   TBranch        *b_Lepton_track_theta;   //!
   TBranch        *b_Lepton_track_qoverp;   //!
   TBranch        *b_Lepton_track_pt;   //!
   TBranch        *b_Lepton_track_eta;   //!
   TBranch        *b_MissingEt_etx;   //!
   TBranch        *b_MissingEt_ety;   //!
   TBranch        *b_MissingEt_et;   //!
   TBranch        *b_MissingEt_sumet;   //!
   TBranch        *b_MissingEt_phi;   //!
   TBranch        *b_Jet_n;   //!
   TBranch        *b_Jet_Barcode;   //!
   TBranch        *b_Jet_flavor_truth_trueflav;   //!
   TBranch        *b_Jet_flavor_truth_BHadronpdg;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_E;   //!
   TBranch        *b_Jet_Et;   //!
   TBranch        *b_Jet_emscale_pt;   //!
   TBranch        *b_Jet_emscale_eta;   //!
   TBranch        *b_Jet_emscale_phi;   //!
   TBranch        *b_Jet_emscale_E;   //!
   TBranch        *b_Jet_flavor_weight;   //!
   TBranch        *b_Jet_flavor_weight_TrackCounting2D;   //!
   TBranch        *b_Jet_flavor_weight_JetProb;   //!
   TBranch        *b_Jet_flavor_weight_IP1D;   //!
   TBranch        *b_Jet_flavor_weight_IP2D;   //!
   TBranch        *b_Jet_flavor_weight_IP3D;   //!
   TBranch        *b_Jet_flavor_weight_SV0;   //!
   TBranch        *b_Jet_flavor_weight_SV1;   //!
   TBranch        *b_Jet_flavor_weight_SV2;   //!
   TBranch        *b_Jet_flavor_weight_JetFitterTag;   //!
   TBranch        *b_Jet_flavor_weight_JetFitterCOMB;   //!
   TBranch        *b_Jet_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_Jet_flavor_weight_SoftMuonTag;   //!
   TBranch        *b_Jet_flavor_weight_SoftElectronTag;   //!
   TBranch        *b_Jet_flavor_truth_dRMinToB;   //!
   TBranch        *b_Jet_flavor_truth_dRMinToC;   //!
   TBranch        *b_Jet_flavor_truth_dRMinToT;   //!
   TBranch        *b_BJet_n;   //!
   TBranch        *b_BJet_Barcode;   //!
   TBranch        *b_BJet_flavor_truth_trueflav;   //!
   TBranch        *b_BJet_flavor_truth_BHadronpdg;   //!
   TBranch        *b_BJet_pt;   //!
   TBranch        *b_BJet_eta;   //!
   TBranch        *b_BJet_phi;   //!
   TBranch        *b_BJet_E;   //!
   TBranch        *b_BJet_Et;   //!
   TBranch        *b_BJet_emscale_pt;   //!
   TBranch        *b_BJet_emscale_eta;   //!
   TBranch        *b_BJet_emscale_phi;   //!
   TBranch        *b_BJet_emscale_E;   //!
   TBranch        *b_BJet_flavor_weight;   //!
   TBranch        *b_BJet_flavor_weight_TrackCounting2D;   //!
   TBranch        *b_BJet_flavor_weight_JetProb;   //!
   TBranch        *b_BJet_flavor_weight_IP1D;   //!
   TBranch        *b_BJet_flavor_weight_IP2D;   //!
   TBranch        *b_BJet_flavor_weight_IP3D;   //!
   TBranch        *b_BJet_flavor_weight_SV0;   //!
   TBranch        *b_BJet_flavor_weight_SV1;   //!
   TBranch        *b_BJet_flavor_weight_SV2;   //!
   TBranch        *b_BJet_flavor_weight_JetFitterTag;   //!
   TBranch        *b_BJet_flavor_weight_JetFitterCOMB;   //!
   TBranch        *b_BJet_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_BJet_flavor_weight_SoftMuonTag;   //!
   TBranch        *b_BJet_flavor_weight_SoftElectronTag;   //!
   TBranch        *b_BJet_flavor_truth_dRMinToB;   //!
   TBranch        *b_BJet_flavor_truth_dRMinToC;   //!
   TBranch        *b_BJet_flavor_truth_dRMinToT;   //!
   TBranch        *b_BJet_index;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_Lepton_inTrigger;   //!
   TBranch        *b_isElectronTriggerPassed;   //!
   TBranch        *b_BJet1_Index;   //!
   TBranch        *b_BJet2_Index;   //!
   TBranch        *b_BestJet_Index;   //!
   TBranch        *b_LastJet_Index;   //!
   TBranch        *b_NotBestJet1_Index;   //!
   TBranch        *b_NotBestJet2_Index;   //!
   TBranch        *b_TriangularCut1;   //!
   TBranch        *b_TriangularCut2;   //!
   TBranch        *b_TriangularCut3;   //!
   TBranch        *b_UJet1_Index;   //!
   TBranch        *b_UJet2_Index;   //!
   TBranch        *b_isBestJet;   //!
   TBranch        *b_isNeutrino;   //!
   TBranch        *b_Aplanarity_AllJetsLepton;   //!
   TBranch        *b_Centrality_AllJetsLepton;   //!
   TBranch        *b_Centrality_Jet1Jet2;   //!
   TBranch        *b_Centrality_LastJetLepton;   //!
   TBranch        *b_CosDeltaPhi_BJet1BJet2;   //!
   TBranch        *b_CosDeltaPhi_BJet1Lepton;   //!
   TBranch        *b_CosDeltaPhi_BJet1UJet1;   //!
   TBranch        *b_CosDeltaPhi_BJet1UJet2;   //!
   TBranch        *b_CosDeltaPhi_BJet2Lepton;   //!
   TBranch        *b_CosDeltaPhi_LeptonMissingEt;   //!
   TBranch        *b_CosDeltaPhi_LeptonNeutrino;   //!
   TBranch        *b_CosDeltaPhi_UJet1Lepton;   //!
   TBranch        *b_CosDeltaPhi_UJet1UJet2;   //!
   TBranch        *b_CosDeltaPhi_UJet2Lepton;   //!
   TBranch        *b_Cos_BJet1RAllJets_AllJets;   //!
   TBranch        *b_Cos_Jet1RAllJets_AllJets;   //!
   TBranch        *b_Cos_Jet1RCM_Zaxis;   //!
   TBranch        *b_Cos_LeptonRW_WRBestJetTop;   //!
   TBranch        *b_Cos_Lepton_UJet1_RBJet1Top;   //!
   TBranch        *b_Cos_Lepton_Zaxis_RBJet1Top;   //!
   TBranch        *b_Cos_Lepton_Zaxis_RBestJetTop;   //!
   TBranch        *b_Cos_NotBestJet1RAllJets_AllJets;   //!
   TBranch        *b_Cos_UJet1RAllJets_AllJets;   //!
   TBranch        *b_DeltaEta_BJet1BJet2;   //!
   TBranch        *b_DeltaEta_BJet1Lepton;   //!
   TBranch        *b_DeltaEta_BJet1UJet1;   //!
   TBranch        *b_DeltaEta_BJet1UJet2;   //!
   TBranch        *b_DeltaEta_BJet2Lepton;   //!
   TBranch        *b_DeltaEta_Jet1Lepton;   //!
   TBranch        *b_DeltaEta_Jet2Lepton;   //!
   TBranch        *b_DeltaEta_Jet3Lepton;   //!
   TBranch        *b_DeltaEta_LeptonNeutrino;   //!
   TBranch        *b_DeltaEta_UJet1Lepton;   //!
   TBranch        *b_DeltaEta_UJet1UJet2;   //!
   TBranch        *b_DeltaEta_UJet2Lepton;   //!
   TBranch        *b_DeltaPhi_BJet1BJet2;   //!
   TBranch        *b_DeltaPhi_BJet1Lepton;   //!
   TBranch        *b_DeltaPhi_BJet1UJet1;   //!
   TBranch        *b_DeltaPhi_BJet1UJet2;   //!
   TBranch        *b_DeltaPhi_BJet2Lepton;   //!
   TBranch        *b_DeltaPhi_Jet1Jet2;   //!
   TBranch        *b_DeltaPhi_Jet1Lepton;   //!
   TBranch        *b_DeltaPhi_Jet1MissingEt;   //!
   TBranch        *b_DeltaPhi_Jet1Neutrino;   //!
   TBranch        *b_DeltaPhi_Jet2Lepton;   //!
   TBranch        *b_DeltaPhi_Jet2MissingEt;   //!
   TBranch        *b_DeltaPhi_Jet2Neutrino;   //!
   TBranch        *b_DeltaPhi_Jet3Lepton;   //!
   TBranch        *b_DeltaPhi_Jet3MissingEt;   //!
   TBranch        *b_DeltaPhi_Jet3Neutrino;   //!
   TBranch        *b_DeltaPhi_LeptonMissingEt;   //!
   TBranch        *b_DeltaPhi_LeptonNeutrino;   //!
   TBranch        *b_DeltaPhi_UJet1Lepton;   //!
   TBranch        *b_DeltaPhi_UJet1UJet2;   //!
   TBranch        *b_DeltaPhi_UJet2Lepton;   //!
   TBranch        *b_DeltaPt_BJet1UJet1;   //!
   TBranch        *b_DeltaPt_BJet1UJet2;   //!
   TBranch        *b_DeltaPt_Jet1Jet2;   //!
   TBranch        *b_DeltaRMin_LeptonAllJets;   //!
   TBranch        *b_DeltaR_BJet1BJet2;   //!
   TBranch        *b_DeltaR_BJet1Lepton;   //!
   TBranch        *b_DeltaR_BJet1UJet1;   //!
   TBranch        *b_DeltaR_BJet1UJet2;   //!
   TBranch        *b_DeltaR_BJet2Lepton;   //!
   TBranch        *b_DeltaR_Jet1Jet2;   //!
   TBranch        *b_DeltaR_Jet1Lepton;   //!
   TBranch        *b_DeltaR_Jet2Lepton;   //!
   TBranch        *b_DeltaR_Jet3Lepton;   //!
   TBranch        *b_DeltaR_UJet1Lepton;   //!
   TBranch        *b_DeltaR_UJet1UJet2;   //!
   TBranch        *b_DeltaR_UJet2Lepton;   //!
   TBranch        *b_H_AllJets;   //!
   TBranch        *b_H_AllJetsLeptonNeutrino;   //!
   TBranch        *b_H_AllJetsMinusBJet1;   //!
   TBranch        *b_H_AllJetsMinusBestJet;   //!
   TBranch        *b_H_Jet1Jet2;   //!
   TBranch        *b_H_Jet1Jet2LeptonMissingEt;   //!
   TBranch        *b_H_Jet1Jet2LeptonNeutrino;   //!
   TBranch        *b_H_Jet1Jet2LeptonNeutrinoMinusBJet1;   //!
   TBranch        *b_Ht_AllJets;   //!
   TBranch        *b_Ht_AllJetsLepton;   //!
   TBranch        *b_Ht_AllJetsLeptonMissingEt;   //!
   TBranch        *b_Ht_AllJetsMinusBJet1;   //!
   TBranch        *b_Ht_AllJetsMinusBestJet;   //!
   TBranch        *b_Ht_Jet1Jet2;   //!
   TBranch        *b_Ht_Jet1Jet2LeptonMissingEt;   //!
   TBranch        *b_Ht_Jet1Jet2LeptonNeutrino;   //!
   TBranch        *b_Ht_Jet1Jet2LeptonNeutrinoMinusBJet1;   //!
   TBranch        *b_Ht_LeptonMissingEt;   //!
   TBranch        *b_Ht_LeptonNeutrino;   //!
   TBranch        *b_InvariantMass_AllJets;   //!
   TBranch        *b_InvariantMass_AllJetsLeptonNeutrino;   //!
   TBranch        *b_InvariantMass_AllJetsMinusBJet1;   //!
   TBranch        *b_InvariantMass_AllJetsMinusBestJet;   //!
   TBranch        *b_InvariantMass_BestJetLeptonNeutrino;   //!
   TBranch        *b_Jet1Jet2LeptonNeutrino_InvariantMass;   //!
   TBranch        *b_Jet1Jet2LeptonNeutrino_Pt;   //!
   TBranch        *b_Jet1Jet2_InvariantMass;   //!
   TBranch        *b_Jet1Jet2_TransverseMass;   //!
   TBranch        *b_Jet1Jet2_p;   //!
   TBranch        *b_Jet1Jet2_pt;   //!
   TBranch        *b_Mass_BJet1BJet2;   //!
   TBranch        *b_MaxJetEta;   //!
   TBranch        *b_MinJetEta;   //!
   TBranch        *b_Neutrino_mass;   //!
   TBranch        *b_Neutrino_pz;   //!
   TBranch        *b_P_AllJets;   //!
   TBranch        *b_Pt_AllJetsMinusBJet1;   //!
   TBranch        *b_Pt_AllJetsMinusBestJet;   //!
   TBranch        *b_Pt_BJet1BJet2;   //!
   TBranch        *b_Sphericity_AllJetsLepton;   //!
   TBranch        *b_TopBJet1_eta;   //!
   TBranch        *b_TopBJet1_mass;   //!
   TBranch        *b_TopBJet1_phi;   //!
   TBranch        *b_TopBJet1_pt;   //!
   TBranch        *b_TopBJet2_eta;   //!
   TBranch        *b_TopBJet2_mass;   //!
   TBranch        *b_TopBJet2_phi;   //!
   TBranch        *b_TopBJet2_pt;   //!
   TBranch        *b_TopBestJet_eta;   //!
   TBranch        *b_TopBestJet_mass;   //!
   TBranch        *b_TopBestJet_phi;   //!
   TBranch        *b_TopBestJet_pt;   //!
   TBranch        *b_TopJet1_eta;   //!
   TBranch        *b_TopJet1_mass;   //!
   TBranch        *b_TopJet1_phi;   //!
   TBranch        *b_TopJet1_pt;   //!
   TBranch        *b_TopJet2_eta;   //!
   TBranch        *b_TopJet2_mass;   //!
   TBranch        *b_TopJet2_phi;   //!
   TBranch        *b_TopJet2_pt;   //!
   TBranch        *b_TopJet3_eta;   //!
   TBranch        *b_TopJet3_mass;   //!
   TBranch        *b_TopJet3_phi;   //!
   TBranch        *b_TopJet3_pt;   //!
   TBranch        *b_W_TransverseMass;   //!
   TBranch        *b_W_charge;   //!
   TBranch        *b_W_eta;   //!
   TBranch        *b_W_mass;   //!
   TBranch        *b_W_phi;   //!
   TBranch        *b_W_pt;   //!
   TBranch        *b_isLeptonTriggerPassed;   //!
   TBranch        *b_isLeptonTriggerMatched;   //!
   TBranch        *b_isNoEleMuOverlap;   //!

   ElectronPreTagTree(TTree *tree=0);
   virtual ~ElectronPreTagTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ElectronPreTagTree_cxx
ElectronPreTagTree::ElectronPreTagTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root");
      if (!f) {
         f = new TFile("/msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root");
      }
      tree = (TTree*)gDirectory->Get("ElectronPreTagTree");

   }
   Init(tree);
}

ElectronPreTagTree::~ElectronPreTagTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ElectronPreTagTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ElectronPreTagTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ElectronPreTagTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Jet_Barcode = 0;
   Jet_flavor_truth_trueflav = 0;
   Jet_flavor_truth_BHadronpdg = 0;
   Jet_pt = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_E = 0;
   Jet_Et = 0;
   Jet_emscale_pt = 0;
   Jet_emscale_eta = 0;
   Jet_emscale_phi = 0;
   Jet_emscale_E = 0;
   Jet_flavor_weight = 0;
   Jet_flavor_weight_TrackCounting2D = 0;
   Jet_flavor_weight_JetProb = 0;
   Jet_flavor_weight_IP1D = 0;
   Jet_flavor_weight_IP2D = 0;
   Jet_flavor_weight_IP3D = 0;
   Jet_flavor_weight_SV0 = 0;
   Jet_flavor_weight_SV1 = 0;
   Jet_flavor_weight_SV2 = 0;
   Jet_flavor_weight_JetFitterTag = 0;
   Jet_flavor_weight_JetFitterCOMB = 0;
   Jet_flavor_weight_JetFitterCOMBNN = 0;
   Jet_flavor_weight_SoftMuonTag = 0;
   Jet_flavor_weight_SoftElectronTag = 0;
   Jet_flavor_truth_dRMinToB = 0;
   Jet_flavor_truth_dRMinToC = 0;
   Jet_flavor_truth_dRMinToT = 0;
   BJet_Barcode = 0;
   BJet_flavor_truth_trueflav = 0;
   BJet_flavor_truth_BHadronpdg = 0;
   BJet_pt = 0;
   BJet_eta = 0;
   BJet_phi = 0;
   BJet_E = 0;
   BJet_Et = 0;
   BJet_emscale_pt = 0;
   BJet_emscale_eta = 0;
   BJet_emscale_phi = 0;
   BJet_emscale_E = 0;
   BJet_flavor_weight = 0;
   BJet_flavor_weight_TrackCounting2D = 0;
   BJet_flavor_weight_JetProb = 0;
   BJet_flavor_weight_IP1D = 0;
   BJet_flavor_weight_IP2D = 0;
   BJet_flavor_weight_IP3D = 0;
   BJet_flavor_weight_SV0 = 0;
   BJet_flavor_weight_SV1 = 0;
   BJet_flavor_weight_SV2 = 0;
   BJet_flavor_weight_JetFitterTag = 0;
   BJet_flavor_weight_JetFitterCOMB = 0;
   BJet_flavor_weight_JetFitterCOMBNN = 0;
   BJet_flavor_weight_SoftMuonTag = 0;
   BJet_flavor_weight_SoftElectronTag = 0;
   BJet_flavor_truth_dRMinToB = 0;
   BJet_flavor_truth_dRMinToC = 0;
   BJet_flavor_truth_dRMinToT = 0;
   BJet_index = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Lepton_loose", &Lepton_loose, &b_Lepton_loose);
   fChain->SetBranchAddress("Lepton_medium", &Lepton_medium, &b_Lepton_medium);
   fChain->SetBranchAddress("Lepton_mediumIso", &Lepton_mediumIso, &b_Lepton_mediumIso);
   fChain->SetBranchAddress("Lepton_tight", &Lepton_tight, &b_Lepton_tight);
   fChain->SetBranchAddress("Lepton_tightIso", &Lepton_tightIso, &b_Lepton_tightIso);
   fChain->SetBranchAddress("Lepton_trackBlayer", &Lepton_trackBlayer, &b_Lepton_trackBlayer);
   fChain->SetBranchAddress("Lepton_trackMatchEoverP", &Lepton_trackMatchEoverP, &b_Lepton_trackMatchEoverP);
   fChain->SetBranchAddress("Lepton_trackTRThits", &Lepton_trackTRThits, &b_Lepton_trackTRThits);
   fChain->SetBranchAddress("Lepton_trackTRTratio", &Lepton_trackTRTratio, &b_Lepton_trackTRTratio);
   fChain->SetBranchAddress("Lepton_trackA0Tight", &Lepton_trackA0Tight, &b_Lepton_trackA0Tight);
   fChain->SetBranchAddress("Lepton_trackMatchEtaTight", &Lepton_trackMatchEtaTight, &b_Lepton_trackMatchEtaTight);
   fChain->SetBranchAddress("Lepton_trackMatchPhi", &Lepton_trackMatchPhi, &b_Lepton_trackMatchPhi);
   fChain->SetBranchAddress("Lepton_conversionMatch", &Lepton_conversionMatch, &b_Lepton_conversionMatch);
   fChain->SetBranchAddress("Lepton_Barcode", &Lepton_Barcode, &b_Lepton_Barcode);
   fChain->SetBranchAddress("Lepton_pdgId", &Lepton_pdgId, &b_Lepton_pdgId);
   fChain->SetBranchAddress("Lepton_charge", &Lepton_charge, &b_Lepton_charge);
   fChain->SetBranchAddress("Lepton_Etcone45", &Lepton_Etcone45, &b_Lepton_Etcone45);
   fChain->SetBranchAddress("Lepton_Etcone20", &Lepton_Etcone20, &b_Lepton_Etcone20);
   fChain->SetBranchAddress("Lepton_Etcone30", &Lepton_Etcone30, &b_Lepton_Etcone30);
   fChain->SetBranchAddress("Lepton_Etcone40", &Lepton_Etcone40, &b_Lepton_Etcone40);
   fChain->SetBranchAddress("Lepton_ptcone20", &Lepton_ptcone20, &b_Lepton_ptcone20);
   fChain->SetBranchAddress("Lepton_ptcone30", &Lepton_ptcone30, &b_Lepton_ptcone30);
   fChain->SetBranchAddress("Lepton_ptcone40", &Lepton_ptcone40, &b_Lepton_ptcone40);
   fChain->SetBranchAddress("Lepton_cluster_pt", &Lepton_cluster_pt, &b_Lepton_cluster_pt);
   fChain->SetBranchAddress("Lepton_cluster_eta", &Lepton_cluster_eta, &b_Lepton_cluster_eta);
   fChain->SetBranchAddress("Lepton_cluster_phi", &Lepton_cluster_phi, &b_Lepton_cluster_phi);
   fChain->SetBranchAddress("Lepton_cluster_E", &Lepton_cluster_E, &b_Lepton_cluster_E);
   fChain->SetBranchAddress("Lepton_cluster_Et", &Lepton_cluster_Et, &b_Lepton_cluster_Et);
   fChain->SetBranchAddress("Lepton_cluster_track_Et", &Lepton_cluster_track_Et, &b_Lepton_cluster_track_Et);
   fChain->SetBranchAddress("Lepton_track_d0", &Lepton_track_d0, &b_Lepton_track_d0);
   fChain->SetBranchAddress("Lepton_track_z0", &Lepton_track_z0, &b_Lepton_track_z0);
   fChain->SetBranchAddress("Lepton_track_phi", &Lepton_track_phi, &b_Lepton_track_phi);
   fChain->SetBranchAddress("Lepton_track_theta", &Lepton_track_theta, &b_Lepton_track_theta);
   fChain->SetBranchAddress("Lepton_track_qoverp", &Lepton_track_qoverp, &b_Lepton_track_qoverp);
   fChain->SetBranchAddress("Lepton_track_pt", &Lepton_track_pt, &b_Lepton_track_pt);
   fChain->SetBranchAddress("Lepton_track_eta", &Lepton_track_eta, &b_Lepton_track_eta);
   fChain->SetBranchAddress("MissingEt_etx", &MissingEt_etx, &b_MissingEt_etx);
   fChain->SetBranchAddress("MissingEt_ety", &MissingEt_ety, &b_MissingEt_ety);
   fChain->SetBranchAddress("MissingEt_et", &MissingEt_et, &b_MissingEt_et);
   fChain->SetBranchAddress("MissingEt_sumet", &MissingEt_sumet, &b_MissingEt_sumet);
   fChain->SetBranchAddress("MissingEt_phi", &MissingEt_phi, &b_MissingEt_phi);
   fChain->SetBranchAddress("Jet_n", &Jet_n, &b_Jet_n);
   fChain->SetBranchAddress("Jet_Barcode", &Jet_Barcode, &b_Jet_Barcode);
   fChain->SetBranchAddress("Jet_flavor_truth_trueflav", &Jet_flavor_truth_trueflav, &b_Jet_flavor_truth_trueflav);
   fChain->SetBranchAddress("Jet_flavor_truth_BHadronpdg", &Jet_flavor_truth_BHadronpdg, &b_Jet_flavor_truth_BHadronpdg);
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_E", &Jet_E, &b_Jet_E);
   fChain->SetBranchAddress("Jet_Et", &Jet_Et, &b_Jet_Et);
   fChain->SetBranchAddress("Jet_emscale_pt", &Jet_emscale_pt, &b_Jet_emscale_pt);
   fChain->SetBranchAddress("Jet_emscale_eta", &Jet_emscale_eta, &b_Jet_emscale_eta);
   fChain->SetBranchAddress("Jet_emscale_phi", &Jet_emscale_phi, &b_Jet_emscale_phi);
   fChain->SetBranchAddress("Jet_emscale_E", &Jet_emscale_E, &b_Jet_emscale_E);
   fChain->SetBranchAddress("Jet_flavor_weight", &Jet_flavor_weight, &b_Jet_flavor_weight);
   fChain->SetBranchAddress("Jet_flavor_weight_TrackCounting2D", &Jet_flavor_weight_TrackCounting2D, &b_Jet_flavor_weight_TrackCounting2D);
   fChain->SetBranchAddress("Jet_flavor_weight_JetProb", &Jet_flavor_weight_JetProb, &b_Jet_flavor_weight_JetProb);
   fChain->SetBranchAddress("Jet_flavor_weight_IP1D", &Jet_flavor_weight_IP1D, &b_Jet_flavor_weight_IP1D);
   fChain->SetBranchAddress("Jet_flavor_weight_IP2D", &Jet_flavor_weight_IP2D, &b_Jet_flavor_weight_IP2D);
   fChain->SetBranchAddress("Jet_flavor_weight_IP3D", &Jet_flavor_weight_IP3D, &b_Jet_flavor_weight_IP3D);
   fChain->SetBranchAddress("Jet_flavor_weight_SV0", &Jet_flavor_weight_SV0, &b_Jet_flavor_weight_SV0);
   fChain->SetBranchAddress("Jet_flavor_weight_SV1", &Jet_flavor_weight_SV1, &b_Jet_flavor_weight_SV1);
   fChain->SetBranchAddress("Jet_flavor_weight_SV2", &Jet_flavor_weight_SV2, &b_Jet_flavor_weight_SV2);
   fChain->SetBranchAddress("Jet_flavor_weight_JetFitterTag", &Jet_flavor_weight_JetFitterTag, &b_Jet_flavor_weight_JetFitterTag);
   fChain->SetBranchAddress("Jet_flavor_weight_JetFitterCOMB", &Jet_flavor_weight_JetFitterCOMB, &b_Jet_flavor_weight_JetFitterCOMB);
   fChain->SetBranchAddress("Jet_flavor_weight_JetFitterCOMBNN", &Jet_flavor_weight_JetFitterCOMBNN, &b_Jet_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("Jet_flavor_weight_SoftMuonTag", &Jet_flavor_weight_SoftMuonTag, &b_Jet_flavor_weight_SoftMuonTag);
   fChain->SetBranchAddress("Jet_flavor_weight_SoftElectronTag", &Jet_flavor_weight_SoftElectronTag, &b_Jet_flavor_weight_SoftElectronTag);
   fChain->SetBranchAddress("Jet_flavor_truth_dRMinToB", &Jet_flavor_truth_dRMinToB, &b_Jet_flavor_truth_dRMinToB);
   fChain->SetBranchAddress("Jet_flavor_truth_dRMinToC", &Jet_flavor_truth_dRMinToC, &b_Jet_flavor_truth_dRMinToC);
   fChain->SetBranchAddress("Jet_flavor_truth_dRMinToT", &Jet_flavor_truth_dRMinToT, &b_Jet_flavor_truth_dRMinToT);
   fChain->SetBranchAddress("BJet_n", &BJet_n, &b_BJet_n);
   fChain->SetBranchAddress("BJet_Barcode", &BJet_Barcode, &b_BJet_Barcode);
   fChain->SetBranchAddress("BJet_flavor_truth_trueflav", &BJet_flavor_truth_trueflav, &b_BJet_flavor_truth_trueflav);
   fChain->SetBranchAddress("BJet_flavor_truth_BHadronpdg", &BJet_flavor_truth_BHadronpdg, &b_BJet_flavor_truth_BHadronpdg);
   fChain->SetBranchAddress("BJet_pt", &BJet_pt, &b_BJet_pt);
   fChain->SetBranchAddress("BJet_eta", &BJet_eta, &b_BJet_eta);
   fChain->SetBranchAddress("BJet_phi", &BJet_phi, &b_BJet_phi);
   fChain->SetBranchAddress("BJet_E", &BJet_E, &b_BJet_E);
   fChain->SetBranchAddress("BJet_Et", &BJet_Et, &b_BJet_Et);
   fChain->SetBranchAddress("BJet_emscale_pt", &BJet_emscale_pt, &b_BJet_emscale_pt);
   fChain->SetBranchAddress("BJet_emscale_eta", &BJet_emscale_eta, &b_BJet_emscale_eta);
   fChain->SetBranchAddress("BJet_emscale_phi", &BJet_emscale_phi, &b_BJet_emscale_phi);
   fChain->SetBranchAddress("BJet_emscale_E", &BJet_emscale_E, &b_BJet_emscale_E);
   fChain->SetBranchAddress("BJet_flavor_weight", &BJet_flavor_weight, &b_BJet_flavor_weight);
   fChain->SetBranchAddress("BJet_flavor_weight_TrackCounting2D", &BJet_flavor_weight_TrackCounting2D, &b_BJet_flavor_weight_TrackCounting2D);
   fChain->SetBranchAddress("BJet_flavor_weight_JetProb", &BJet_flavor_weight_JetProb, &b_BJet_flavor_weight_JetProb);
   fChain->SetBranchAddress("BJet_flavor_weight_IP1D", &BJet_flavor_weight_IP1D, &b_BJet_flavor_weight_IP1D);
   fChain->SetBranchAddress("BJet_flavor_weight_IP2D", &BJet_flavor_weight_IP2D, &b_BJet_flavor_weight_IP2D);
   fChain->SetBranchAddress("BJet_flavor_weight_IP3D", &BJet_flavor_weight_IP3D, &b_BJet_flavor_weight_IP3D);
   fChain->SetBranchAddress("BJet_flavor_weight_SV0", &BJet_flavor_weight_SV0, &b_BJet_flavor_weight_SV0);
   fChain->SetBranchAddress("BJet_flavor_weight_SV1", &BJet_flavor_weight_SV1, &b_BJet_flavor_weight_SV1);
   fChain->SetBranchAddress("BJet_flavor_weight_SV2", &BJet_flavor_weight_SV2, &b_BJet_flavor_weight_SV2);
   fChain->SetBranchAddress("BJet_flavor_weight_JetFitterTag", &BJet_flavor_weight_JetFitterTag, &b_BJet_flavor_weight_JetFitterTag);
   fChain->SetBranchAddress("BJet_flavor_weight_JetFitterCOMB", &BJet_flavor_weight_JetFitterCOMB, &b_BJet_flavor_weight_JetFitterCOMB);
   fChain->SetBranchAddress("BJet_flavor_weight_JetFitterCOMBNN", &BJet_flavor_weight_JetFitterCOMBNN, &b_BJet_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("BJet_flavor_weight_SoftMuonTag", &BJet_flavor_weight_SoftMuonTag, &b_BJet_flavor_weight_SoftMuonTag);
   fChain->SetBranchAddress("BJet_flavor_weight_SoftElectronTag", &BJet_flavor_weight_SoftElectronTag, &b_BJet_flavor_weight_SoftElectronTag);
   fChain->SetBranchAddress("BJet_flavor_truth_dRMinToB", &BJet_flavor_truth_dRMinToB, &b_BJet_flavor_truth_dRMinToB);
   fChain->SetBranchAddress("BJet_flavor_truth_dRMinToC", &BJet_flavor_truth_dRMinToC, &b_BJet_flavor_truth_dRMinToC);
   fChain->SetBranchAddress("BJet_flavor_truth_dRMinToT", &BJet_flavor_truth_dRMinToT, &b_BJet_flavor_truth_dRMinToT);
   fChain->SetBranchAddress("BJet_index", &BJet_index, &b_BJet_index);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("Lepton_inTrigger", &Lepton_inTrigger, &b_Lepton_inTrigger);
   fChain->SetBranchAddress("isElectronTriggerPassed", &isElectronTriggerPassed, &b_isElectronTriggerPassed);
   fChain->SetBranchAddress("BJet1_Index", &BJet1_Index, &b_BJet1_Index);
   fChain->SetBranchAddress("BJet2_Index", &BJet2_Index, &b_BJet2_Index);
   fChain->SetBranchAddress("BestJet_Index", &BestJet_Index, &b_BestJet_Index);
   fChain->SetBranchAddress("LastJet_Index", &LastJet_Index, &b_LastJet_Index);
   fChain->SetBranchAddress("NotBestJet1_Index", &NotBestJet1_Index, &b_NotBestJet1_Index);
   fChain->SetBranchAddress("NotBestJet2_Index", &NotBestJet2_Index, &b_NotBestJet2_Index);
   fChain->SetBranchAddress("TriangularCut1", &TriangularCut1, &b_TriangularCut1);
   fChain->SetBranchAddress("TriangularCut2", &TriangularCut2, &b_TriangularCut2);
   fChain->SetBranchAddress("TriangularCut3", &TriangularCut3, &b_TriangularCut3);
   fChain->SetBranchAddress("UJet1_Index", &UJet1_Index, &b_UJet1_Index);
   fChain->SetBranchAddress("UJet2_Index", &UJet2_Index, &b_UJet2_Index);
   fChain->SetBranchAddress("isBestJet", &isBestJet, &b_isBestJet);
   fChain->SetBranchAddress("isNeutrino", &isNeutrino, &b_isNeutrino);
   fChain->SetBranchAddress("Aplanarity_AllJetsLepton", &Aplanarity_AllJetsLepton, &b_Aplanarity_AllJetsLepton);
   fChain->SetBranchAddress("Centrality_AllJetsLepton", &Centrality_AllJetsLepton, &b_Centrality_AllJetsLepton);
   fChain->SetBranchAddress("Centrality_Jet1Jet2", &Centrality_Jet1Jet2, &b_Centrality_Jet1Jet2);
   fChain->SetBranchAddress("Centrality_LastJetLepton", &Centrality_LastJetLepton, &b_Centrality_LastJetLepton);
   fChain->SetBranchAddress("CosDeltaPhi_BJet1BJet2", &CosDeltaPhi_BJet1BJet2, &b_CosDeltaPhi_BJet1BJet2);
   fChain->SetBranchAddress("CosDeltaPhi_BJet1Lepton", &CosDeltaPhi_BJet1Lepton, &b_CosDeltaPhi_BJet1Lepton);
   fChain->SetBranchAddress("CosDeltaPhi_BJet1UJet1", &CosDeltaPhi_BJet1UJet1, &b_CosDeltaPhi_BJet1UJet1);
   fChain->SetBranchAddress("CosDeltaPhi_BJet1UJet2", &CosDeltaPhi_BJet1UJet2, &b_CosDeltaPhi_BJet1UJet2);
   fChain->SetBranchAddress("CosDeltaPhi_BJet2Lepton", &CosDeltaPhi_BJet2Lepton, &b_CosDeltaPhi_BJet2Lepton);
   fChain->SetBranchAddress("CosDeltaPhi_LeptonMissingEt", &CosDeltaPhi_LeptonMissingEt, &b_CosDeltaPhi_LeptonMissingEt);
   fChain->SetBranchAddress("CosDeltaPhi_LeptonNeutrino", &CosDeltaPhi_LeptonNeutrino, &b_CosDeltaPhi_LeptonNeutrino);
   fChain->SetBranchAddress("CosDeltaPhi_UJet1Lepton", &CosDeltaPhi_UJet1Lepton, &b_CosDeltaPhi_UJet1Lepton);
   fChain->SetBranchAddress("CosDeltaPhi_UJet1UJet2", &CosDeltaPhi_UJet1UJet2, &b_CosDeltaPhi_UJet1UJet2);
   fChain->SetBranchAddress("CosDeltaPhi_UJet2Lepton", &CosDeltaPhi_UJet2Lepton, &b_CosDeltaPhi_UJet2Lepton);
   fChain->SetBranchAddress("Cos_BJet1RAllJets_AllJets", &Cos_BJet1RAllJets_AllJets, &b_Cos_BJet1RAllJets_AllJets);
   fChain->SetBranchAddress("Cos_Jet1RAllJets_AllJets", &Cos_Jet1RAllJets_AllJets, &b_Cos_Jet1RAllJets_AllJets);
   fChain->SetBranchAddress("Cos_Jet1RCM_Zaxis", &Cos_Jet1RCM_Zaxis, &b_Cos_Jet1RCM_Zaxis);
   fChain->SetBranchAddress("Cos_LeptonRW_WRBestJetTop", &Cos_LeptonRW_WRBestJetTop, &b_Cos_LeptonRW_WRBestJetTop);
   fChain->SetBranchAddress("Cos_Lepton_UJet1_RBJet1Top", &Cos_Lepton_UJet1_RBJet1Top, &b_Cos_Lepton_UJet1_RBJet1Top);
   fChain->SetBranchAddress("Cos_Lepton_Zaxis_RBJet1Top", &Cos_Lepton_Zaxis_RBJet1Top, &b_Cos_Lepton_Zaxis_RBJet1Top);
   fChain->SetBranchAddress("Cos_Lepton_Zaxis_RBestJetTop", &Cos_Lepton_Zaxis_RBestJetTop, &b_Cos_Lepton_Zaxis_RBestJetTop);
   fChain->SetBranchAddress("Cos_NotBestJet1RAllJets_AllJets", &Cos_NotBestJet1RAllJets_AllJets, &b_Cos_NotBestJet1RAllJets_AllJets);
   fChain->SetBranchAddress("Cos_UJet1RAllJets_AllJets", &Cos_UJet1RAllJets_AllJets, &b_Cos_UJet1RAllJets_AllJets);
   fChain->SetBranchAddress("DeltaEta_BJet1BJet2", &DeltaEta_BJet1BJet2, &b_DeltaEta_BJet1BJet2);
   fChain->SetBranchAddress("DeltaEta_BJet1Lepton", &DeltaEta_BJet1Lepton, &b_DeltaEta_BJet1Lepton);
   fChain->SetBranchAddress("DeltaEta_BJet1UJet1", &DeltaEta_BJet1UJet1, &b_DeltaEta_BJet1UJet1);
   fChain->SetBranchAddress("DeltaEta_BJet1UJet2", &DeltaEta_BJet1UJet2, &b_DeltaEta_BJet1UJet2);
   fChain->SetBranchAddress("DeltaEta_BJet2Lepton", &DeltaEta_BJet2Lepton, &b_DeltaEta_BJet2Lepton);
   fChain->SetBranchAddress("DeltaEta_Jet1Lepton", &DeltaEta_Jet1Lepton, &b_DeltaEta_Jet1Lepton);
   fChain->SetBranchAddress("DeltaEta_Jet2Lepton", &DeltaEta_Jet2Lepton, &b_DeltaEta_Jet2Lepton);
   fChain->SetBranchAddress("DeltaEta_Jet3Lepton", &DeltaEta_Jet3Lepton, &b_DeltaEta_Jet3Lepton);
   fChain->SetBranchAddress("DeltaEta_LeptonNeutrino", &DeltaEta_LeptonNeutrino, &b_DeltaEta_LeptonNeutrino);
   fChain->SetBranchAddress("DeltaEta_UJet1Lepton", &DeltaEta_UJet1Lepton, &b_DeltaEta_UJet1Lepton);
   fChain->SetBranchAddress("DeltaEta_UJet1UJet2", &DeltaEta_UJet1UJet2, &b_DeltaEta_UJet1UJet2);
   fChain->SetBranchAddress("DeltaEta_UJet2Lepton", &DeltaEta_UJet2Lepton, &b_DeltaEta_UJet2Lepton);
   fChain->SetBranchAddress("DeltaPhi_BJet1BJet2", &DeltaPhi_BJet1BJet2, &b_DeltaPhi_BJet1BJet2);
   fChain->SetBranchAddress("DeltaPhi_BJet1Lepton", &DeltaPhi_BJet1Lepton, &b_DeltaPhi_BJet1Lepton);
   fChain->SetBranchAddress("DeltaPhi_BJet1UJet1", &DeltaPhi_BJet1UJet1, &b_DeltaPhi_BJet1UJet1);
   fChain->SetBranchAddress("DeltaPhi_BJet1UJet2", &DeltaPhi_BJet1UJet2, &b_DeltaPhi_BJet1UJet2);
   fChain->SetBranchAddress("DeltaPhi_BJet2Lepton", &DeltaPhi_BJet2Lepton, &b_DeltaPhi_BJet2Lepton);
   fChain->SetBranchAddress("DeltaPhi_Jet1Jet2", &DeltaPhi_Jet1Jet2, &b_DeltaPhi_Jet1Jet2);
   fChain->SetBranchAddress("DeltaPhi_Jet1Lepton", &DeltaPhi_Jet1Lepton, &b_DeltaPhi_Jet1Lepton);
   fChain->SetBranchAddress("DeltaPhi_Jet1MissingEt", &DeltaPhi_Jet1MissingEt, &b_DeltaPhi_Jet1MissingEt);
   fChain->SetBranchAddress("DeltaPhi_Jet1Neutrino", &DeltaPhi_Jet1Neutrino, &b_DeltaPhi_Jet1Neutrino);
   fChain->SetBranchAddress("DeltaPhi_Jet2Lepton", &DeltaPhi_Jet2Lepton, &b_DeltaPhi_Jet2Lepton);
   fChain->SetBranchAddress("DeltaPhi_Jet2MissingEt", &DeltaPhi_Jet2MissingEt, &b_DeltaPhi_Jet2MissingEt);
   fChain->SetBranchAddress("DeltaPhi_Jet2Neutrino", &DeltaPhi_Jet2Neutrino, &b_DeltaPhi_Jet2Neutrino);
   fChain->SetBranchAddress("DeltaPhi_Jet3Lepton", &DeltaPhi_Jet3Lepton, &b_DeltaPhi_Jet3Lepton);
   fChain->SetBranchAddress("DeltaPhi_Jet3MissingEt", &DeltaPhi_Jet3MissingEt, &b_DeltaPhi_Jet3MissingEt);
   fChain->SetBranchAddress("DeltaPhi_Jet3Neutrino", &DeltaPhi_Jet3Neutrino, &b_DeltaPhi_Jet3Neutrino);
   fChain->SetBranchAddress("DeltaPhi_LeptonMissingEt", &DeltaPhi_LeptonMissingEt, &b_DeltaPhi_LeptonMissingEt);
   fChain->SetBranchAddress("DeltaPhi_LeptonNeutrino", &DeltaPhi_LeptonNeutrino, &b_DeltaPhi_LeptonNeutrino);
   fChain->SetBranchAddress("DeltaPhi_UJet1Lepton", &DeltaPhi_UJet1Lepton, &b_DeltaPhi_UJet1Lepton);
   fChain->SetBranchAddress("DeltaPhi_UJet1UJet2", &DeltaPhi_UJet1UJet2, &b_DeltaPhi_UJet1UJet2);
   fChain->SetBranchAddress("DeltaPhi_UJet2Lepton", &DeltaPhi_UJet2Lepton, &b_DeltaPhi_UJet2Lepton);
   fChain->SetBranchAddress("DeltaPt_BJet1UJet1", &DeltaPt_BJet1UJet1, &b_DeltaPt_BJet1UJet1);
   fChain->SetBranchAddress("DeltaPt_BJet1UJet2", &DeltaPt_BJet1UJet2, &b_DeltaPt_BJet1UJet2);
   fChain->SetBranchAddress("DeltaPt_Jet1Jet2", &DeltaPt_Jet1Jet2, &b_DeltaPt_Jet1Jet2);
   fChain->SetBranchAddress("DeltaRMin_LeptonAllJets", &DeltaRMin_LeptonAllJets, &b_DeltaRMin_LeptonAllJets);
   fChain->SetBranchAddress("DeltaR_BJet1BJet2", &DeltaR_BJet1BJet2, &b_DeltaR_BJet1BJet2);
   fChain->SetBranchAddress("DeltaR_BJet1Lepton", &DeltaR_BJet1Lepton, &b_DeltaR_BJet1Lepton);
   fChain->SetBranchAddress("DeltaR_BJet1UJet1", &DeltaR_BJet1UJet1, &b_DeltaR_BJet1UJet1);
   fChain->SetBranchAddress("DeltaR_BJet1UJet2", &DeltaR_BJet1UJet2, &b_DeltaR_BJet1UJet2);
   fChain->SetBranchAddress("DeltaR_BJet2Lepton", &DeltaR_BJet2Lepton, &b_DeltaR_BJet2Lepton);
   fChain->SetBranchAddress("DeltaR_Jet1Jet2", &DeltaR_Jet1Jet2, &b_DeltaR_Jet1Jet2);
   fChain->SetBranchAddress("DeltaR_Jet1Lepton", &DeltaR_Jet1Lepton, &b_DeltaR_Jet1Lepton);
   fChain->SetBranchAddress("DeltaR_Jet2Lepton", &DeltaR_Jet2Lepton, &b_DeltaR_Jet2Lepton);
   fChain->SetBranchAddress("DeltaR_Jet3Lepton", &DeltaR_Jet3Lepton, &b_DeltaR_Jet3Lepton);
   fChain->SetBranchAddress("DeltaR_UJet1Lepton", &DeltaR_UJet1Lepton, &b_DeltaR_UJet1Lepton);
   fChain->SetBranchAddress("DeltaR_UJet1UJet2", &DeltaR_UJet1UJet2, &b_DeltaR_UJet1UJet2);
   fChain->SetBranchAddress("DeltaR_UJet2Lepton", &DeltaR_UJet2Lepton, &b_DeltaR_UJet2Lepton);
   fChain->SetBranchAddress("H_AllJets", &H_AllJets, &b_H_AllJets);
   fChain->SetBranchAddress("H_AllJetsLeptonNeutrino", &H_AllJetsLeptonNeutrino, &b_H_AllJetsLeptonNeutrino);
   fChain->SetBranchAddress("H_AllJetsMinusBJet1", &H_AllJetsMinusBJet1, &b_H_AllJetsMinusBJet1);
   fChain->SetBranchAddress("H_AllJetsMinusBestJet", &H_AllJetsMinusBestJet, &b_H_AllJetsMinusBestJet);
   fChain->SetBranchAddress("H_Jet1Jet2", &H_Jet1Jet2, &b_H_Jet1Jet2);
   fChain->SetBranchAddress("H_Jet1Jet2LeptonMissingEt", &H_Jet1Jet2LeptonMissingEt, &b_H_Jet1Jet2LeptonMissingEt);
   fChain->SetBranchAddress("H_Jet1Jet2LeptonNeutrino", &H_Jet1Jet2LeptonNeutrino, &b_H_Jet1Jet2LeptonNeutrino);
   fChain->SetBranchAddress("H_Jet1Jet2LeptonNeutrinoMinusBJet1", &H_Jet1Jet2LeptonNeutrinoMinusBJet1, &b_H_Jet1Jet2LeptonNeutrinoMinusBJet1);
   fChain->SetBranchAddress("Ht_AllJets", &Ht_AllJets, &b_Ht_AllJets);
   fChain->SetBranchAddress("Ht_AllJetsLepton", &Ht_AllJetsLepton, &b_Ht_AllJetsLepton);
   fChain->SetBranchAddress("Ht_AllJetsLeptonMissingEt", &Ht_AllJetsLeptonMissingEt, &b_Ht_AllJetsLeptonMissingEt);
   fChain->SetBranchAddress("Ht_AllJetsMinusBJet1", &Ht_AllJetsMinusBJet1, &b_Ht_AllJetsMinusBJet1);
   fChain->SetBranchAddress("Ht_AllJetsMinusBestJet", &Ht_AllJetsMinusBestJet, &b_Ht_AllJetsMinusBestJet);
   fChain->SetBranchAddress("Ht_Jet1Jet2", &Ht_Jet1Jet2, &b_Ht_Jet1Jet2);
   fChain->SetBranchAddress("Ht_Jet1Jet2LeptonMissingEt", &Ht_Jet1Jet2LeptonMissingEt, &b_Ht_Jet1Jet2LeptonMissingEt);
   fChain->SetBranchAddress("Ht_Jet1Jet2LeptonNeutrino", &Ht_Jet1Jet2LeptonNeutrino, &b_Ht_Jet1Jet2LeptonNeutrino);
   fChain->SetBranchAddress("Ht_Jet1Jet2LeptonNeutrinoMinusBJet1", &Ht_Jet1Jet2LeptonNeutrinoMinusBJet1, &b_Ht_Jet1Jet2LeptonNeutrinoMinusBJet1);
   fChain->SetBranchAddress("Ht_LeptonMissingEt", &Ht_LeptonMissingEt, &b_Ht_LeptonMissingEt);
   fChain->SetBranchAddress("Ht_LeptonNeutrino", &Ht_LeptonNeutrino, &b_Ht_LeptonNeutrino);
   fChain->SetBranchAddress("InvariantMass_AllJets", &InvariantMass_AllJets, &b_InvariantMass_AllJets);
   fChain->SetBranchAddress("InvariantMass_AllJetsLeptonNeutrino", &InvariantMass_AllJetsLeptonNeutrino, &b_InvariantMass_AllJetsLeptonNeutrino);
   fChain->SetBranchAddress("InvariantMass_AllJetsMinusBJet1", &InvariantMass_AllJetsMinusBJet1, &b_InvariantMass_AllJetsMinusBJet1);
   fChain->SetBranchAddress("InvariantMass_AllJetsMinusBestJet", &InvariantMass_AllJetsMinusBestJet, &b_InvariantMass_AllJetsMinusBestJet);
   fChain->SetBranchAddress("InvariantMass_BestJetLeptonNeutrino", &InvariantMass_BestJetLeptonNeutrino, &b_InvariantMass_BestJetLeptonNeutrino);
   fChain->SetBranchAddress("Jet1Jet2LeptonNeutrino_InvariantMass", &Jet1Jet2LeptonNeutrino_InvariantMass, &b_Jet1Jet2LeptonNeutrino_InvariantMass);
   fChain->SetBranchAddress("Jet1Jet2LeptonNeutrino_Pt", &Jet1Jet2LeptonNeutrino_Pt, &b_Jet1Jet2LeptonNeutrino_Pt);
   fChain->SetBranchAddress("Jet1Jet2_InvariantMass", &Jet1Jet2_InvariantMass, &b_Jet1Jet2_InvariantMass);
   fChain->SetBranchAddress("Jet1Jet2_TransverseMass", &Jet1Jet2_TransverseMass, &b_Jet1Jet2_TransverseMass);
   fChain->SetBranchAddress("Jet1Jet2_p", &Jet1Jet2_p, &b_Jet1Jet2_p);
   fChain->SetBranchAddress("Jet1Jet2_pt", &Jet1Jet2_pt, &b_Jet1Jet2_pt);
   fChain->SetBranchAddress("Mass_BJet1BJet2", &Mass_BJet1BJet2, &b_Mass_BJet1BJet2);
   fChain->SetBranchAddress("MaxJetEta", &MaxJetEta, &b_MaxJetEta);
   fChain->SetBranchAddress("MinJetEta", &MinJetEta, &b_MinJetEta);
   fChain->SetBranchAddress("Neutrino_mass", &Neutrino_mass, &b_Neutrino_mass);
   fChain->SetBranchAddress("Neutrino_pz", &Neutrino_pz, &b_Neutrino_pz);
   fChain->SetBranchAddress("P_AllJets", &P_AllJets, &b_P_AllJets);
   fChain->SetBranchAddress("Pt_AllJetsMinusBJet1", &Pt_AllJetsMinusBJet1, &b_Pt_AllJetsMinusBJet1);
   fChain->SetBranchAddress("Pt_AllJetsMinusBestJet", &Pt_AllJetsMinusBestJet, &b_Pt_AllJetsMinusBestJet);
   fChain->SetBranchAddress("Pt_BJet1BJet2", &Pt_BJet1BJet2, &b_Pt_BJet1BJet2);
   fChain->SetBranchAddress("Sphericity_AllJetsLepton", &Sphericity_AllJetsLepton, &b_Sphericity_AllJetsLepton);
   fChain->SetBranchAddress("TopBJet1_eta", &TopBJet1_eta, &b_TopBJet1_eta);
   fChain->SetBranchAddress("TopBJet1_mass", &TopBJet1_mass, &b_TopBJet1_mass);
   fChain->SetBranchAddress("TopBJet1_phi", &TopBJet1_phi, &b_TopBJet1_phi);
   fChain->SetBranchAddress("TopBJet1_pt", &TopBJet1_pt, &b_TopBJet1_pt);
   fChain->SetBranchAddress("TopBJet2_eta", &TopBJet2_eta, &b_TopBJet2_eta);
   fChain->SetBranchAddress("TopBJet2_mass", &TopBJet2_mass, &b_TopBJet2_mass);
   fChain->SetBranchAddress("TopBJet2_phi", &TopBJet2_phi, &b_TopBJet2_phi);
   fChain->SetBranchAddress("TopBJet2_pt", &TopBJet2_pt, &b_TopBJet2_pt);
   fChain->SetBranchAddress("TopBestJet_eta", &TopBestJet_eta, &b_TopBestJet_eta);
   fChain->SetBranchAddress("TopBestJet_mass", &TopBestJet_mass, &b_TopBestJet_mass);
   fChain->SetBranchAddress("TopBestJet_phi", &TopBestJet_phi, &b_TopBestJet_phi);
   fChain->SetBranchAddress("TopBestJet_pt", &TopBestJet_pt, &b_TopBestJet_pt);
   fChain->SetBranchAddress("TopJet1_eta", &TopJet1_eta, &b_TopJet1_eta);
   fChain->SetBranchAddress("TopJet1_mass", &TopJet1_mass, &b_TopJet1_mass);
   fChain->SetBranchAddress("TopJet1_phi", &TopJet1_phi, &b_TopJet1_phi);
   fChain->SetBranchAddress("TopJet1_pt", &TopJet1_pt, &b_TopJet1_pt);
   fChain->SetBranchAddress("TopJet2_eta", &TopJet2_eta, &b_TopJet2_eta);
   fChain->SetBranchAddress("TopJet2_mass", &TopJet2_mass, &b_TopJet2_mass);
   fChain->SetBranchAddress("TopJet2_phi", &TopJet2_phi, &b_TopJet2_phi);
   fChain->SetBranchAddress("TopJet2_pt", &TopJet2_pt, &b_TopJet2_pt);
   fChain->SetBranchAddress("TopJet3_eta", &TopJet3_eta, &b_TopJet3_eta);
   fChain->SetBranchAddress("TopJet3_mass", &TopJet3_mass, &b_TopJet3_mass);
   fChain->SetBranchAddress("TopJet3_phi", &TopJet3_phi, &b_TopJet3_phi);
   fChain->SetBranchAddress("TopJet3_pt", &TopJet3_pt, &b_TopJet3_pt);
   fChain->SetBranchAddress("W_TransverseMass", &W_TransverseMass, &b_W_TransverseMass);
   fChain->SetBranchAddress("W_charge", &W_charge, &b_W_charge);
   fChain->SetBranchAddress("W_eta", &W_eta, &b_W_eta);
   fChain->SetBranchAddress("W_mass", &W_mass, &b_W_mass);
   fChain->SetBranchAddress("W_phi", &W_phi, &b_W_phi);
   fChain->SetBranchAddress("W_pt", &W_pt, &b_W_pt);
   fChain->SetBranchAddress("isLeptonTriggerPassed", &isLeptonTriggerPassed, &b_isLeptonTriggerPassed);
   fChain->SetBranchAddress("isLeptonTriggerMatched", &isLeptonTriggerMatched, &b_isLeptonTriggerMatched);
   fChain->SetBranchAddress("isNoEleMuOverlap", &isNoEleMuOverlap, &b_isNoEleMuOverlap);
   Notify();
}

Bool_t ElectronPreTagTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ElectronPreTagTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ElectronPreTagTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ElectronPreTagTree_cxx
