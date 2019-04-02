/******************************************************************************
 * HistogrammingWtDiLepTopology.hpp                                           *
 *                                                                            *
 * Fills a tree with specific variables.                                      *
 * Should be called after preselection cuts in bin/analysis                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingWtDiLepTopology()    -- Parameterized Constructor          *
 *    ~HistogrammingWtDiLepTopology()   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingWtDiLepTopology"*
 *                                                                            *
 * Private Data Members of this class                                         *
 *        TFile *_skimFile                                                    *
 *        TTree *_skimTree;                                                   *
 *        Lots of variables, histograms                                       *
 *                                                                            *
 * History                                                                    *
 *  arguments: - pointer to input particles class                             *
 *             - write: whether to write out the NN skim or not               *
 *             - name: name of the NNvars skim file                           *
 *             - doSpin: whether to include top spin correlation vars in NN   *
 *                                                                            *
 * History                                                                    *
 *      Jan 2011, First Version for Wt MVA analysis   Huaqiao ZHANG @ CERN    *
 *****************************************************************************/

#ifndef HistogrammingWtDiLepTopology_h
#define HistogrammingWtDiLepTopology_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
//#include "SingleTopRootAnalysis/Base/Dictionary/ElecSF.hpp" 
//#include "SingleTopRootAnalysis/Base/Dictionary/MuonSFv09.hpp" 
#include "TLorentzVector.h"
#include "TMath.h"

#include <string>
#include <vector>
#include <sstream>

#include "TLorentzVector.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TVectorD.h"
//#include "TMinuit.h"


using namespace std;


class HistogrammingWtDiLepTopology : public HistoCut
{

public:

  HistogrammingWtDiLepTopology(EventContainer *obj,TFile *skimFilez, int whichstream=-1);

  // destructor
  ~HistogrammingWtDiLepTopology();

  // methods:
  void BookHistogram();
  bool Apply();

  inline std::string GetCutName (void) { return "HistogrammingWtDiLepTopology"; };

private:
  TFile *_skimFile;
  TTree *_skimTree;
  std::vector<Double_t> _bweightvec;
  double xy_to_phi(double x, double y);
  double dPhi(double phi1,double phi2);
  void soulveZtautau(TLorentzVector lep1, TLorentzVector lep2, double METX, double METY, TLorentzVector & nv1, TLorentzVector & nv2);
  bool getWtSys(vector<TLorentzVector> & W_neu, vector<TLorentzVector> & T_W_neu, TLorentzVector T_b, TLorentzVector T_W_lep, TLorentzVector W_lep,double MEx, double MEy);
  //------------------------------------------------------
  // Event Weight
  //------------------------------------------------------
  Float_t _EventWeight;
  Float_t _PileupWeight; //pileUp weight per event
  Int_t _RunNumber;
  Int_t _EventNumber;
  Float_t mulbn;
  Int_t Pvtx_n;
  Float_t actualIntPerXing;
  Float_t averageIntPerXing;
  Int_t _EventType;
  Int_t _HFOR;
  Bool_t _isLArbadjet;
  Bool_t _isLArbadelectron;
  Bool_t _isLArbad;
  //MuonSFv09 MUSFobj;
  //ElecSF ElecSFobj;

  Int_t _theStream;
  Int_t _stream;
  Int_t _fireElecTrig;
  Int_t _fireMuonTrig;

  Int_t _NMuonTight;
  Int_t _NMuonVeto;
  Int_t _NMuonAll;

  Int_t _NElectronTight;
  Int_t _NElectronVeto;
  Int_t _NElectronAll;
  Int_t pdgid;

  Int_t _NWDecayLeptons;
  Int_t _NWDecayLightQuarks;

  Float_t _TaggingWeight; //tagging weight per event
  Float_t _TaggingWeight_Bup;//tagging weight per event applying the b SF up shift
  Float_t _TaggingWeight_Bdown;//tagging weight per event applying the b SF down shift
  Float_t _TaggingWeight_Lqup;//tagging weight per event applying the Lq SF up shift
  Float_t _TaggingWeight_Lqdown;//tagging weight per event applying the Lq SF down shift
  //------------------------------------------------------
  // System Mass
  //------------------------------------------------------
  Float_t _InvariantMass_AllJets;
  Float_t _InvariantMass_Jet1Jet2Lepton;
  Float_t _InvariantMass_AllJetsLepton;
  Float_t _InvariantMass_AllJets_MinusBTaggedJet;
  Float_t _InvariantMass_Jet1Jet2;
  Float_t _InvariantMass_Lep1Lep2;
  Float_t _InvariantMass_METcoLinearLep1Lep2;
  Float_t _TransverseMass_Jet1Jet2;
  Float_t _InvariantMass_Lep1Jet1;
  Float_t _InvariantMass_Lep2Jet1;
  Float_t _InvariantMass_Lep1Jet2;
  Float_t _InvariantMass_Lep2Jet2;
  Float_t _InvariantMass_MaxLepJet1;
  Float_t _InvariantMass_MinLepJet1;
  Float_t _Shat;

  //------------------------------------------------------
  // Reco System 
  //------------------------------------------------------
  Int_t  _reco_n;
  Float_t  _reco_mWt     [8];
  Float_t  _reco_Wt_pT   [8];
  Float_t  _reco_Wt_pz   [8];
  Float_t  _reco_WnzTnz  [8];
  Float_t  _reco_top_pT  [8];
  Float_t  _reco_W_pT    [8];
  Float_t  _reco_topW_pT [8];
  Float_t  _reco_top_pz  [8];
  Float_t  _reco_W_pz    [8];
  Float_t  _reco_topW_pz [8];


  Float_t  _M_T1;
  Float_t  _M_T2;
  Float_t  _M_T3;
  Float_t  _M_T4;
  Float_t  _M_T5;

  //------------------------------------------------------
  // System 

  //------------------------------------------------------
  Float_t _pT_sys;
  Float_t _pT_sys_sig;
  Float_t _pT_sys_lepsJet1;
  Float_t _eta_sys_lepsJet1;
  //------------------------------------------------------
  // H
  //------------------------------------------------------
  Float_t _H_Jet1Jet2Lepton;
  Float_t _H_Jet1Jet2Lepton_MinusBJet;
  Float_t _H_AllJetsLepton;
  Float_t _H_AllJets;
  Float_t _H_AllJets_MinusBestJet;
  Float_t _H_AllJets_MinusBTaggedJet;
  Float_t _H_Jet1Jet2;

  //------------------------------------------------------
  // Ht
  //------------------------------------------------------
  Float_t _HT;
  Float_t _HT_Jet1Jet2LeptonMET;
  Float_t _HT_LeptonMET;
  Float_t _HT_Jet1Jet2Lepton_MinusBJet;
  Float_t _HT_LastJet;
  Float_t _P_LastJet;
  Float_t _HT_AllJets_MinusBTaggedJet;
  Float_t _HT_Jet1Jet2;
  Float_t _HT_AllJetsLeptonMET;
  Float_t _HT_Jet1LeptonMET;
  Float_t _HT_AllJets;
  Float_t _HT_RHT;  //_HT_Jet1LeptonMET/_HT_AllJetsLeptonMET
  Float_t _HT_AllLeps;

  //------------------------------------------------------
  // Jet
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
  Float_t _Pt_Jet1Jet2;
  Float_t _DeltaPt_Jet1Jet2;
  Float_t _Pt_AllJets_MinusBTaggedJet;
  int _NJets;
  int _NTaggedJets;
  int _NUntaggedJets;

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

  Float_t _JetFitterCombNNWeight1;
  Float_t _JetFitterCombNNWeight2;
  Float_t _JetFitterCombNNWeight3;
  Float_t _JetFitterCombNNWeight4;

  Int_t _JetTruthLabel;
  //------------------------------------------------------
  // Missing Et
  //------------------------------------------------------
  Float_t _MET;
  Float_t _METx;
  Float_t _METy;
  Float_t _METPhi;

  //------------------------------------------------------
  // Lepton
  //------------------------------------------------------
  Float_t _LeadingLeptonID;
  Float_t _LeadingLeptonPt;
  Float_t _LeadingLeptonCharge;
  Float_t _LeadingLeptonEta;
  Float_t _LeadingLeptonPhi;
  Float_t _LeadingLeptonEtCone20;
  Float_t _LeadingLeptonEtCone30;
  Float_t _LeadingLeptonEtCone40;
  Float_t _LeadingLeptonPtCone30;
  Float_t _LeadingLeptonDetEta;
  Float_t _LeadingLeptonDetEt;
  Float_t _LeadingLeptonDetE;
  Float_t _LeadingLeptonDetPhi;
  Int_t   _LeadingLeptonIsTight;
  Int_t   _LeadingIsRobustMedium;
  Int_t   _LeadingIsemTrackBlayer;
  Int_t   _LeadingIsemTrackMatchEoverP;
  Int_t   _LeadingLeptonTriggerIn;

  Float_t _SubLeadingLeptonID;
  Float_t _SubLeadingLeptonPt;
  Float_t _SubLeadingLeptonCharge;
  Float_t _SubLeadingLeptonEta;
  Float_t _SubLeadingLeptonPhi;
  Float_t _SubLeadingLeptonEtCone20;
  Float_t _SubLeadingLeptonEtCone30;
  Float_t _SubLeadingLeptonEtCone40;
  Float_t _SubLeadingLeptonPtCone30;
  Float_t _SubLeadingLeptonDetEta;
  Float_t _SubLeadingLeptonDetEt;
  Float_t _SubLeadingLeptonDetE;
  Float_t _SubLeadingLeptonDetPhi;
  Int_t   _SubLeadingLeptonIsTight;
  Int_t   _SubLeadingIsRobustMedium;
  Int_t   _SubLeadingIsemTrackBlayer;
  Int_t   _SubLeadingIsemTrackMatchEoverP;
  Int_t   _SubLeadingLeptonTriggerIn;

  //------------------------------------------------------
  // UnTagged Jet Histograms
  //------------------------------------------------------
  Float_t _LeadingUntaggedJetPt;
  Float_t _SecondUntaggedJetPt;
  Float_t _LeadingUntaggedJetEta;
  Float_t _DeltaPt_tag_untag;
  Float_t _DeltaRLeadingUntaggedJetLepton;
  //------------------------------------------------------
  // b-Tagged Jet
  //------------------------------------------------------
  Float_t _LeadingBTaggedJetPt;
  Float_t _SecondBTaggedJetPt;
  Float_t _LeadingBTaggedJetEta;
  Float_t _SecondBTaggedJetEta;
  Float_t _DeltaRLeadingBTaggedJetLepton;
  Float_t _DeltaRSecondBTaggedJetLepton;
  Float_t _DeltaRBTaggedJet1BTaggedJet2;
  //------------------------------------------------------
  // Angular Variable
  //------------------------------------------------------
  Float_t _Cos_UntaggedJetLepton_BTaggedTop;
  Float_t _Cos_Jet1AllJets_AllJets;
  Float_t _Cos_BTaggedJetAllJets_AllJets;
  Float_t _Cos_UntaggedJetAllJets_AllJets;
  Float_t _Cos_NotBestJetAllJets_AllJets;
  Float_t _Cos_LeptonQZ_BestTop;
  Float_t _Cos_wHelicityThetaStar;
  Float_t _Cos_Jet1ThetaStar;
  Float_t _DeltaRJet1Jet2;
  Float_t _DeltaPhiLep1Lep2;
  Float_t _DeltaPhiLeadingJetMissingEt;
  Float_t _DeltaPhiLeadingLeptonMissingEt;
  Float_t _DeltaPhiSubLeadingLeptonMissingEt;
  Float_t _DeltaPhiLepton1Lepton2;
  Float_t _DeltaPhiLeadingLeptonLeadingJet;
  Float_t _DeltaPhiSubLeadingLeptonLeadingJet;
  Float_t _DeltaMaxPhiLeptonLeadingJet;
  Float_t _DeltaMinPhiLeptonLeadingJet;
  Float_t _DeltaEtaLeadingLeptonLeadingJet;
  Float_t _DeltaEtaSubLeadingLeptonLeadingJet;
  Float_t _DeltaMaxEtaLeptonLeadingJet;
  Float_t _DeltaMinEtaLeptonLeadingJet;
  Float_t _DeltaRLeadingLeptonLeadingJet;
  Float_t _DeltaRSubLeadingLeptonLeadingJet;
  Float_t _DeltaMaxRLeptonLeadingJet;
  Float_t _DeltaMinRLeptonLeadingJet;
  Float_t _DeltaEtaLep1Lep2;
  Float_t _DeltaEta_SLep1Jet1_Lep2;
  Float_t _DeltaEta_SLep1Lep2_Jet1;
  Float_t _DeltaEta_SLep2Jet1_Lep1;
  Float_t _DeltaPhi_SLep1Jet1_Lep2;
  Float_t _DeltaPhi_SLep1Lep2_Jet1;
  Float_t _DeltaPhi_SLep2Jet1_Lep1;
  Float_t _DeltaRLep1Lep2;
  Float_t _DeltaR_SLep1Jet1_Lep2;
  Float_t _DeltaR_SLep1Lep2_Jet1;
  Float_t _DeltaR_SLep2Jet1_Lep1;

  //------------------------------------------------------
  // Reconstructed variables
  //------------------------------------------------------

  Float_t _SystemLep1Jet1_E;
  Float_t _SystemLep1Jet1_eta;
  Float_t _SystemLep1Jet1_phi;
  Float_t _SystemLep1Jet1_pt;
  Float_t _SystemLep1Lep2_E;
  Float_t _SystemLep1Lep2_eta;
  Float_t _SystemLep1Lep2_phi;
  Float_t _SystemLep1Lep2_pt;
  Float_t _SystemLep2Jet1_E;
  Float_t _SystemLep2Jet1_eta;
  Float_t _SystemLep2Jet1_phi;
  Float_t _SystemLep2Jet1_pt;

  //------------------------------------------------------
  // Global variable
  //------------------------------------------------------
  Float_t _ThrustEta;
  Float_t _ThrustPhi;
  Float_t _Thrust;

  Float_t _A_Thrust_Hemi;
  Float_t _A_LeadingJet_Hemi;

  //------------------------------------------------------
  Float_t Aplanarity_AllJetsLepton;
  Float_t Sphericity_AllJetsLepton;
  //------------------------------------------------------


  //------------------------------------------------------
  // Dilepton SF variables
  //------------------------------------------------------

  Float_t _DilepSFNomNom;
  Float_t _DilepSFUpNom;
  Float_t _DilepSFDownNom;
  Float_t _DilepSFNomUp;
  Float_t _DilepSFNomDown;
  Float_t _DilepSFUpUp;
  Float_t _DilepSFDownDown;


  ////////////////////////////////////////////////////////
  // Histograms declaration
  ////////////////////////////////////////////////////////

  myTH1F* _hP_LastJet;
  myTH1F* _hWPt;
  myTH1F* _hEventWeight;
  myTH1F* _hRunNumber;
  myTH1F* _hEventNumber;
  myTH1F* _hEventType;
  myTH1F* _hHT_LastJet;
  myTH1F* _hAbsLogEventWeight;
  myTH1F* _hHFOR;

  myTH1F* _hNWDecayLeptons;
  myTH1F* _hNWDecayLightQuarks;

  myTH1F* _hTaggingWeight;
  myTH1F* _hTaggingWeight_Bup;
  myTH1F* _hTaggingWeight_Bdown;
  myTH1F* _hTaggingWeight_Lqup;
  myTH1F* _hTaggingWeight_Lqdown;

  //------------------------------------------------------
  // System Mass Histograms
  //------------------------------------------------------
  myTH1F* _hInvariantMass_AllJets;
  myTH1F* _hInvariantMass_Jet1Jet2Lepton;
  myTH1F* _hInvariantMass_AllJetsLepton;
  myTH1F* _hInvariantMass_AllJets_MinusBTaggedJet;
  myTH1F* _hInvariantMass_Jet1Jet2;
  myTH1F* _hTransverseMass_Jet1Jet2;
  myTH1F* _hShat;

  //------------------------------------------------------
  // H Histograms
  //------------------------------------------------------
  myTH1F* _hH_Jet1Jet2Lepton;
  myTH1F* _hH_Jet1Jet2Lepton_MinusBJet;
  myTH1F* _hH_AllJetsLepton;
  myTH1F* _hH_AllJets;
  myTH1F* _hH_AllJets_MinusBestJet;
  myTH1F* _hH_AllJets_MinusBTaggedJet;
  myTH1F* _hH_Jet1Jet2;
  //------------------------------------------------------
  // Ht Histograms
  //------------------------------------------------------
  myTH1F* _hHT;  // the sum of the ETs of the lepton, MET, and all jets.
  myTH1F* _hHT_Jet1Jet2LeptonMET;
  myTH1F* _hHT_Jet1LeptonMET;
  myTH1F* _hHT_LeptonMET;
  myTH1F* _hHT_Jet1Jet2Lepton_MinusBJet;
  myTH1F* _hHT_AllJetsLeptonMET;
  myTH1F* _hHT_AllJets;
  myTH1F* _hHT_AllJets_MinusBTaggedJet;
  myTH1F* _hHT_Jet1Jet2;
  myTH1F* _hHT_RHT;

  //------------------------------------------------------
  // Jet Histograms
  //------------------------------------------------------
  myTH1F* _hJet1Pt;
  myTH1F* _hJet2Pt;
  myTH1F* _hJet3Pt;
  myTH1F* _hJet4Pt;
  myTH1F* _hJet1Eta;
  myTH1F* _hJet2Eta;
  myTH1F* _hJet3Eta;
  myTH1F* _hJet4Eta;
  myTH1F* _hJet1Phi;
  myTH1F* _hJet2Phi;
  myTH1F* _hJet3Phi;
  myTH1F* _hJet4Phi;
  myTH1F* _hPt_Jet1Jet2;
  myTH1F* _hDeltaPt_Jet1Jet2;
  myTH1F* _hNJets;
  myTH1F* _hNGoodJets;
  myTH1F* _hNTaggedJets;
  myTH1F* _hNUntaggedJets;
  myTH1F* _hJet12_Centrality;
  myTH1F* _hAllJetsLepton_Centrality;
  myTH1F* _hLastJetLepton_Centrality;
  myTH1F* _hPt_AllJets_MinusBTaggedJet;
  myTH1F* _hMaxJetEta;
  myTH1F* _hMinJetEta;

  myTH1F* _hJetProbWeight1;
  myTH1F* _hJetProbWeight2;
  myTH1F* _hJetProbWeight3;
  myTH1F* _hJetProbWeight4;

  myTH1F* _hSV0Weight1;
  myTH1F* _hSV0Weight2;
  myTH1F* _hSV0Weight3;
  myTH1F* _hSV0Weight4;

  myTH1F* _hJetFitterCombNNWeight1;
  myTH1F* _hJetFitterCombNNWeight2;
  myTH1F* _hJetFitterCombNNWeight3;
  myTH1F* _hJetFitterCombNNWeight4;

  myTH1F* _hJetTruthLabel;
  //------------------------------------------------------
  // Missing Et Histograms
  //------------------------------------------------------
  myTH1F* _hMET;
  //------------------------------------------------------
  // Lepton Histograms
  //------------------------------------------------------
  myTH1F* _hLeadingLeptonPt;

  myTH1F* _hLeadingLeptonEta;
  myTH1F* _hLeadingLeptonPhi;
  myTH1F* _hLeadingLeptonIsTight;
  myTH1F* _hLeadingLeptonEtCone20;
  myTH1F* _hLeadingLeptonEtCone30;
  myTH1F* _hLeadingLeptonEtCone40;
  myTH1F* _hLeadingLeptonPtCone30;
  myTH1F* _hSubLeadingLeptonPt;

  myTH1F* _hSubLeadingLeptonEta;
  myTH1F* _hSubLeadingLeptonPhi;
  myTH1F* _hSubLeadingLeptonIsTight;
  myTH1F* _hSubLeadingLeptonEtCone20;
  myTH1F* _hSubLeadingLeptonEtCone30;
  myTH1F* _hSubLeadingLeptonEtCone40;
  myTH1F* _hSubLeadingLeptonPtCone30;
  myTH1F* _hIsemTrackMatchEoverP;
  //------------------------------------------------------
  // UnTagged Jet Histograms
  //------------------------------------------------------
  myTH1F* _hLeadingUntaggedJetPt;
  myTH1F* _hSecondUntaggedJetPt;
  myTH1F* _hLeadingUntaggedJetEta;
  myTH1F* _hDeltaPt_tag_untag;
  myTH1F* _hDeltaRLeadingUntaggedJetLepton;
  //------------------------------------------------------
  // b-Tagged Jet Histograms
  //------------------------------------------------------
  myTH1F* _hLeadingBTaggedJetPt;
  myTH1F* _hSecondBTaggedJetPt;
  myTH1F* _hLeadingBTaggedJetEta;
  myTH1F* _hSecondBTaggedJetEta;
  myTH1F* _hDeltaRLeadingBTaggedJetLepton;
  myTH1F* _hDeltaRSecondBTaggedJetLepton;
  myTH1F* _hDeltaRBTaggedJet1BTaggedJet2;
  //------------------------------------------------------
  // Angular Variable Histograms
  //------------------------------------------------------
  myTH1F* _hCos_UntaggedJetLepton_BTaggedTop;
  myTH1F* _hCos_Jet1AllJets_AllJets;
  myTH1F* _hCos_BTaggedJetAllJets_AllJets;
  myTH1F* _hCos_UntaggedJetAllJets_AllJets;
  myTH1F* _hCos_NotBestJetAllJets_AllJets;
  myTH1F* _hCos_LeptonQZ_BestTop;
  myTH1F* _hCos_wHelicityThetaStar;
  myTH1F* _hCos_Jet1ThetaStar;
  myTH1F* _hDeltaRJet1Jet2;
  myTH1F* _hDeltaPhiLeadingJetMissingEt;
  myTH1F* _hDeltaPhiLeadingLeptonMissingEt;
  myTH1F* _hDeltaPhiSubLeadingLeptonMissingEt;
  myTH1F* _hDeltaPhiLepton1Lepton2;
  myTH1F* _hDeltaPhiLeadingLeptonLeadingJet;
  myTH1F* _hDeltaPhiSubLeadingLeptonLeadingJet;
  myTH1F* _hDeltaMaxPhiLeptonLeadingJet;
  myTH1F* _hDeltaMinPhiLeptonLeadingJet;
  myTH1F* _hDeltaEtaLeadingLeptonLeadingJet;
  myTH1F* _hDeltaEtaSubLeadingLeptonLeadingJet;
  myTH1F* _hDeltaMaxEtaLeptonLeadingJet;
  myTH1F* _hDeltaMinEtaLeptonLeadingJet;
  myTH1F* _hDeltaRLeadingLeptonLeadingJet;
  myTH1F* _hDeltaRSubLeadingLeptonLeadingJet;
  myTH1F* _hDeltaMaxRLeptonLeadingJet;
  myTH1F* _hDeltaMinRLeptonLeadingJet;
  //------------------------------------------------------
  // Global variable
  //------------------------------------------------------
  myTH1F* _hThrustEta;
  myTH1F* _hThrustPhi;
  myTH1F* _hThrust;

  myTH1F* _hA_Thrust_Hemi;
  myTH1F* _hA_LeadingJet_Hemi;

};


#endif


