/******************************************************************************
 * HistogrammingWtDiLepTopologyCutPtEta.hpp                                           *
 *                                                                            *
 * Fills a tree with specific variables.                                      *
 * Should be called after preselection cuts in bin/analysis                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingWtDiLepTopologyCutPtEta()    -- Parameterized Constructor          *
 *    ~HistogrammingWtDiLepTopologyCutPtEta()   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingWtDiLepTopologyCutPtEta"    *
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
 *      Jan 2011, First Version for Wt MVA analysis   Huaqiao ZHANG @ CERN    *
 *****************************************************************************/

#ifndef HistogrammingWtDiLepTopologyCutPtEta_h
#define HistogrammingWtDiLepTopologyCutPtEta_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;


class HistogrammingWtDiLepTopologyCutPtEta : public HistoCut 
{

public:

  HistogrammingWtDiLepTopologyCutPtEta(EventContainer *obj,TFile *skimFilez);
  
  // destructor
  ~HistogrammingWtDiLepTopologyCutPtEta();

  // methods:
  void BookHistogram();
  bool Apply();
  
  inline std::string GetCutName (void) { return "HistogrammingWtDiLepTopologyCutPtEta"; };

private:
  TFile *_skimFile;
  TTree *_skimTree;
  std::vector<Double_t> _bweightvec;
  double xy_to_phi(double x, double y);
  double dPhi(double phi1,double phi2);

  //------------------------------------------------------
  // Event Weight
  //------------------------------------------------------
  Float_t _EventWeight;
  Int_t _RunNumber;
  Int_t _EventNumber;
  Int_t _EventType;
  Int_t _HFOR;
 
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
  Float_t _TransverseMass_Jet1Jet2;
  Float_t _Shat;
    
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

  Int_t _JetTruthLabel;
  //------------------------------------------------------
  // Missing Et 
  //------------------------------------------------------
  Float_t _MET;
  Float_t _METx;
  Float_t _METy;
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
  Int_t   _LeadingLeptonIsTight_CutPtEta;
  Float_t _LeadingLeptonEtCone20_CutPtEta;
  Float_t _LeadingLeptonEtCone30_CutPtEta;
  Float_t _LeadingLeptonEtCone40_CutPtEta;
  Float_t _LeadingLeptonPtCone30_CutPtEta;
  Float_t _LeadingLeptonEtRatioCone30_CutPtEta;
  Float_t _LeadingLeptonPtRatioCone30_CutPtEta;
  Float_t _LeadingLeptonJetDeltaRMin_CutPtEta;
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
  Int_t   _SubLeadingLeptonIsTight_CutPtEta;
  Float_t _SubLeadingLeptonEtCone20_CutPtEta;
  Float_t _SubLeadingLeptonEtCone30_CutPtEta;
  Float_t _SubLeadingLeptonEtCone40_CutPtEta;
  Float_t _SubLeadingLeptonPtCone30_CutPtEta;
  Float_t _SubLeadingLeptonEtRatioCone30_CutPtEta;
  Float_t _SubLeadingLeptonPtRatioCone30_CutPtEta;
  Float_t _SubLeadingLeptonJetDeltaRMin_CutPtEta;
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

  //------------------------------------------------------
  // Global variable
  //------------------------------------------------------
  Float_t _ThrustEta;
  Float_t _ThrustPhi;
  Float_t _Thrust;

  Float_t _A_Thrust_Hemi;
  Float_t _A_LeadingJet_Hemi;

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


