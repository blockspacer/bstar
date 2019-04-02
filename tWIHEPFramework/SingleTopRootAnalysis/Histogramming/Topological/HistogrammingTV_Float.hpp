/******************************************************************************
 * HistogrammingTV_Float.hpp                                                  *
 *                                                                            *
 * Fills a tree with specific variables.                                      *
 * Should be called after preselection cuts in bin/analysis                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTV_Float()              -- Parameterized Constructor       *
 *    ~HistogrammingTV_Float()             -- Destructor                      *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTV_Float"    *
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

#ifndef HistogrammingTV_Float_h
#define HistogrammingTV_Float_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/ScaleWjets.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;


class HistogrammingTV_Float : public HistoCut 
{

public:

  HistogrammingTV_Float(EventContainer *obj,TFile *skimFilez);
  
  // destructor
  ~HistogrammingTV_Float();

  // methods:
  void BookHistogram();
  bool Apply();
  
  inline std::string GetCutName (void) { return "HistogrammingTV_Float"; };

private:
  TFile *_skimFile;
  TTree *_skimTree;
  std::vector<Double_t> _bweightvec;

  // skim tree variables
  Float_t _MuonSF;
  Float_t _MuonSFError;

 ScaleWjets* scpt10;
  ScaleWjets* scpt20;
  ScaleWjets* scopt2;
  ScaleWjets* scopt3;
  TLorentzVector *jet;

 Float_t jetweight_opt2;
  Float_t jetweight_opt3;
  Float_t jetweight_pt10;
  Float_t jetweight_pt20;


  Float_t actualIntPerXing;
  Float_t averageIntPerXing;

  //------------------------------------------------------
  // Event Weight
  //------------------------------------------------------
  Float_t _EventWeight;
  Int_t _RunNumber;
  Int_t _EventNumber;
  Int_t _HFOR;
 
  Int_t _NMuonTight;
  Int_t _NMuonVeto;
  Int_t _NMuonAll;
 Int_t pdgid;
  Int_t _Pvtx_n;

  Int_t _NWDecayLeptons;
  Int_t _NWDecayLightQuarks;

 Float_t _PileupWeight; //tagging weight per event
  Float_t _TaggingWeight; //tagging weight per event
  Float_t _TaggingWeight_Bup;//tagging weight per event applying the b SF up shift
  Float_t _TaggingWeight_Bdown;//tagging weight per event applying the b SF down shift
  Float_t _TaggingWeight_Lqup;//tagging weight per event applying the Lq SF up shift
  Float_t _TaggingWeight_Lqdown;//tagging weight per event applying the Lq SF down shift
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

  Float_t _DeltaEtaJet1Jet2;
  Float_t _DeltaEtaBJet1Lepton;

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
 Float_t _LeptonPhi;
 Float_t _ElectronDetEta;
 Float_t _ElectronDetEt;
 Float_t _ElectronDetPhi;
 Float_t _ElectronDetE;
  //  Float_t _QTimesEtaLepton;
  Float_t _LeptonChargePDGID;
  //Float_t _EtCone20;
  //------- 
  Int_t _LeptonIsTight;
  Float_t _LeptonEtCone20;
  Float_t _LeptonEtCone30;
  Float_t _LeptonEtCone40;
  Float_t _LeptonPtCone30;
  Float_t _LeptonEtRatioCone30;
  Float_t _LeptonPtRatioCone30;
  Int_t  _ElectronQuality;
  Int_t  _IsRobustMedium;
  Int_t _IsemTrackBlayer;
  Int_t _IsemTrackMatchEoverP;
  Int_t _LeptonIsTight_CutPtEta;
  Float_t _LeptonEtCone20_CutPtEta;
  Float_t _LeptonEtCone30_CutPtEta;
  Float_t _LeptonEtCone40_CutPtEta;
  Float_t _LeptonPtCone30_CutPtEta;
  Float_t _LeptonEtRatioCone30_CutPtEta;
  Float_t _LeptonPtRatioCone30_CutPtEta;
  Float_t _LeptonJetDeltaRMin_CutPtEta;
  Int_t  _ElectronQuality_CutPtEta;
  Int_t  _IsRobustMedium_CutPtEta;
  Int_t _IsemTrackBlayer_CutPtEta;
  Int_t _IsemTrackMatchEoverP_CutPtEta;
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


  ////////////////////////////////////////////////////////
  // Histograms declaration 
  ////////////////////////////////////////////////////////

  myTH1F* _hP_LastJet;
  myTH1F* _hWPt;
  myTH1F* _hEventWeight;
  myTH1F* _hRunNumber;
  myTH1F* _hEventNumber;
  myTH1F* _hHT_LastJet;
  myTH1F* _hAbsLogEventWeight;
  myTH1F* _hHFOR;
  myTH1F* _hPvtx_n;
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
  myTH1F* _hInvariantMass_Jet1Jet2LeptonMET;
  myTH1F* _hInvariantMass_AllJetsLeptonMET;
  myTH1F* _hInvariantMass_LeptonMETBestJet;
  myTH1F* _hInvariantMass_AllJets_MinusBestJet;
  myTH1F* _hInvariantMass_AllJets_MinusBTaggedJet;
  myTH1F* _hInvariantMass_Jet1Jet2;
  myTH1F* _hTransverseMass_Jet1Jet2;
  myTH1F* _hShat;
  // myTH1F* _hNeutrino_Mass;
  myTH1F* _hNeutrino_Pz;
  myTH1F* _hHT;  // the sum of the ETs of the lepton, MET, and all jets.

  //------------------------------------------------------
  // H Histograms
  //------------------------------------------------------
  myTH1F* _hH_Jet1Jet2LeptonMET;
  myTH1F* _hH_Jet1Jet2LeptonMET_MinusBJet;
  myTH1F* _hH_AllJetsLeptonMET;
  myTH1F* _hH_AllJets;
  myTH1F* _hH_AllJets_MinusBestJet;
  myTH1F* _hH_AllJets_MinusBTaggedJet;
  myTH1F* _hH_Jet1Jet2;
    
  //------------------------------------------------------
  // Ht Histograms
  //------------------------------------------------------
  myTH1F* _hHT_Jet1Jet2LeptonMET;
  myTH1F* _hHT_LeptonMET;
  myTH1F* _hHT_Jet1Jet2LeptonMET_MinusBJet;
  myTH1F* _hHT_AllJetsLeptonMET;
  myTH1F* _hHT_AllJetsLepton;
  myTH1F* _hHT_AllJets;
  myTH1F* _hHT_AllJets_MinusBestJet;
  myTH1F* _hHT_AllJets_MinusBTaggedJet;
  myTH1F* _hHT_Jet1Jet2;
 
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
  myTH1F* _hJet1Pt_NotBest;
  myTH1F* _hJet2Pt_NotBest;
  myTH1F* _hPt_Jet1Jet2;
  myTH1F* _hPt_AllJets_MinusBestJet;
  myTH1F* _hPt_AllJets_MinusBTaggedJet;
  myTH1F* _hDeltaPt_Jet1Jet2;
  myTH1F* _hNJets;
  myTH1F* _hNGoodJets;
  myTH1F* _hNTaggedJets;
  myTH1F* _hNUntaggedJets;
  myTH1F* _hBestJetPt;
  myTH1F* _hDeltaRJet1Jet2;
  myTH1F* _hDeltaPhiLeadingJetLepton;
  myTH1F* _hDeltaPhiLeadingJetMissingEt;
  myTH1F* _hDeltaRLeadingJetLepton;
  myTH1F* _hDeltaPhiLeptonMissingEt;
  myTH1F* _hJet12_Centrality;
  myTH1F* _hAllJetsLepton_Centrality;
  myTH1F* _hLastJetLepton_Centrality;
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
  //myTH1F* _hMETPt;
  myTH1F* _hMET;
  //------------------------------------------------------
  // W Histograms
  //------------------------------------------------------
  myTH1F* _hWTransverseMass;
   //myTH1F* _hWPt;
  myTH1F* _hWCharge;
     
  //------------------------------------------------------
  // Top Histograms
  //------------------------------------------------------
  myTH1F* _hBestTopMass;
  myTH1F* _hBTaggedTopMass;
  myTH1F* _hLeadingJetTopMass;
    
  //------------------------------------------------------
  // Lepton Histograms
  //------------------------------------------------------
  myTH1F* _hLeptonPt;
  myTH1F* _hLeptonChargePDGID;
  myTH1F* _hLeptonEta;
  myTH1F* _hLeptonPhi;
  myTH1F* _hElectronDetEta;
  myTH1F* _hElectronDetEt;
  myTH1F* _hElectronDetPhi;
  myTH1F* _hElectronDetE;
  // myTH1F* _hQTimesEtaLepton;
  //myTH1F* _hEtCone20;
  myTH1F* _hLeptonIsTight;
  myTH1F* _hLeptonEtCone20;
   myTH1F* _hLeptonEtCone30;
   myTH1F* _hLeptonEtCone40;
   myTH1F* _hLeptonPtCone30;
   myTH1F*  _hElectronQuality;
   myTH1F*  _hIsRobustMedium;
   myTH1F* _hIsemTrackBlayer;
   myTH1F* _hIsemTrackMatchEoverP;
  //------------------------------------------------------
  // UnTagged Jet Histograms
  //------------------------------------------------------
  myTH1F* _hLeadingUntaggedJetPt;
  myTH1F* _hSecondUntaggedJetPt;
  myTH1F* _hLeadingUntaggedJetEta;
  //myTH1F* _hQTimesEta;
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


};


#endif


