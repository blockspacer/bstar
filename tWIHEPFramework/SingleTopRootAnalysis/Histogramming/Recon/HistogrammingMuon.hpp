/******************************************************************************
 * HistogrammingMuon.hpp                                                      *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Must pass a paramter describing the type of muon you want to plot          *
 *   Choices are All, Veto, Tight, and UnIsolated                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMuon()              -- Parameterized Constructor           *
 *    ~HistogrammingMuon()             -- Destructor                          *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMuon"        *
 *                                                                            *
 * Private Data Members of HistogrammingMuon class                            *
 * - lots of histograms of muon propreties                                    *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      17 Jul 2007 - J. Holzbauer added MuonInJet histogram                  *
 *      18 Jan 2007 - P. Ryan - modified for all/tight/veto/unisolated        *
 *****************************************************************************/

#ifndef HistogrammingMuon_h
#define HistogrammingMuon_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMuon : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMuon(EventContainer *obj, TString muonType);
  
  // Destructor
  ~HistogrammingMuon();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMuon"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of muons

  myTH1F* _hPtObj1;      // Muon 1 PT
  myTH1F* _hEtaObj1;     // Muon 1 Eta
  myTH1F* _hPhiObj1;     // Muon 1 Phi
  myTH1F* _hRelIsoObj1;  // Muon 1 relIso
  myTH2F* _hPtEtaObj1;   // Muon 1 Pt-Eta
  myTH2F* _hPtPhiObj1;   // Muon 1 Pt-Eta
  myTH2F* _hEtaPhiObj1;  // Muon 1 Pt-Eta

  myTH1F* _hPtObj2;      // Muon 2 PT
  myTH1F* _hEtaObj2;     // Muon 2 Eta
  myTH1F* _hPhiObj2;     // Muon 2 Phi
  myTH1F* _hRelIsoObj2;  // Muon 2 relIso
  myTH2F* _hPtEtaObj2;   // Muon 2 Pt-Eta
  myTH2F* _hPtPhiObj2;   // Muon 2 Pt-Eta
  myTH2F* _hEtaPhiObj2;  // Muon 2 Pt-Eta

  myTH1F* _hRPCEtaObj1; // Muon 1 RPC Eta
  myTH1F* _hRPCPhiObj1; // Muon 1 RPC Phi
  myTH1F* _hTGCEtaObj1; // Muon 1 TGC Eta
  myTH1F* _hTGCPhiObj1; // Muon 1 TGC Phi

  myTH1F* _hRPCEtaObj2; // Muon 2 RPC Eta
  myTH1F* _hRPCPhiObj2; // Muon 2 RPC Phi
  myTH1F* _hTGCEtaObj2; // Muon 2 TGC Eta
  myTH1F* _hTGCPhiObj2; // Muon 2 TGC Phi

  myTH1F *_hPt;         // Pt of all muons
  myTH1F *_hEta;        // Eta of all muons
  myTH1F *_hPhi;        // Phi of all muons
  myTH1F *_hRelIso;     // relIso of all muons
  myTH2F* _hPtEta;      // Muon all Pt-Eta
  myTH2F* _hPtPhi;      // Muon all Pt-Eta
  myTH2F* _hEtaPhi;     // Muon all Pt-Eta

  myTH1F *_hCharge;     // Charge of all muons 
  myTH1F *_hMuonInJet;  // MuonInjet
  myTH1F *_hEtCone20;     // calorimeter ET in a cone for all muons 
  myTH1F *_hEtCone30;     // calorimeter ET in a cone for all muons 
  myTH1F *_hPtCone30;     // inner dectector track sum Pt in a cone for all muons 
  myTH1F *_hEtCone30OverPt;     // calorimeter ET in a cone 30 over Pt for all muons 
  myTH1F *_hPtCone30OverPt;     // inner dectector track sum Pt in a cone 30 over Pt for all muons 

  myTH1F* _hRPCEta; // Muon All RPC Eta
  myTH1F* _hRPCPhi; // Muon All RPC Phi
  myTH1F* _hTGCEta; // Muon All TGC Eta
  myTH1F* _hTGCPhi; // Muon All TGC Phi

 myTH1F *_hJetDeltaPhi;
 myTH1F *_hJetDeltaR;

  TString muonType; // Type of muon being plotted

};


#endif







