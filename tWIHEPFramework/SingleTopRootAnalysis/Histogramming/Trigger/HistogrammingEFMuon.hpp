/******************************************************************************
 * HistogrammingEFMuon.hpp                                                    *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingEFMuon()              -- Parameterized Constructor         *
 *    ~HistogrammingEFMuon()             -- Destructor                        *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingEFMuon"      *
 *                                                                            *
 * Private Data Members of HistogrammingEFMuon class                          *
 * - lots of histograms of muon propreties                                    *
 *                                                                            *
 * History                                                                    *
 *      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                     *
 *****************************************************************************/

#ifndef HistogrammingEFMuon_h
#define HistogrammingEFMuon_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingEFMuon : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingEFMuon(EventContainer *obj);
  
  // Destructor
  ~HistogrammingEFMuon();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingEFMuon"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:  
  ////////////////////////////////////////////////////////////////
  // Number of Objects
  ////////////////////////////////////////////////////////////////
  myTH1F* _hNObj;     // Histogram of number of muons

  ////////////////////////////////////////////////////////////////
  // Kinematic Quantities - Individual
  ////////////////////////////////////////////////////////////////
  myTH1F* _hPtObj1;   // Muon 1 PT
  myTH1F* _hEtaObj1;  // Muon 1 Eta
  myTH1F* _hPhiObj1;  // Muon 1 Phi

  myTH1F* _hPtObj2;   // Muon 2 PT
  myTH1F* _hEtaObj2;  // Muon 2 Eta
  myTH1F* _hPhiObj2;  // Muon 2 Phi

  myTH1F* _hPtObj3;   // Muon 3 PT
  myTH1F* _hEtaObj3;  // Muon 3 Eta
  myTH1F* _hPhiObj3;  // Muon 3 Phi

  myTH1F* _hPtObj4;   // Muon 4 PT
  myTH1F* _hEtaObj4;  // Muon 4 Eta
  myTH1F* _hPhiObj4;  // Muon 4 Phi

  myTH1F* _hPtObj5;   // Muon 5 PT
  myTH1F* _hEtaObj5;  // Muon 5 Eta
  myTH1F* _hPhiObj5;  // Muon 5 Phi

  myTH1F* _hPtObj6;   // Muon 6 PT
  myTH1F* _hEtaObj6;  // Muon 6 Eta
  myTH1F* _hPhiObj6;  // Muon 6 Phi

  myTH1F* _hPtObj7;   // Muon 7 PT
  myTH1F* _hEtaObj7;  // Muon 7 Eta
  myTH1F* _hPhiObj7;  // Muon 7 Phi

  myTH1F* _hPtObj8;   // Muon 8 PT
  myTH1F* _hEtaObj8;  // Muon 8 Eta
  myTH1F* _hPhiObj8;  // Muon 8 Phi

  ////////////////////////////////////////////////////////////////
  // Delta R - Individual 
  ////////////////////////////////////////////////////////////////
  myTH1F *_hDeltaRMuonObj1;             // Muon 1 DeltaR between EF and Reconstructed Muons
  myTH1F *_hDeltaRMCMuonObj1;           // Muon 1 DeltaR between EF and MC muons
  myTH1F *_hDeltaRLevel1Obj1;           // Muon 1 DeltaR between EF and EF muons
  myTH1F *_hDeltaRLevel2Obj1;        // Muon 1 DeltaR between EF and Level2 muons
  myTH1F *_hDeltaRJetObj1;              // Muon 1 DeltaR between EF and Recon jets

  myTH1F *_hDeltaRMuonObj2;             // Muon 2 DeltaR between EF and Reconstructed Muons
  myTH1F *_hDeltaRMCMuonObj2;           // Muon 2 DeltaR between EF and MC muons
  myTH1F *_hDeltaRLevel1Obj2;           // Muon 2 DeltaR between EF and EF muons
  myTH1F *_hDeltaRLevel2Obj2;        // Muon 2 DeltaR between EF and Level2 muons
  myTH1F *_hDeltaRJetObj2;              // Muon 2 DeltaR between EF and Recon jets

  myTH1F *_hDeltaRMuonObj3;             // Muon 3 DeltaR between EF and Reconstructed Muons
  myTH1F *_hDeltaRMCMuonObj3;           // Muon 3 DeltaR between EF and MC muons
  myTH1F *_hDeltaRLevel1Obj3;           // Muon 3 DeltaR between EF and EF muons
  myTH1F *_hDeltaRLevel2Obj3;        // Muon 3 DeltaR between EF and Level2 muons
  myTH1F *_hDeltaRJetObj3;              // Muon 3 DeltaR between EF and Recon jets

  myTH1F *_hDeltaRMuonObj4;             // Muon 4 DeltaR between EF and Reconstructed Muons
  myTH1F *_hDeltaRMCMuonObj4;           // Muon 4 DeltaR between EF and MC muons
  myTH1F *_hDeltaRLevel1Obj4;           // Muon 4 DeltaR between EF and EF muons
  myTH1F *_hDeltaRLevel2Obj4;        // Muon 4 DeltaR between EF and Level2 muons
  myTH1F *_hDeltaRJetObj4;              // Muon 4 DeltaR between EF and Recon jets

  myTH1F *_hDeltaRMuonObj5;             // Muon 5 DeltaR between EF and Reconstructed Muons
  myTH1F *_hDeltaRMCMuonObj5;           // Muon 5 DeltaR between EF and MC muons
  myTH1F *_hDeltaRLevel1Obj5;           // Muon 5 DeltaR between EF and EF muons
  myTH1F *_hDeltaRLevel2Obj5;        // Muon 5 DeltaR between EF and Level2 muons
  myTH1F *_hDeltaRJetObj5;              // Muon 5 DeltaR between EF and Recon jets

  myTH1F *_hDeltaRMuonObj6;             // Muon 6 DeltaR between EF and Reconstructed Muons
  myTH1F *_hDeltaRMCMuonObj6;           // Muon 6 DeltaR between EF and MC muons
  myTH1F *_hDeltaRLevel1Obj6;           // Muon 6 DeltaR between EF and EF muons
  myTH1F *_hDeltaRLevel2Obj6;        // Muon 6 DeltaR between EF and Level2 muons
  myTH1F *_hDeltaRJetObj6;              // Muon 6 DeltaR between EF and Recon jets

  myTH1F *_hDeltaRMuonObj7;             // Muon 7 DeltaR between EF and Reconstructed Muons
  myTH1F *_hDeltaRMCMuonObj7;           // Muon 7 DeltaR between EF and MC muons
  myTH1F *_hDeltaRLevel1Obj7;           // Muon 7 DeltaR between EF and EF muons
  myTH1F *_hDeltaRLevel2Obj7;        // Muon 7 DeltaR between EF and Level2 muons
  myTH1F *_hDeltaRJetObj7;              // Muon 7 DeltaR between EF and Recon jets

  myTH1F *_hDeltaRMuonObj8;             // Muon 8 DeltaR between EF and Reconstructed Muons
  myTH1F *_hDeltaRMCMuonObj8;           // Muon 8 DeltaR between EF and MC muons
  myTH1F *_hDeltaRLevel1Obj8;           // Muon 8 DeltaR between EF and EF muons
  myTH1F *_hDeltaRLevel2Obj8;        // Muon 8 DeltaR between EF and Level2 muons
  myTH1F *_hDeltaRJetObj8;              // Muon 8 DeltaR between EF and Recon jets

  ////////////////////////////////////////////////////////////////
  // Eta-Phi Space - Individual
  ////////////////////////////////////////////////////////////////
  myTH2F *_hEtaPhiObj1;          // Muon 1 2-d eta-phi plot
  myTH2F *_hEtaPhiObj2;          // Muon 2 2-d eta-phi plot
  myTH2F *_hEtaPhiObj3;          // Muon 3 2-d eta-phi plot
  myTH2F *_hEtaPhiObj4;          // Muon 4 2-d eta-phi plot
  myTH2F *_hEtaPhiObj5;          // Muon 5 2-d eta-phi plot
  myTH2F *_hEtaPhiObj6;          // Muon 6 2-d eta-phi plot
  myTH2F *_hEtaPhiObj7;          // Muon 7 2-d eta-phi plot
  myTH2F *_hEtaPhiObj8;          // Muon 8 2-d eta-phi plot


  ////////////////////////////////////////////////////////////////
  // Kinematic Qunatities - All
  ////////////////////////////////////////////////////////////////
  myTH1F *_hPt;       // Pt  of all muons
  myTH1F *_hEta;      // Eta of all muons
  myTH1F *_hPhi;      // Phi of all muons
  myTH1F *_hCharge;   // Charge of all muons 

  ////////////////////////////////////////////////////////////////
  // Delta R - All
  ////////////////////////////////////////////////////////////////
  myTH1F *_hDeltaRMuon;             // DeltaR between EF and Reconstructed Muons
  myTH1F *_hDeltaRMCMuon;           // DeltaR between EF and MC muons
  myTH1F *_hDeltaRLevel1;           // DeltaR between EF and Recon muons
  myTH1F *_hDeltaRLevel2;           // DeltaR between EF and Level2 muons
  myTH1F *_hDeltaRJet;              // DeltaR between EF and Recon jets

  ////////////////////////////////////////////////////////////////
  // Eta-Phi Space - All
  ////////////////////////////////////////////////////////////////
  myTH2F *_hEtaPhi;          // 2-d eta-phi plot


};


#endif
