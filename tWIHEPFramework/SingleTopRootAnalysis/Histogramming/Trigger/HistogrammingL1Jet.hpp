/******************************************************************************
 * HistogrammingL1Jet.hpp                                                     *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingL1Jet()              -- Parameterized Constructor          *
 *    ~HistogrammingL1Jet()             -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingL1Jet"       *
 *                                                                            *
 * Private Data Members of HistogrammingL1Jet class                           *
 * - lots of histograms of jet propreties                                     *
 *                                                                            *
 * History                                                                    *
 *      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                     *
 *****************************************************************************/

#ifndef HistogrammingL1Jet_h
#define HistogrammingL1Jet_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingL1Jet : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingL1Jet(EventContainer *obj);
  
  // Destructor
  ~HistogrammingL1Jet();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingL1Jet"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:  
  ////////////////////////////////////////////////////////////////
  // Number of Objects
  ////////////////////////////////////////////////////////////////
  myTH1F* _hNObj;     // Histogram of number of jets

  ////////////////////////////////////////////////////////////////
  // Kinematic Quantities - Individual
  ////////////////////////////////////////////////////////////////
  myTH1F* _hPtObj1;   // Jet 1 PT
  myTH1F* _hEtaObj1;  // Jet 1 Eta
  myTH1F* _hPhiObj1;  // Jet 1 Phi

  myTH1F* _hPtObj2;   // Jet 2 PT
  myTH1F* _hEtaObj2;  // Jet 2 Eta
  myTH1F* _hPhiObj2;  // Jet 2 Phi

  myTH1F* _hPtObj3;   // Jet 3 PT
  myTH1F* _hEtaObj3;  // Jet 3 Eta
  myTH1F* _hPhiObj3;  // Jet 3 Phi

  myTH1F* _hPtObj4;   // Jet 4 PT
  myTH1F* _hEtaObj4;  // Jet 4 Eta
  myTH1F* _hPhiObj4;  // Jet 4 Phi

  myTH1F* _hPtObj5;   // Jet 5 PT
  myTH1F* _hEtaObj5;  // Jet 5 Eta
  myTH1F* _hPhiObj5;  // Jet 5 Phi

  myTH1F* _hPtObj6;   // Jet 6 PT
  myTH1F* _hEtaObj6;  // Jet 6 Eta
  myTH1F* _hPhiObj6;  // Jet 6 Phi

  myTH1F* _hPtObj7;   // Jet 7 PT
  myTH1F* _hEtaObj7;  // Jet 7 Eta
  myTH1F* _hPhiObj7;  // Jet 7 Phi

  myTH1F* _hPtObj8;   // Jet 8 PT
  myTH1F* _hEtaObj8;  // Jet 8 Eta
  myTH1F* _hPhiObj8;  // Jet 8 Phi

  ////////////////////////////////////////////////////////////////
  // Delta R - Individual 
  ////////////////////////////////////////////////////////////////
  myTH1F *_hDeltaRJetObj1;             // Jet 1 DeltaR between L1 and Reconstructed Jets
  myTH1F *_hDeltaRMCJetObj1;           // Jet 1 DeltaR between L1 and MC jets
  myTH1F *_hDeltaRLevel2Obj1;           // Jet 1 DeltaR between L1 and L2 jets
  myTH1F *_hDeltaREvtFilterObj1;        // Jet 1 DeltaR between L1 and EvtFilter jets

  myTH1F *_hDeltaRJetObj2;             // Jet 2 DeltaR between L1 and Reconstructed Jets
  myTH1F *_hDeltaRMCJetObj2;           // Jet 2 DeltaR between L1 and MC jets
  myTH1F *_hDeltaRLevel2Obj2;           // Jet 2 DeltaR between L1 and L2 jets
  myTH1F *_hDeltaREvtFilterObj2;        // Jet 2 DeltaR between L1 and EvtFilter jets

  myTH1F *_hDeltaRJetObj3;             // Jet 3 DeltaR between L1 and Reconstructed Jets
  myTH1F *_hDeltaRMCJetObj3;           // Jet 3 DeltaR between L1 and MC jets
  myTH1F *_hDeltaRLevel2Obj3;           // Jet 3 DeltaR between L1 and L2 jets
  myTH1F *_hDeltaREvtFilterObj3;        // Jet 3 DeltaR between L1 and EvtFilter jets

  myTH1F *_hDeltaRJetObj4;             // Jet 4 DeltaR between L1 and Reconstructed Jets
  myTH1F *_hDeltaRMCJetObj4;           // Jet 4 DeltaR between L1 and MC jets
  myTH1F *_hDeltaRLevel2Obj4;           // Jet 4 DeltaR between L1 and L2 jets
  myTH1F *_hDeltaREvtFilterObj4;        // Jet 4 DeltaR between L1 and EvtFilter jets

  myTH1F *_hDeltaRJetObj5;             // Jet 5 DeltaR between L1 and Reconstructed Jets
  myTH1F *_hDeltaRMCJetObj5;           // Jet 5 DeltaR between L1 and MC jets
  myTH1F *_hDeltaRLevel2Obj5;           // Jet 5 DeltaR between L1 and L2 jets
  myTH1F *_hDeltaREvtFilterObj5;        // Jet 5 DeltaR between L1 and EvtFilter jets

  myTH1F *_hDeltaRJetObj6;             // Jet 6 DeltaR between L1 and Reconstructed Jets
  myTH1F *_hDeltaRMCJetObj6;           // Jet 6 DeltaR between L1 and MC jets
  myTH1F *_hDeltaRLevel2Obj6;           // Jet 6 DeltaR between L1 and L2 jets
  myTH1F *_hDeltaREvtFilterObj6;        // Jet 6 DeltaR between L1 and EvtFilter jets

  myTH1F *_hDeltaRJetObj7;             // Jet 7 DeltaR between L1 and Reconstructed Jets
  myTH1F *_hDeltaRMCJetObj7;           // Jet 7 DeltaR between L1 and MC jets
  myTH1F *_hDeltaRLevel2Obj7;           // Jet 7 DeltaR between L1 and L2 jets
  myTH1F *_hDeltaREvtFilterObj7;        // Jet 7 DeltaR between L1 and EvtFilter jets

  myTH1F *_hDeltaRJetObj8;             // Jet 8 DeltaR between L1 and Reconstructed Jets
  myTH1F *_hDeltaRMCJetObj8;           // Jet 8 DeltaR between L1 and MC jets
  myTH1F *_hDeltaRLevel2Obj8;           // Jet 8 DeltaR between L1 and L2 jets
  myTH1F *_hDeltaREvtFilterObj8;        // Jet 8 DeltaR between L1 and EvtFilter jets

  ////////////////////////////////////////////////////////////////
  // Eta-Phi Space - Individual
  ////////////////////////////////////////////////////////////////
  myTH2F *_hEtaPhiObj1;          // Jet 1 2-d eta-phi plot
  myTH2F *_hEtaPhiObj2;          // Jet 2 2-d eta-phi plot
  myTH2F *_hEtaPhiObj3;          // Jet 3 2-d eta-phi plot
  myTH2F *_hEtaPhiObj4;          // Jet 4 2-d eta-phi plot
  myTH2F *_hEtaPhiObj5;          // Jet 5 2-d eta-phi plot
  myTH2F *_hEtaPhiObj6;          // Jet 6 2-d eta-phi plot
  myTH2F *_hEtaPhiObj7;          // Jet 7 2-d eta-phi plot
  myTH2F *_hEtaPhiObj8;          // Jet 8 2-d eta-phi plot

  ////////////////////////////////////////////////////////////////
  // RoI - Individual 
  ////////////////////////////////////////////////////////////////
  // RoI
  myTH1F *_hRoIObj1;             // RoI object 1
  myTH1F *_hRoIObj2;             // RoI object 2
  myTH1F *_hRoIObj3;             // RoI object 3
  myTH1F *_hRoIObj4;             // RoI object 4
  myTH1F *_hRoIObj5;             // RoI object 5
  myTH1F *_hRoIObj6;             // RoI object 6
  myTH1F *_hRoIObj7;             // RoI object 7
  myTH1F *_hRoIObj8;             // RoI object 8

  // RoI SectorType
  myTH1F *_hRoISectorTypeObj1;       // RoI SectorType object 1
  myTH1F *_hRoISectorTypeObj2;       // RoI SectorType object 2
  myTH1F *_hRoISectorTypeObj3;       // RoI SectorType object 3
  myTH1F *_hRoISectorTypeObj4;       // RoI SectorType object 4
  myTH1F *_hRoISectorTypeObj5;       // RoI SectorType object 5
  myTH1F *_hRoISectorTypeObj6;       // RoI SectorType object 6
  myTH1F *_hRoISectorTypeObj7;       // RoI SectorType object 7
  myTH1F *_hRoISectorTypeObj8;       // RoI SectorType object 8

  ////////////////////////////////////////////////////////////////
  // Kinematic Qunatities - All
  ////////////////////////////////////////////////////////////////
  myTH1F *_hPt;       // Pt  of all jets
  myTH1F *_hEta;      // Eta of all jets
  myTH1F *_hPhi;      // Phi of all jets
  myTH1F *_hCharge;   // Charge of all jets 

  ////////////////////////////////////////////////////////////////
  // Delta R - All
  ////////////////////////////////////////////////////////////////
  myTH1F *_hDeltaRJet;             // DeltaR between L1 and Reconstructed Jets
  myTH1F *_hDeltaRMCJet;           // DeltaR between L1 and MC jets
  myTH1F *_hDeltaRLevel2;           // DeltaR between L1 and Recon jets
  myTH1F *_hDeltaREvtFilter;        // DeltaR between L1 and EvtFilter jets

  ////////////////////////////////////////////////////////////////
  // Eta-Phi Space - All
  ////////////////////////////////////////////////////////////////
  myTH2F *_hEtaPhi;          // 2-d eta-phi plot

  ////////////////////////////////////////////////////////////////
  // RoI - All
  ////////////////////////////////////////////////////////////////
  myTH1F *_hRoI;             // RoI 
  myTH1F *_hRoISectorType;   // RoI SectorType

};


#endif




