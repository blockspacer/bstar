/******************************************************************************
 * HistogrammingL2Electron.hpp                                                *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingL2Electron()              -- Parameterized Constructor     *
 *    ~HistogrammingL2Electron()             -- Destructor                    *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingL2Electron"  *
 *                                                                            *
 * Private Data Members of HistogrammingL2Electron class                      *
 * - lots of histograms of electron propreties                                *
 *                                                                            *
 * History                                                                    *
 *      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                     *
 *****************************************************************************/

#ifndef HistogrammingL2Electron_h
#define HistogrammingL2Electron_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingL2Electron : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingL2Electron(EventContainer *obj);
  
  // Destructor
  ~HistogrammingL2Electron();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingL2Electron"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:  
  ////////////////////////////////////////////////////////////////
  // Number of Objects
  ////////////////////////////////////////////////////////////////
  myTH1F* _hNObj;     // Histogram of number of electrons

  ////////////////////////////////////////////////////////////////
  // Kinematic Quantities - Individual
  ////////////////////////////////////////////////////////////////
  myTH1F* _hPtObj1;   // Electron 1 PT
  myTH1F* _hEtaObj1;  // Electron 1 Eta
  myTH1F* _hPhiObj1;  // Electron 1 Phi

  myTH1F* _hPtObj2;   // Electron 2 PT
  myTH1F* _hEtaObj2;  // Electron 2 Eta
  myTH1F* _hPhiObj2;  // Electron 2 Phi

  myTH1F* _hPtObj3;   // Electron 3 PT
  myTH1F* _hEtaObj3;  // Electron 3 Eta
  myTH1F* _hPhiObj3;  // Electron 3 Phi

  myTH1F* _hPtObj4;   // Electron 4 PT
  myTH1F* _hEtaObj4;  // Electron 4 Eta
  myTH1F* _hPhiObj4;  // Electron 4 Phi

  myTH1F* _hPtObj5;   // Electron 5 PT
  myTH1F* _hEtaObj5;  // Electron 5 Eta
  myTH1F* _hPhiObj5;  // Electron 5 Phi

  myTH1F* _hPtObj6;   // Electron 6 PT
  myTH1F* _hEtaObj6;  // Electron 6 Eta
  myTH1F* _hPhiObj6;  // Electron 6 Phi

  myTH1F* _hPtObj7;   // Electron 7 PT
  myTH1F* _hEtaObj7;  // Electron 7 Eta
  myTH1F* _hPhiObj7;  // Electron 7 Phi

  myTH1F* _hPtObj8;   // Electron 8 PT
  myTH1F* _hEtaObj8;  // Electron 8 Eta
  myTH1F* _hPhiObj8;  // Electron 8 Phi

  ////////////////////////////////////////////////////////////////
  // Delta R - Individual 
  ////////////////////////////////////////////////////////////////
  myTH1F *_hDeltaRElectronObj1;             // Electron 1 DeltaR between L2 and Reconstructed Electrons
  myTH1F *_hDeltaRMCElectronObj1;           // Electron 1 DeltaR between L2 and MC electrons
  myTH1F *_hDeltaRLevel1Obj1;           // Electron 1 DeltaR between L2 and L2 electrons
  myTH1F *_hDeltaREvtFilterObj1;        // Electron 1 DeltaR between L2 and EvtFilter electrons
  myTH1F *_hDeltaRJetObj1;              // Electron 1 DeltaR between L2 and Recon jets

  myTH1F *_hDeltaRElectronObj2;             // Electron 2 DeltaR between L2 and Reconstructed Electrons
  myTH1F *_hDeltaRMCElectronObj2;           // Electron 2 DeltaR between L2 and MC electrons
  myTH1F *_hDeltaRLevel1Obj2;           // Electron 2 DeltaR between L2 and L2 electrons
  myTH1F *_hDeltaREvtFilterObj2;        // Electron 2 DeltaR between L2 and EvtFilter electrons
  myTH1F *_hDeltaRJetObj2;              // Electron 2 DeltaR between L2 and Recon jets

  myTH1F *_hDeltaRElectronObj3;             // Electron 3 DeltaR between L2 and Reconstructed Electrons
  myTH1F *_hDeltaRMCElectronObj3;           // Electron 3 DeltaR between L2 and MC electrons
  myTH1F *_hDeltaRLevel1Obj3;           // Electron 3 DeltaR between L2 and L2 electrons
  myTH1F *_hDeltaREvtFilterObj3;        // Electron 3 DeltaR between L2 and EvtFilter electrons
  myTH1F *_hDeltaRJetObj3;              // Electron 3 DeltaR between L2 and Recon jets

  myTH1F *_hDeltaRElectronObj4;             // Electron 4 DeltaR between L2 and Reconstructed Electrons
  myTH1F *_hDeltaRMCElectronObj4;           // Electron 4 DeltaR between L2 and MC electrons
  myTH1F *_hDeltaRLevel1Obj4;           // Electron 4 DeltaR between L2 and L2 electrons
  myTH1F *_hDeltaREvtFilterObj4;        // Electron 4 DeltaR between L2 and EvtFilter electrons
  myTH1F *_hDeltaRJetObj4;              // Electron 4 DeltaR between L2 and Recon jets

  myTH1F *_hDeltaRElectronObj5;             // Electron 5 DeltaR between L2 and Reconstructed Electrons
  myTH1F *_hDeltaRMCElectronObj5;           // Electron 5 DeltaR between L2 and MC electrons
  myTH1F *_hDeltaRLevel1Obj5;           // Electron 5 DeltaR between L2 and L2 electrons
  myTH1F *_hDeltaREvtFilterObj5;        // Electron 5 DeltaR between L2 and EvtFilter electrons
  myTH1F *_hDeltaRJetObj5;              // Electron 5 DeltaR between L2 and Recon jets

  myTH1F *_hDeltaRElectronObj6;             // Electron 6 DeltaR between L2 and Reconstructed Electrons
  myTH1F *_hDeltaRMCElectronObj6;           // Electron 6 DeltaR between L2 and MC electrons
  myTH1F *_hDeltaRLevel1Obj6;           // Electron 6 DeltaR between L2 and L2 electrons
  myTH1F *_hDeltaREvtFilterObj6;        // Electron 6 DeltaR between L2 and EvtFilter electrons
  myTH1F *_hDeltaRJetObj6;              // Electron 6 DeltaR between L2 and Recon jets

  myTH1F *_hDeltaRElectronObj7;             // Electron 7 DeltaR between L2 and Reconstructed Electrons
  myTH1F *_hDeltaRMCElectronObj7;           // Electron 7 DeltaR between L2 and MC electrons
  myTH1F *_hDeltaRLevel1Obj7;           // Electron 7 DeltaR between L2 and L2 electrons
  myTH1F *_hDeltaREvtFilterObj7;        // Electron 7 DeltaR between L2 and EvtFilter electrons
  myTH1F *_hDeltaRJetObj7;              // Electron 7 DeltaR between L2 and Recon jets

  myTH1F *_hDeltaRElectronObj8;             // Electron 8 DeltaR between L2 and Reconstructed Electrons
  myTH1F *_hDeltaRMCElectronObj8;           // Electron 8 DeltaR between L2 and MC electrons
  myTH1F *_hDeltaRLevel1Obj8;           // Electron 8 DeltaR between L2 and L2 electrons
  myTH1F *_hDeltaREvtFilterObj8;        // Electron 8 DeltaR between L2 and EvtFilter electrons
  myTH1F *_hDeltaRJetObj8;              // Electron 8 DeltaR between L2 and Recon jets

  ////////////////////////////////////////////////////////////////
  // Eta-Phi Space - Individual
  ////////////////////////////////////////////////////////////////
  myTH2F *_hEtaPhiObj1;          // Electron 1 2-d eta-phi plot
  myTH2F *_hEtaPhiObj2;          // Electron 2 2-d eta-phi plot
  myTH2F *_hEtaPhiObj3;          // Electron 3 2-d eta-phi plot
  myTH2F *_hEtaPhiObj4;          // Electron 4 2-d eta-phi plot
  myTH2F *_hEtaPhiObj5;          // Electron 5 2-d eta-phi plot
  myTH2F *_hEtaPhiObj6;          // Electron 6 2-d eta-phi plot
  myTH2F *_hEtaPhiObj7;          // Electron 7 2-d eta-phi plot
  myTH2F *_hEtaPhiObj8;          // Electron 8 2-d eta-phi plot

  ////////////////////////////////////////////////////////////////
  // Kinematic Qunatities - All
  ////////////////////////////////////////////////////////////////
  myTH1F *_hPt;       // Pt  of all electrons
  myTH1F *_hEta;      // Eta of all electrons
  myTH1F *_hPhi;      // Phi of all electrons
  myTH1F *_hCharge;   // Charge of all electrons 

  ////////////////////////////////////////////////////////////////
  // Delta R - All
  ////////////////////////////////////////////////////////////////
  myTH1F *_hDeltaRElectron;             // DeltaR between L2 and Reconstructed Electrons
  myTH1F *_hDeltaRMCElectron;           // DeltaR between L2 and MC electrons
  myTH1F *_hDeltaRLevel1;           // DeltaR between L2 and Recon electrons
  myTH1F *_hDeltaREvtFilter;        // DeltaR between L2 and EvtFilter electrons
  myTH1F *_hDeltaRJet;              // DeltaR between L2 and Recon jets

  ////////////////////////////////////////////////////////////////
  // Eta-Phi Space - All
  ////////////////////////////////////////////////////////////////
  myTH2F *_hEtaPhi;          // 2-d eta-phi plot


};


#endif
