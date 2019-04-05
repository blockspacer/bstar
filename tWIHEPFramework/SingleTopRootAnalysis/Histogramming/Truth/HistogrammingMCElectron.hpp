/******************************************************************************
 * HistogrammingMCElectron.hpp                                                *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMCElectron()              -- Parameterized Constructor     *
 *    ~HistogrammingMCElectron()             -- Destructor                    *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMCElectron"  *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of MCElectron properties                             *
 *                                                                            *
 * History                                                                    *
 *      07 Feb 2007 - Created by J. Holzbauer to plot MCMuon Histograms       *
 *      25 Apr 2007 - Modified by J. Holzbauer for DeltaR, PDG ID plots       *
 *****************************************************************************/

#ifndef HistogrammingMCElectron_h
#define HistogrammingMCElectron_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"



class HistogrammingMCElectron : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMCElectron(EventContainer *obj);
  
  // Destructor
  ~HistogrammingMCElectron();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMCElectron"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of MCElectrons
  myTH1F* _hPtObj1;   // MCElectron 1 PT
  myTH1F* _hPtBlowupObj1;   // MCElectron 1 PT, region 0-100GeV.
  myTH1F* _hPtObj2;   // MCElectron 2 PT
  myTH1F* _hEtaObj1;  // MCElectron 1 Eta
  myTH1F* _hEtaObj2;  // MCElectron 2 Eta
  myTH1F* _hPhiObj1;  // MCElectron 1 Phi
  myTH1F* _hPhiObj2;  // MCElectron 2 Phi

  myTH1F *_hPt;       // Pt of all MCElectrons
  myTH1F *_hEta;      // Eta of all MCElectrons
  myTH1F *_hPhi;      // Phi of all MCElectrons
  myTH1F *_hCharge;   // Charge of all MCElectrons 
  myTH2F *_hPtEta;   //pt vs eta of all electons.
  myTH1F *_hElPdgId;   // PDG ID of all MCElectrons 
  myTH1F *_hdeltaRtree;  //DeltaR between MC and Reco Muons, from tree


};


#endif






