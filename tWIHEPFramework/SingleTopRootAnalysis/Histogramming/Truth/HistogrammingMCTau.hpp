/******************************************************************************
 * HistogrammingMCTau.hpp                                                     *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMCTau()              -- Parameterized Constructor          *
 *    ~HistogrammingMCTau()             -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMCTau"       *
 *                                                                            *
 * Private Data Members of this class                                         *
 *   - lots of histograms of MCTau properties                                 *
 *                                                                            *
 * History                                                                    *
 *      07 Feb 2007 - Created by J. Holzbauer to plot MCMuon Histograms       *
 *      25 Apr 2007 - Modified by J. Holzbauer for DeltaR, PDG ID plots       *
 *****************************************************************************/

#ifndef HistogrammingMCTau_h
#define HistogrammingMCTau_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMCTau : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMCTau(EventContainer *obj);
  
  // Destructor
  ~HistogrammingMCTau();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMCTau"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of MCTaus
  myTH1F* _hPtObj1;   // MCTau 1 PT
  myTH1F* _hPtObj2;   // MCTau 2 PT
  myTH1F* _hEtaObj1;  // MCTau 1 Eta
  myTH1F* _hEtaObj2;  // MCTau 2 Eta
  myTH1F* _hPhiObj1;  // MCTau 1 Phi
  myTH1F* _hPhiObj2;  // MCTau 2 Phi

  myTH1F *_hPt;       // Pt of all MCTaus
  myTH1F *_hEta;      // Eta of all MCTaus
  myTH1F *_hPhi;      // Phi of all MCTaus
  myTH1F *_hCharge;   // Charge of all MCTaus 
  myTH1F *_hTauPdgId;   // PDG ID of all MCTaus 
  myTH1F *_hdeltaRtree;  //DeltaR between MC and Reco Muons, from tree

};


#endif






