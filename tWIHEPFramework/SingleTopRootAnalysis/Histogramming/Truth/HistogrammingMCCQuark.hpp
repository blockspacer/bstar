/*****************************************************************************
 * HistogrammingMCCQuark.hpp                                                   *
 *                                                                           *
 * Books and fills histograms                                                *
 * Used for events passing cuts applied in other classes                     *
 *                                                                           *
 * Derived from HistoCut which is in turn derived from BaseCut               *
 *                                                                           *
 * Public Member Functions of AnalysisMain class                             *
 *    HistogrammingMCCQuark()              -- Parameterized Constructor        *
 *    ~HistogrammingMCCQuark()             -- Destructor                       *
 *    BookHistogram()                   -- Book histograms                   *
 *    Apply()                           -- Fill histograms only (No Cuts)    *
 *    GetCutName()                      -- Returns "HistogrammingMCQuark"      *
 *                                                                           *
 * Private Data Members of HistogrammingMCQuark class                          *
 * - lots of histograms of muon propreties                                   *
 *                                                                           *
 * History                                                                   *
 *     17 Jul 2007- Created by J. Holzbauer, based on HistogrammingMCBQuark    *
 *****************************************************************************/

#ifndef HistogrammingMCCQuark_h
#define HistogrammingMCCQuark_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMCCQuark : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMCCQuark(EventContainer *obj);
  
  // Destructor
  ~HistogrammingMCCQuark();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMCCQuark"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of jets
 
  myTH1F* _hPtObj1;   // Quark 1 PT
  myTH1F* _hEtaObj1;  // Quark 1 Eta
  myTH1F* _hPhiObj1;  // Quark 1 Phi

  myTH1F* _hEtaObj2;  // Quark 2 Eta
  myTH1F* _hPtObj2;   // Quark 2 PT
  myTH1F* _hPhiObj2;  // Quark 2 Phi

  myTH1F* _hPtObj3;   // Quark 3 PT
  myTH1F* _hEtaObj3;  // Quark 3 Eta
  myTH1F* _hPhiObj3;  // Quark 3 Phi

  myTH1F* _hPtObj4;   // Quark 4 PT
  myTH1F* _hEtaObj4;  // Quark 4 Eta
  myTH1F* _hPhiObj4;  // Quark 4 Phi

  myTH1F* _hPtObj5;   // Quark 5 PT
  myTH1F* _hEtaObj5;  // Quark 5 Eta
  myTH1F* _hPhiObj5;  // Quark 5 Phi

  myTH1F* _hPtObj6;   // Quark 6 PT
  myTH1F* _hEtaObj6;  // Quark 6 Eta
  myTH1F* _hPhiObj6;  // Quark 6 Phi

  myTH1F* _hPtObj7;   // Quark 7 PT
  myTH1F* _hEtaObj7;  // Quark 7 Eta
  myTH1F* _hPhiObj7;  // Quark 7 Phi

  myTH1F *_hPt;       // Pt of all jets
  myTH1F *_hEta;      // Eta of all jets
  myTH1F *_hPhi;      // Phi of all jets
  myTH1F *_hCharge;   // Charge of all jets 
  myTH1F *_hQuarkPdgId; //PDG ID of all jets
  myTH1F *_hdeltaRtree;  //DeltaR between MC and Reco Muons, from tree

};


#endif



