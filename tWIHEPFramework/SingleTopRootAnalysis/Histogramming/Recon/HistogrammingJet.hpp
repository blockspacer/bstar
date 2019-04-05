/******************************************************************************
 * HistogrammingJet.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingJet()              -- Parameterized Constructor            *
 *    ~HistogrammingJet()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingJet"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *   - lots of histograms of jet properties                                   *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *****************************************************************************/

#ifndef HistogrammingJet_h
#define HistogrammingJet_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingJet : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingJet(EventContainer *obj);
  
  // Destructor
  ~HistogrammingJet();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingJet"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of jets

  myTH1F* _hPtObj1;   // Jet 1 PT
  myTH1F* _hEtaObj1;  // Jet 1 Eta
  myTH1F* _hPhiObj1;  // Jet 1 Phi
  myTH1F* _hTagObj1;  // Jet 1 Tag

  myTH1F* _hEtaObj2;  // Jet 2 Eta
  myTH1F* _hPtObj2;   // Jet 2 PT
  myTH1F* _hPhiObj2;  // Jet 2 Phi
  myTH1F* _hTagObj2;  // Jet 2 Tag

  myTH1F* _hPtObj3;   // Jet 3 PT
  myTH1F* _hEtaObj3;  // Jet 3 Eta
  myTH1F* _hPhiObj3;  // Jet 3 Phi
  myTH1F* _hTagObj3;  // Jet 3 Tag

  myTH1F* _hPtObj4;   // Jet 4 PT
  myTH1F* _hEtaObj4;  // Jet 4 Eta
  myTH1F* _hPhiObj4;  // Jet 4 Phi
  myTH1F* _hTagObj4;  // Jet 4 Tag

  myTH1F* _hPtObj5;   // Jet 5 PT
  myTH1F* _hEtaObj5;  // Jet 5 Eta
  myTH1F* _hPhiObj5;  // Jet 5 Phi
  myTH1F* _hTagObj5;  // Jet 5 Tag

  myTH1F* _hPtObj6;   // Jet 6 PT
  myTH1F* _hEtaObj6;  // Jet 6 Eta
  myTH1F* _hPhiObj6;  // Jet 6 Phi
  myTH1F* _hTagObj6;  // Jet 6 Tag

  myTH1F* _hPtObj7;   // Jet 7 PT
  myTH1F* _hEtaObj7;  // Jet 7 Eta
  myTH1F* _hPhiObj7;  // Jet 7 Phi
  myTH1F* _hTagObj7;  // Jet 7 Tag

  myTH1F *_hPt;       // Pt of all jets
  myTH1F *_hEta;      // Eta of all jets
  myTH1F *_hPhi;      // Phi of all jets
  myTH1F *_hCharge;   // Charge of all jets 
  myTH1F *_hEtEM0;    // ET in EM layer 0 for all jets
  myTH1F *_hTruthLabel;  

 

};


#endif







