/******************************************************************************
 * HistogrammingTau.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTau()              -- Parameterized Constructor            *
 *    ~HistogrammingTau()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTau"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *   - lots of histograms of tau properties                                   *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *****************************************************************************/

#ifndef HistogrammingTau_h
#define HistogrammingTau_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingTau : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingTau(EventContainer *obj);
  
  // Destructor
  ~HistogrammingTau();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingTau"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of taus
  myTH1F* _hPtObj1;   // Tau 1 PT
  myTH1F* _hPtObj2;   // Tau 2 PT
  myTH1F* _hEtaObj1;  // Tau 1 Eta
  myTH1F* _hEtaObj2;  // Tau 2 Eta
  myTH1F* _hPhiObj1;  // Tau 1 Phi
  myTH1F* _hPhiObj2;  // Tau 2 Phi

  myTH1F *_hPt;       // Pt of all taus
  myTH1F *_hEta;      // Eta of all taus
  myTH1F *_hPhi;      // Phi of all taus
  myTH1F *_hCharge;   // Charge of all taus 
  myTH1F *_hNumTracks;   // number of tracks for all taus 



};


#endif







