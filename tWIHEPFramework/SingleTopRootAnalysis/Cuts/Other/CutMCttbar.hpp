/******************************************************************************
 * CutMCttbar.hpp                                                             *
 *                                                                            *
 * Select events based on ttbar decay chain and books and fills histograms    *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    CutMCttbar()              -- Parameterized Constructor                  *
 *    ~CutMCttbar()             -- Destructor                                 *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "CutMCttbar"               *
 *                                                                            *
 * Private Data Members of CutMCttbar class                                   *
 * - lots of histograms of top propreties                                     *
 *                                                                            *
 * History                                                                    *
 * 19 July 2007 - Created by R. Schwienhorst                                  *
 *****************************************************************************/

#ifndef CutMCttbar_h
#define CutMCttbar_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutMCttbar : public HistoCut 
{

public:

  // Parameterized Constructor
  CutMCttbar(EventContainer *obj,Bool_t doLepJets=kTRUE, Bool_t doDiLep=kFALSE);
  
  // Destructor
  ~CutMCttbar();

  // Get the name describing the cut
  inline std::string GetCutName() { return "CutMCttbar"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // the variables storing which events we pass
  Bool_t _doLepJets; // keep lepton+jets events
  Bool_t _doDiLep;   // keep di-lepton events

  // Histograms declarations 
  myTH1F* _hNTops;    // Histogram of number of top quarks
  myTH1F* _hNWs;      // Histogram of number of W bosons
  myTH1F* _hNLepps;   // Histogram of number of positive leptons from Ws
  myTH1F* _hNLepms;   // Histogram of number of negative leptons from Ws
  

};


#endif




