/*********************************************************************************
 * HistogrammingTriggersPassed.hpp                                               *
 *                                                                               *
 * Books and fills histograms                                                    *
 * Plot histograms showing which triggers are passed                             *
 * Shows which triggers can give 100% efficiency                                 *
 *                                                                               *
 * Derived from HistoCut which is in turn derived from BaseCut                   *
 *                                                                               *
 * Public Member Functions of AnalysisMain class                                 *
 *    HistogrammingTriggersPassed()              -- Parameterized Constructor    *
 *    ~HistogrammingTriggersPassed()             -- Destructor                   *
 *    BookHistogram()                   -- Book histograms                       *
 *    Apply()                           -- Fill histograms only (No Cuts)        *
 *    GetCutName()                      -- Returns "HistogrammingTriggersPassed" *
 *                                                                               *
 * Private Data Members of this class                                            *
 * - lots of histograms                                                          *
 *                                                                               *
 * History                                                                       *
 *      11 June 2010 - Created by P. Ryan                                        *
 ********************************************************************************/

#ifndef HistogrammingTriggersPassed_h
#define HistogrammingTriggersPassed_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingTriggersPassed : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingTriggersPassed(EventContainer *obj);
  
  // Destructor
  ~HistogrammingTriggersPassed();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingTriggersPassed"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:  
  ////////////////////////////////////////////////////////////////
  // Jet Triggers
  ////////////////////////////////////////////////////////////////
  myTH1F* _hJet10;
  myTH1F* _hJet20;
  myTH1F* _hJet40;
  myTH1F* _hJet80;
  myTH1F* _hJet140;
  myTH1F* _hJet200;

  ////////////////////////////////////////////////////////////////
  // Missing ET Triggers
  ////////////////////////////////////////////////////////////////
  myTH1F* _hxe20;
  myTH1F* _hxe30unbiased;

  ////////////////////////////////////////////////////////////////
  // Jet + Missing ET Triggers (none as of now)
  ////////////////////////////////////////////////////////////////
  myTH1F* _h2Jet40;
  myTH1F* _hJet80_xe30unbiased;
  myTH1F* _h2Jet40_xe30unbiased;
  myTH1F* _hJet80_xe20;
  myTH1F* _h2Jet40_xe20;

  ////////////////////////////////////////////////////////////////
  // No trigger
  ////////////////////////////////////////////////////////////////
  myTH1F* _hNone;

  // tmp
  Int_t count_L1_J55;
  Int_t count_L2_j90;
  Int_t count_EF_j140;
  Int_t count_All;

};
  
#endif






