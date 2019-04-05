/******************************************************************************
 * CutBkgdTriggerMatch.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutBkgdTriggerMatch class                              *
 *    CutBkgdTriggerMatch()                    -- Parameterized Constructor          *
 *    ~CutBkgdTriggerMatch()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutBkgdTriggerMatch"             *
 *                                                                            *
 * Private Data Members of CutBkgdTriggerMatch                                       *
 *    myTH1F* _hBkgdTriggerMatchBefore         -- Hist of BkgdTriggerMatch before cuts      *
 *    myTH1F* _hBkgdTriggerMatchAfter          -- Hist of BkgdTriggerMatch after cuts       *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutBkgdTriggerMatch_h
#define CutBkgdTriggerMatch_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutBkgdTriggerMatch : public HistoCut 
{

public:

  // Parameterized Constructor
  CutBkgdTriggerMatch(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutBkgdTriggerMatch();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutBkgdTriggerMatch"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hBkgdTriggerMatchBefore;    // Histogram Missing Et
  myTH1F* _hBkgdTriggerMatchAfter;     // Histogram Missing Et

};


#endif







