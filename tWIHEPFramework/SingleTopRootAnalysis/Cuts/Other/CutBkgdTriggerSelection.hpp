/******************************************************************************
 * CutBkgdTriggerSelection.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutBkgdTriggerSelection class                              *
 *    CutBkgdTriggerSelection()                    -- Parameterized Constructor          *
 *    ~CutBkgdTriggerSelection()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutBkgdTriggerSelection"             *
 *                                                                            *
 * Private Data Members of CutBkgdTriggerSelection                                       *
 *    myTH1F* _hBkgdTriggerSelectionBefore         -- Hist of BkgdTriggerSelection before cuts      *
 *    myTH1F* _hBkgdTriggerSelectionAfter          -- Hist of BkgdTriggerSelection after cuts       *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutBkgdTriggerSelection_h
#define CutBkgdTriggerSelection_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutBkgdTriggerSelection : public HistoCut 
{

public:

  // Parameterized Constructor
  CutBkgdTriggerSelection(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutBkgdTriggerSelection();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutBkgdTriggerSelection"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hBkgdTriggerSelectionBefore;    // Histogram Missing Et
  myTH1F* _hBkgdTriggerSelectionAfter;     // Histogram Missing Et

};


#endif







