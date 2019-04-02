/******************************************************************************
 * CutTriggerSelection.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTriggerSelection class                              *
 *    CutTriggerSelection()                    -- Parameterized Constructor          *
 *    ~CutTriggerSelection()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutTriggerSelection"             *
 *                                                                            *
 * Private Data Members of CutTriggerSelection                                       *
 *    myTH1F* _hTriggerSelectionBefore         -- Hist of TriggerSelection before cuts      *
 *    myTH1F* _hTriggerSelectionAfter          -- Hist of TriggerSelection after cuts       *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutTriggerSelection_h
#define CutTriggerSelection_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <sstream>

class CutTriggerSelection : public HistoCut 
{

public:

  // Parameterized Constructor
  CutTriggerSelection(EventContainer *EventContainerObj, int whichtrigger = -1);
  
  // Destructor
  ~CutTriggerSelection();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutTriggerSelection"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hTriggerSelectionBefore;    // Histogram Missing Et
  myTH1F* _hTriggerSelectionAfter;     // Histogram Missing Et
  int _whichtrigger;
  TString _triggerChannel;

};


#endif







