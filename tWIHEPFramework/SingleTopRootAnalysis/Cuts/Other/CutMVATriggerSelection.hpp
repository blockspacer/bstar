/******************************************************************************
 * CutMVATriggerSelection.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMVATriggerSelection class                              *
 *    CutMVATriggerSelection()                    -- Parameterized Constructor          *
 *    ~CutMVATriggerSelection()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutMVATriggerSelection"             *
 *                                                                            *
 * Private Data Members of CutMVATriggerSelection                                       *
 *    myTH1F* _hMVATriggerSelectionBefore         -- Hist of MVATriggerSelection before cuts      *
 *    myTH1F* _hMVATriggerSelectionAfter          -- Hist of MVATriggerSelection after cuts       *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutMVATriggerSelection_h
#define CutMVATriggerSelection_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutMVATriggerSelection : public HistoCut 
{

public:

  // Parameterized Constructor
  CutMVATriggerSelection(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutMVATriggerSelection();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutMVATriggerSelection"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hMVATriggerSelectionBefore;    // Histogram Missing Et
  myTH1F* _hMVATriggerSelectionAfter;     // Histogram Missing Et

};


#endif







