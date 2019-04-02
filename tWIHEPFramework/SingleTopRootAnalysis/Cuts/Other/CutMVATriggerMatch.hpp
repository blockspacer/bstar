/******************************************************************************
 * CutMVATriggerMatch.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMVATriggerMatch class                              *
 *    CutMVATriggerMatch()                    -- Parameterized Constructor          *
 *    ~CutMVATriggerMatch()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutMVATriggerMatch"             *
 *                                                                            *
 * Private Data Members of CutMVATriggerMatch                                       *
 *    myTH1F* _hMVATriggerMatchBefore         -- Hist of MVATriggerMatch before cuts      *
 *    myTH1F* _hMVATriggerMatchAfter          -- Hist of MVATriggerMatch after cuts       *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutMVATriggerMatch_h
#define CutMVATriggerMatch_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutMVATriggerMatch : public HistoCut 
{

public:

  // Parameterized Constructor
  CutMVATriggerMatch(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutMVATriggerMatch();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutMVATriggerMatch"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hMVATriggerMatchBefore;    // Histogram Missing Et
  myTH1F* _hMVATriggerMatchAfter;     // Histogram Missing Et

};


#endif







