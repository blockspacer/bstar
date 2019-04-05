/******************************************************************************
 * CutHFOR.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutHFOR class                              *
 *    CutHFOR()                    -- Parameterized Constructor          *
 *    ~CutHFOR()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutHFOR"             *
 *                                                                            *
 * Private Data Members of CutHFOR                                       *
 *    myTH1F* _hHFORBefore         -- Hist of HFOR before cuts      *
 *    myTH1F* _hHFORAfter          -- Hist of HFOR after cuts       *
 *                                                                            *
 * History                                                                    *
 *      10 Dec 10: Created J. Holzbauer                                       *
 *****************************************************************************/

#ifndef CutHFOR_h
#define CutHFOR_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutHFOR : public HistoCut 
{

public:

  // Parameterized Constructor
  CutHFOR(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutHFOR();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutHFOR"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hHFORBefore;
  myTH1F* _hHFORAfter; 

};


#endif







