/******************************************************************************
 * HistogrammingNone.cpp                                                      *
 *                                                                            *
 * Used as a placeholder.  Does not histogram anything                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingNone()              -- Parameterized Constructor           *
 *    ~HistogrammingNone()             -- Destructor                          *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingNone"        *
 *                                                                            *
 * Private Data Members of this class                                         *
 * - lots of histograms of muon properties                                    *
 *                                                                            *
 * History                                                                    *
 *      30 Mar 2009 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef HistogrammingNone_h
#define HistogrammingNone_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingNone : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingNone(EventContainer *obj);
  
  // Destructor
  ~HistogrammingNone();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingNone"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNone;     // Histogram of nothing

};


#endif







