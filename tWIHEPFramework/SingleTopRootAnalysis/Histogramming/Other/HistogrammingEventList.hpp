/******************************************************************************
 * HistogrammingEventList.cpp                                                 *
 *                                                                            *
 * Format write out selected Events: runNumber	EventNumber                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingEventList()              -- Parameterized Constructor      *
 *    ~HistogrammingEventList()             -- Destructor                     *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingEventList"   *
 *                                                                            *
 * Private Data Members of this class                                         *
 * - lots of histograms of muon properties                                    *
 *                                                                            *
 * History                                                                    *
 *      16 Mar 2011 - Created by Huaqiao ZHANG @ CERN                         *
 *****************************************************************************/

#ifndef HistogrammingEventList_h
#define HistogrammingEventList_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include <fstream>

class HistogrammingEventList : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingEventList(EventContainer *obj, std::string evtFile);
  
  // Destructor
  ~HistogrammingEventList();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingEventList"; };

  // Book Histograms
  void BookHistogram();

  ofstream evtListFile;

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNone;     // Histogram of nothing

};


#endif







