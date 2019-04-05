/******************************************************************************
 * HistogrammingNJets.hpp                                                     *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingNJets()              -- Parameterized Constructor          *
 *    ~HistogrammingNJets()             -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingNJets"       *
 *                                                                            *
 * Private Data Members of AnalysisMain class                                 *
 *    myTH1F* _hNJets                   -- Hist of Num of jets                *
 *                                                                            *
 * History                                                                    *
 *      19 Oct 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      10 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      13 Nov 2006 - RS - add jet PT histograms                              *
 *****************************************************************************/

#ifndef HistogrammingNJets_h
#define HistogrammingNJets_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"



class HistogrammingNJets : public HistoCut 
{

public:

  // Parameterized Donstructor
  HistogrammingNJets(EventContainer *obj);
  
  // Destructor
  ~HistogrammingNJets();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingNJets"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNJets;  // Histogram of events passing cuts
  myTH1F* _hPtJet1;  // Jet 1 PT
  myTH1F* _hPtJet2;  // Jet 2 PT
  myTH1F* _hPtJet3;  // Jet 3 PT
  myTH1F* _hPtJet4;  // Jet 4 PT

};


#endif







