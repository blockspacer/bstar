/******************************************************************************
 * HistogrammingMET.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMET()              -- Parameterized Constructor            *
 *    ~HistogrammingMET()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMET"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - histograms of MET                                                       *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      20 Nov 2006 - Modified by Bernard Pope                                *
 *      21 Mar 2007 - RS: Fill from event container, add sumET, mex, mey      *
 *****************************************************************************/

#ifndef HistogrammingMET_h
#define HistogrammingMET_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMET : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMET(EventContainer *obj);
  
  // Destructor
  ~HistogrammingMET();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMET"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hMET;  // Histogram of MET
  myTH1F* _hMEX;  // Histogram of MEX
  myTH1F* _hMEY;  // Histogram of MEY
  myTH1F* _hMETPhi;  // Histogram of MET phi

  //Corrected met distributions
  myTH1F* _hMET_xy; // Histogram of xy corrected MET
  myTH1F* _hMEX_xy; // Histogram of xy corrected MEX
  myTH1F* _hMEY_xy; // Histogram of xy corrected MEY
  myTH1F* _hMETPhi_xy; // Histogram of xy corrected MET phi

  myTH1F* _hSumEt;  // Histogram of SumEt

};


#endif







