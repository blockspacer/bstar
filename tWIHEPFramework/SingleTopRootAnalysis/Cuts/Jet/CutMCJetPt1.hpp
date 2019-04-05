/******************************************************************************
 * CutMCJetPt1.hpp                                                            *
 *                                                                            *
 * Cuts on pT of 1st True Jet                                                 *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMCJetPt1 class                               *
 *    CutMCJetPt1()                       -- Parameterized Constructor        *
 *    ~CutMCJetPt1()                      -- Destructor                       *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutMCJetPt1"              *
 *                                                                            *
 * Private Data Members of CutMCJetPt1                                        *
 *    myTH1F* _hJetPt1Before            -- Hist of Num of jets before cuts    *
 *    myTH1F* _hJetPt1After             -- Hist of Num of jets before cuts    *
 *    Double_t _JetPt1Min               -- Minimum Jet 1 Pt                   *
 *    Double_t _JetPt1Max               -- Maximum Jet 1 Pt                   *
 *                                                                            *
 * History                                                                    *
 *      26 Mar 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutMCJetPt1_h
#define CutMCJetPt1_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutMCJetPt1 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutMCJetPt1(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutMCJetPt1();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutMCJetPt1"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetPt1Before;    // Histogram Pt of jet 1 before the cut
  myTH1F* _hJetPt1After;     // Histogram Pt of jet 1 after the cut

  // Cut parameters
  Double_t _JetPt1Min;       // Minimum Pt of jet 1
  Double_t _JetPt1Max;       // Maximum Pt of jet 1

};


#endif







