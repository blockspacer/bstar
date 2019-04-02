/******************************************************************************
 * CutJetPt2.hpp                                                              *
 *                                                                            *
 * Cuts on pT of 2nd Jet                                                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetPt2 class                                 *
 *    CutJetPt2()                       -- Parameterized Constructor          *
 *    ~CutJetPt2()                      -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetPt2"                *
 *                                                                            *
 * Private Data Members of CutJetPt2                                          *
 *    myTH1F* _hJetPt2Before            -- Hist of Num of jets before cuts    *
 *    myTH1F* _hJetPt2After             -- Hist of Num of jets before cuts    *
 *    Double_t _JetPt2Min               -- Minimum Jet 2 Pt                   *
 *    Double_t _JetPt2Max               -- Maximum Jet 2 Pt                   *
 *                                                                            *
 * History                                                                    *
 *      14 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutJetPt2_h
#define CutJetPt2_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutJetPt2 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutJetPt2(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutJetPt2();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutJetPt2"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetPt2Before;    // Histogram Pt of jet 2 before the cut
  myTH1F* _hJetPt2After;     // Histogram Pt of jet 2 after the cut

  // Cut parameters
  Double_t _JetPt2Min;       // Minimum Pt of jet 2
  Double_t _JetPt2Max;       // Maximum Pt of jet 2

};


#endif







