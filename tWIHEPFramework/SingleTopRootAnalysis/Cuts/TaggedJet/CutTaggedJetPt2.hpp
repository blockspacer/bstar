/******************************************************************************
 * CutTaggedJetPt2.hpp                                                        *
 *                                                                            *
 * Cuts on pT of 2st Tagged Jet                                               *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTaggedJetPt2 class                           *
 *    CutTaggedJetPt2()                 -- Parameterized Constructor          *
 *    ~CutTaggedJetPt2()                -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutTaggedJetPt2"          *
 *                                                                            *
 * Private Data Members of CutTaggedJetPt2                                    *
 *    myTH1F* _hJetPt2Before            -- Hist of Num of jets before cuts    *
 *    myTH1F* _hJetPt2After             -- Hist of Num of jets before cuts    *
 *    Double_t _JetPt2Min               -- Minimum Jet 2 Pt                   *
 *    Double_t _JetPt2Max               -- Maximum Jet 2 Pt                   *
 *                                                                            *
 * History                                                                    *
 *      31 Jul 2007- Created by J. Holzbauer, passes events with <1 TaggedJets*
 *****************************************************************************/

#ifndef CutTaggedJetPt2_h
#define CutTaggedJetPt2_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutTaggedJetPt2 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutTaggedJetPt2(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutTaggedJetPt2();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutTaggedJetPt2"; };

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






