/******************************************************************************
 * CutJetPt4.hpp                                                              *
 *                                                                            *
 * Cuts on pT of 4th Jet                                                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetPt4 class                                 *
 *    CutJetPt4()                       -- Parameterized Constructor          *
 *    ~CutJetPt4()                      -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetPt4"                *
 *                                                                            *
 * Private Data Members of CutJetPt4                                          *
 *    myTH1F* _hJetPt4Before            -- Hist of Num of jets before cuts    *
 *    myTH1F* _hJetPt4After             -- Hist of Num of jets before cuts    *
 *    Double_t _JetPt4Min               -- Minimum Jet 4 Pt                   *
 *    Double_t _JetPt4Max               -- Maximum Jet 4 Pt                   *
 *                                                                            *
 * History                                                                    *
 *      11 Jul 2007 -- Created by J. Holzbauer, based on CutJetPt2, by P. Ryan*
 *****************************************************************************/

#ifndef CutJetPt4_h
#define CutJetPt4_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutJetPt4 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutJetPt4(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutJetPt4();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutJetPt4"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetPt4Before;    // Histogram Pt of jet 4 before the cut
  myTH1F* _hJetPt4After;     // Histogram Pt of jet 4 after the cut

  // Cut parameters
  Double_t _JetPt4Min;       // Minimum Pt of jet 4
  Double_t _JetPt4Max;       // Maximum Pt of jet 4

};


#endif
