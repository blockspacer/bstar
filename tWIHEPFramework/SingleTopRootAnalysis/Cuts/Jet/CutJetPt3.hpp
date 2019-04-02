/******************************************************************************
 * CutJetPt3.hpp                                                              *
 *                                                                            *
 * Cuts on pT of 3rd Jet                                                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetPt3 class                                 *
 *    CutJetPt3()                       -- Parameterized Constructor          *
 *    ~CutJetPt3()                      -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetPt3"                *
 *                                                                            *
 * Private Data Members of CutJetPt3                                          *
 *    myTH1F* _hJetPt3Before            -- Hist of Num of jets before cuts    *
 *    myTH1F* _hJetPt3After             -- Hist of Num of jets before cuts    *
 *    Double_t _JetPt3Min               -- Minimum Jet 3 Pt                   *
 *    Double_t _JetPt3Max               -- Maximum Jet 3 Pt                   *
 *                                                                            *
 * History                                                                    *
 *      11 Jul 2007 -- Created by J. Holzbauer, based on CutJetPt2, by P. Ryan*
 *****************************************************************************/

#ifndef CutJetPt3_h
#define CutJetPt3_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutJetPt3 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutJetPt3(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutJetPt3();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutJetPt3"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetPt3Before;    // Histogram Pt of jet 3 before the cut
  myTH1F* _hJetPt3After;     // Histogram Pt of jet 3 after the cut

  // Cut parameters
  Double_t _JetPt3Min;       // Minimum Pt of jet 3
  Double_t _JetPt3Max;       // Maximum Pt of jet 3

};


#endif






