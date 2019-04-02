/******************************************************************************
 * CutIsolatedMuonPt1.hpp                                                     *
 *                                                                            *
 * Cuts on pT of 1st Isolated Muon                                            *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutIsolatedMuonPt1 class                        *
 *    CutIsolatedMuonPt1()                 -- Parameterized Constructor       *
 *    ~CutIsolatedMuonPt1()                -- Destructor                      *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutIsolatedMuonPt1"       *
 *                                                                            *
 * Private Data Members of CutIsolatedMuonPt1                                 *
 *    myTH1F* _hMuonPt1Before            -- Hist of Num of muons before cuts  *
 *    myTH1F* _hMuonPt1After             -- Hist of Num of muons before cuts  *
 *    Double_t _MuonPt1Min               -- Minimum Muon 1 Pt                 *
 *    Double_t _MuonPt1Max               -- Maximum Muon 1 Pt                 *
 *                                                                            *
 * History                                                                    *
 *      17 Jul 2007 - created by J. Holzbauer based on CutMuonPt1             *
 *****************************************************************************/

#ifndef CutIsolatedMuonPt1_h
#define CutIsolatedMuonPt1_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutIsolatedMuonPt1 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutIsolatedMuonPt1(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutIsolatedMuonPt1();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutIsolatedMuonPt1"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hMuonPt1Before;    // Histogram Pt of muon 1 before the cut
  myTH1F* _hMuonPt1After;     // Histogram Pt of muon 1 after the cut

  // Cut parameters
  Double_t _MuonPt1Min;       // Minimum Pt of muon 1
  Double_t _MuonPt1Max;       // Maximum Pt of muon 1

};


#endif






