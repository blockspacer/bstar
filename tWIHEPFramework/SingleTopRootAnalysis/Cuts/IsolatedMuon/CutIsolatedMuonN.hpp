/******************************************************************************
 * CutIsolatedMuonN.hpp                                                       *
 *                                                                            *
 * Cuts on number of isolated muons in an event                               *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutIsolatedMuonN class                          *
 *    CutIsolatedMuonN()                   -- Parameterized Constructor       *
 *    ~CutIsolatedMuonN()                  -- Destructor                      *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutIsolatedMuonN"         *
 *                                                                            *
 * Private Data Members of CutIsolatedMuonN                                   *
 *    myTH1F* _hMuonNumberBefore         -- Hist of Num of muons before cuts  *
 *    myTH1F* _hMuonNumberAfter          -- Hist of Num of muons before cuts  *
 *    Int_t _MuonNumberMin;              -- Minimum number of muons to require*
 *    Int_t _MuonNumberMax;              -- Maximum number of muons allowed   *
 *                                                                            *
 * History                                                                    *
 *      17 Jul 2007 - created by J. Holzbauer based on CutMuonN               *
 *****************************************************************************/

#ifndef CutIsolatedMuonN_h
#define CutIsolatedMuonN_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutIsolatedMuonN : public HistoCut 
{

public:

  // Parameterized Constructor
  CutIsolatedMuonN(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutIsolatedMuonN();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutIsolatedMuonN"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hMuonNumberBefore;    // Histogram number of muons before the cut
  myTH1F* _hMuonNumberAfter;     // Histogram number of muons after the cut

  // Cut parameters
  Int_t _MuonNumberMin;           // Minimum number of muons to require
  Int_t _MuonNumberMax;           // Maximum number of muons allowed

};


#endif







