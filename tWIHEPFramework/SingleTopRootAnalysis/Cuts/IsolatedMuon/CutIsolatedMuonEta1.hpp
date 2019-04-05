/******************************************************************************
 * CutIsolatedMuonEta1.cpp                                                            *
 *                                                                            *
 * Cuts on Muon quantities                                                    *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutIsolatedMuonEta1 class                               *
 *    CutIsolatedMuonEta1()                  -- Parameterized Constructor             *
 *    ~CutIsolatedMuonEta1()                 -- Destructor                            *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutIsolatedMuonEta1"                 *
 *                                                                            *
 * Private Data Members of CutIsolatedMuonEta1 class                                  *
 *    myTH1F* _hIsolatedMuonEta1Before;      -- Hist e 1 Eta before cut               *
 *    myTH1F* _hIsolatedMuonEta1After;       -- Hist e 1 Eta of jets after cut        *
 *                                                                            *
 *    Double_t _IsolatedMuonEta1Min;         -- Minimum Muon 1 Eta                    *
 *    Double_t _IsolatedMuonEta1Max;         -- Maximum Muon 1 Eta                    *
 *                                                                            *
 * History                                                                    *
 *      16 Aug 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/

#ifndef CutIsolatedMuonEta1_h
#define CutIsolatedMuonEta1_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutIsolatedMuonEta1 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutIsolatedMuonEta1(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutIsolatedMuonEta1();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutIsolatedMuonEta1"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hMuonEta1Before;     // Histogram Muon 1 Eta before the cut
  myTH1F* _hMuonEta1After;      // Histogram Muon 1 Eta of jets after the cut

  // Cut parameters
  Int_t _MuonNumberMin;         // Minimum Muon Number
  Int_t _MuonNumberMax;         // Maximum Muon Number

  Double_t _MuonEta1Min;        // Minimum Muon 1 Eta
  Double_t _MuonEta1Max;        // Maximum Muon 1 Eta

};


#endif







