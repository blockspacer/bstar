/******************************************************************************
 * CutIsolatedMuonEtaAll.hpp                                                          *
 *                                                                            *
 * Apply eta cut to ALL muons in event (not just 1st or 2nd in pT)            *
 * Select only events in which muon falls within eta region specified         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutIsolatedMuonEtaAll class                             *
 *    CutIsolatedMuonEtaAll()               -- Parameterized Constructor              *
 *    ~CutIsolatedMuonEtaAll()              -- Destructor                             *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutIsolatedMuonEtaAll"            *
 *                                                                            *
 * Private Data Members of CutIsolatedMuonEtaAll class                                *
 *    myTH1F* _hMuonEtaAllBefore;   -- Hist e Eta before cut                  *
 *    myTH1F* _hMuonEtaAllAfter;    -- Hist e Eta of jets after cut           *
 *                                                                            *
 *    Double_t _MuonEtaAllMin;      -- Minimum Muon Eta                       *
 *    Double_t _MuonEtaAllMax;      -- Maximum Muon Eta                       *
 *                                                                            *
 * History                                                                    *
 *      22 Aug 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutIsolatedMuonEtaAll_h
#define CutIsolatedMuonEtaAll_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutIsolatedMuonEtaAll : public HistoCut 
{

public:

  // Parameterized Constructor
  CutIsolatedMuonEtaAll(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutIsolatedMuonEtaAll();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutIsolatedMuonEtaAll"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hMuonEtaAllBefore;    // Histogram Muon 1 Eta before the cut
  myTH1F* _hMuonEtaAllAfter;     // Histogram Muon 1 Eta of jets after the cut

  // Cut parameters
  Int_t _MuonNumberMin;          // Minimum Muon Number
  Int_t _MuonNumberMax;          // Maximum Muon Number

  Double_t _MuonEtaAllMin;       // Minimum Muon 1 Eta
  Double_t _MuonEtaAllMax;       // Maximum Muon 1 Eta

};


#endif






