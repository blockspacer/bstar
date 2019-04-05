/******************************************************************************
 * CutL1MuonAccept.hpp                                                        *
 *                                                                            *
 * Cuts on muon triggers                                                      *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMuonPt1 class                                *
 *    CutL1MuonAccept()                    -- Parameterized Constructor       *
 *    ~CutL1MuonAccept()                   -- Destructor                      *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutL1MuonAccept"             *
 *                                                                            *
 * Private Data Members of CutL1MuonAccept class                              *
 *    myTH1F* _hL1MuonBefore;        -- Hist before trigger cut               *
 *    myTH1F* _hL1MuonAfter;         -- Hist after trigger cut                *
 *                                                                            *
 *    Double_t _L1MuonName;          -- Name of trigger cut                   *
 *                                                                            *
 * History                                                                    *
 *      11 Apr 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutL1MuonAccept_h
#define CutL1MuonAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutL1MuonAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutL1MuonAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutL1MuonAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutL1MuonAccept_" + _L1MuonName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hL1MuonBefore;      // Histogram Muon Before Trigger Cut
  myTH1F* _hL1MuonAfter;      // Histogram Muon Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _L1MuonName;         // Muon trigger name
  
};


#endif
