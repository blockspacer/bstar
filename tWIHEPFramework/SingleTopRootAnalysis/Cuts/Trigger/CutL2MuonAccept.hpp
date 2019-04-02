/******************************************************************************
 * CutL2MuonAccept.hpp                                                        *
 *                                                                            *
 * Cuts on muon triggers at 2nd level                                         *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMuonPt1 class                                *
 *    CutL2MuonAccept()                    -- Parameterized Constructor       *
 *    ~CutL2MuonAccept()                   -- Destructor                      *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutL2MuonAccept"             *
 *                                                                            *
 * Private Data Members of CutL2MuonAccept class                              *
 *    myTH1F* _hL1MuonBefore;        -- Hist before trigger cut               *
 *    myTH1F* _hL1MuonAfter;         -- Hist after trigger cut                *
 *                                                                            *
 *    Double_t _L1MuonName;          -- Name of trigger cut                   *
 *                                                                            *
 * History                                                                    *
 *      14 May 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutL2MuonAccept_h
#define CutL2MuonAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutL2MuonAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutL2MuonAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutL2MuonAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutL2MuonAccept_" + _L2MuonName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hL2MuonBefore;      // Histogram Muon Before Trigger Cut
  myTH1F* _hL2MuonAfter;      // Histogram Muon Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _L2MuonName;         // Muon trigger name
  
};


#endif
