/******************************************************************************
 * CutEFMuonAccept.hpp                                                        *
 *                                                                            *
 * Cuts on muon triggers at EF level                                          *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMuonPt1 class                                *
 *    CutEFMuonAccept()                    -- Parameterized Constructor       *
 *    ~CutEFMuonAccept()                   -- Destructor                      *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutEFMuonAccept"             *
 *                                                                            *
 * Private Data Members of CutEFMuonAccept class                              *
 *    myTH1F* _hL1MuonBefore;        -- Hist before trigger cut               *
 *    myTH1F* _hL1MuonAfter;         -- Hist after trigger cut                *
 *                                                                            *
 *    Double_t _L1MuonName;          -- Name of trigger cut                   *
 *                                                                            *
 * History                                                                    *
 *      20 Jun 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutEFMuonAccept_h
#define CutEFMuonAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutEFMuonAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutEFMuonAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutEFMuonAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutEFMuonAccept_" + _EFMuonName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hEFMuonBefore;      // Histogram Muon Before Trigger Cut
  myTH1F* _hEFMuonAfter;      // Histogram Muon Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _EFMuonName;         // Muon trigger name
  
};


#endif
