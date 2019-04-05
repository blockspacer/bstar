/******************************************************************************
 * CutL1TauAccept.hpp                                                         *
 *                                                                            *
 * Cuts on tau triggers                                                       *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTauPt1 class                                 *
 *    CutL1TauAccept()                    -- Parameterized Constructor        *
 *    ~CutL1TauAccept()                   -- Destructor                       *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutL1TauAccept"              *
 *                                                                            *
 * Private Data Members of CutL1TauAccept class                               *
 *    myTH1F* _hL1TauBefore;        -- Hist before trigger cut                *
 *    myTH1F* _hL1TauAfter;         -- Hist after trigger cut                 *
 *                                                                            *
 *    Double_t _L1TauName;          -- Name of trigger cut                    *
 *                                                                            *
 * History                                                                    *
 *      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                     *
 *****************************************************************************/

#ifndef CutL1TauAccept_h
#define CutL1TauAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutL1TauAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutL1TauAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutL1TauAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutL1TauAccept_" + _L1TauName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hL1TauBefore;      // Histogram Tau Before Trigger Cut
  myTH1F* _hL1TauAfter;      // Histogram Tau Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _L1TauName;         // Tau trigger name
  
};


#endif
