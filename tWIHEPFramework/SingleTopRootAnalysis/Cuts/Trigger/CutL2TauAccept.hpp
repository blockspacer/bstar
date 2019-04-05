/******************************************************************************
 * CutL2TauAccept.hpp                                                         *
 *                                                                            *
 * Cuts on tau triggers                                                       *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTauPt1 class                                 *
 *    CutL2TauAccept()                    -- Parameterized Constructor        *
 *    ~CutL2TauAccept()                   -- Destructor                       *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutL2TauAccept"              *
 *                                                                            *
 * Private Data Members of CutL2TauAccept class                               *
 *    myTH1F* _hL2TauBefore;        -- Hist before trigger cut                *
 *    myTH1F* _hL2TauAfter;         -- Hist after trigger cut                 *
 *                                                                            *
 *    Double_t _L2TauName;          -- Name of trigger cut                    *
 *                                                                            *
 * History                                                                    *
 *      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                     *
 *****************************************************************************/

#ifndef CutL2TauAccept_h
#define CutL2TauAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutL2TauAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutL2TauAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutL2TauAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutL2TauAccept_" + _L2TauName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hL2TauBefore;      // Histogram Tau Before Trigger Cut
  myTH1F* _hL2TauAfter;      // Histogram Tau Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _L2TauName;         // Tau trigger name
  
};


#endif
