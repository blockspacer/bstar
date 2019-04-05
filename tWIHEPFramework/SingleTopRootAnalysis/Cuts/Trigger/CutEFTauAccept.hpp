/******************************************************************************
 * CutEFTauAccept.hpp                                                         *
 *                                                                            *
 * Cuts on tau triggers                                                       *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTauPt1 class                                 *
 *    CutEFTauAccept()                    -- Parameterized Constructor        *
 *    ~CutEFTauAccept()                   -- Destructor                       *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutEFTauAccept"              *
 *                                                                            *
 * Private Data Members of CutEFTauAccept class                               *
 *    myTH1F* _hEFTauBefore;        -- Hist before trigger cut                *
 *    myTH1F* _hEFTauAfter;         -- Hist after trigger cut                 *
 *                                                                            *
 *    Double_t _EFTauName;          -- Name of trigger cut                    *
 *                                                                            *
 * History                                                                    *
 *      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                     *
 *****************************************************************************/

#ifndef CutEFTauAccept_h
#define CutEFTauAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutEFTauAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutEFTauAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutEFTauAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutEFTauAccept_" + _EFTauName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hEFTauBefore;      // Histogram Tau Before Trigger Cut
  myTH1F* _hEFTauAfter;      // Histogram Tau Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _EFTauName;         // Tau trigger name
  
};


#endif
