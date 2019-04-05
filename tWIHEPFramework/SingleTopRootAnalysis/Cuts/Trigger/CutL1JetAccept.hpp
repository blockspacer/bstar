/******************************************************************************
 * CutL1JetAccept.hpp                                                         *
 *                                                                            *
 * Cuts on jet triggers                                                       *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetPt1 class                                 *
 *    CutL1JetAccept()                    -- Parameterized Constructor        *
 *    ~CutL1JetAccept()                   -- Destructor                       * 
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutL1JetAccept"              *
 *                                                                            *
 * Private Data Members of CutL1JetAccept class                               *
 *    myTH1F* _hL1JetBefore;        -- Hist before trigger cut                *
 *    myTH1F* _hL1JetAfter;         -- Hist after trigger cut                 *
 *                                                                            *
 *    Double_t _L1JetName;          -- Name of trigger cut                    *
 *                                                                            *
 * History                                                                    *
*      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                      *
 *****************************************************************************/

#ifndef CutL1JetAccept_h
#define CutL1JetAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutL1JetAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutL1JetAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutL1JetAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutL1JetAccept_" + _L1JetName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hL1JetBefore;      // Histogram Jet Before Trigger Cut
  myTH1F* _hL1JetAfter;      // Histogram Jet Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _L1JetName;         // Jet trigger name
  
};


#endif
