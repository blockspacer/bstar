/******************************************************************************
 * CutL2JetAccept.hpp                                                         *
 *                                                                            *
 * Cuts on jet triggers                                                       *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetPt1 class                                 *
 *    CutL2JetAccept()                    -- Parameterized Constructor        *
 *    ~CutL2JetAccept()                   -- Destructor                       *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutL2JetAccept"              *
 *                                                                            *
 * Private Data Members of CutL2JetAccept class                               *
 *    myTH1F* _hL2JetBefore;        -- Hist before trigger cut                *
 *    myTH1F* _hL2JetAfter;         -- Hist after trigger cut                 *
 *                                                                            *
 *    Double_t _L2JetName;          -- Name of trigger cut                    *
 *                                                                            *
 * History                                                                    *
*      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                      *
 *****************************************************************************/

#ifndef CutL2JetAccept_h
#define CutL2JetAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutL2JetAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutL2JetAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutL2JetAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutL2JetAccept_" + _L2JetName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hL2JetBefore;      // Histogram Jet Before Trigger Cut
  myTH1F* _hL2JetAfter;      // Histogram Jet Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _L2JetName;         // Jet trigger name
  
};


#endif
