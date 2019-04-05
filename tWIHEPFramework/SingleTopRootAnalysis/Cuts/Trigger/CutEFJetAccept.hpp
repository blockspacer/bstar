/******************************************************************************
 * CutEFJetAccept.hpp                                                         *
 *                                                                            *
 * Cuts on jet triggers                                                       *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetPt1 class                                 *
 *    CutEFJetAccept()                    -- Parameterized Constructor        *
 *    ~CutEFJetAccept()                   -- Destructor                       *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutEFJetAccept"              *
 *                                                                            *
 * Private Data Members of CutEFJetAccept class                               *
 *    myTH1F* _hEFJetBefore;        -- Hist before trigger cut                *
 *    myTH1F* _hEFJetAfter;         -- Hist after trigger cut                 *
 *                                                                            *
 *    Double_t _EFJetName;          -- Name of trigger cut                    *
 *                                                                            *
 * History                                                                    *
*      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                      *
 *****************************************************************************/

#ifndef CutEFJetAccept_h
#define CutEFJetAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutEFJetAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutEFJetAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutEFJetAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutEFJetAccept_" + _EFJetName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hEFJetBefore;      // Histogram Jet Before Trigger Cut
  myTH1F* _hEFJetAfter;      // Histogram Jet Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _EFJetName;         // Jet trigger name
  
};


#endif
