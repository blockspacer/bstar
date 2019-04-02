/******************************************************************************
 * CutEFElectronAccept.hpp                                                    *
 *                                                                            *
 * Cuts on electron triggers                                                  *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutElectronPt1 class                            *
 *    CutEFElectronAccept()                    -- Parameterized Constructor   *
 *    ~CutEFElectronAccept()                   -- Destructor                  *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutEFElectronAccept"         *
 *                                                                            *
 * Private Data Members of CutEFElectronAccept class                          *
 *    myTH1F* _hEFElectronBefore;        -- Hist before trigger cut           *
 *    myTH1F* _hEFElectronAfter;         -- Hist after trigger cut            *
 *                                                                            *
 *    Double_t _EFElectronName;          -- Name of trigger cut               *
 *                                                                            *
 * History                                                                    *
 *      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                     *
 *****************************************************************************/

#ifndef CutEFElectronAccept_h
#define CutEFElectronAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutEFElectronAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutEFElectronAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutEFElectronAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutEFElectronAccept_" + _EFElectronName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hEFElectronBefore;      // Histogram Electron Before Trigger Cut
  myTH1F* _hEFElectronAfter;      // Histogram Electron Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _EFElectronName;         // Electron trigger name
  
};


#endif
