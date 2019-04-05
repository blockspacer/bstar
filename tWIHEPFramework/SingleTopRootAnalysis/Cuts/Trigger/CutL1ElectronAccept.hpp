/******************************************************************************
 * CutL1ElectronAccept.hpp                                                    *
 *                                                                            *
 * NOTE: Cuts on EM triggers                                                  *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutElectronPt1 class                            *
 *    CutL1ElectronAccept()                    -- Parameterized Constructor   *
 *    ~CutL1ElectronAccept()                   -- Destructor                  *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutL1ElectronAccept"         *
 *                                                                            *
 * Private Data Members of CutL1ElectronAccept class                          *
 *    myTH1F* _hL1ElectronBefore;        -- Hist before trigger cut           *
 *    myTH1F* _hL1ElectronAfter;         -- Hist after trigger cut            *
 *                                                                            *
 *    Double_t _L1ElectronName;          -- Name of trigger cut               *
 *                                                                            *
 * History                                                                    *
 *      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                     * 
 *****************************************************************************/

#ifndef CutL1ElectronAccept_h
#define CutL1ElectronAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutL1ElectronAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutL1ElectronAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutL1ElectronAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutL1ElectronAccept_" + _L1ElectronName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hL1ElectronBefore;      // Histogram Electron Before Trigger Cut
  myTH1F* _hL1ElectronAfter;      // Histogram Electron Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _L1ElectronName;         // Electron trigger name
  
};


#endif
