/******************************************************************************
 * CutL2ElectronAccept.hpp                                                    *
 *                                                                            *
 * Cuts on electron triggers                                                  *
 * Cut can be listed in config file or passed to constructor                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutElectronPt1 class                            *
 *    CutL2ElectronAccept()                    -- Parameterized Constructor   *
 *    ~CutL2ElectronAccept()                   -- Destructor                  *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutL2ElectronAccept"         *
 *                                                                            *
 * Private Data Members of CutL2ElectronAccept class                          *
 *    myTH1F* _hL2ElectronBefore;        -- Hist before trigger cut           *
 *    myTH1F* _hL2ElectronAfter;         -- Hist after trigger cut            *
 *                                                                            *
 *    Double_t _L2ElectronName;          -- Name of trigger cut               *
 *                                                                            *
 * History                                                                    *
 *      25 Jun 2007 - Created by J. Holzbauer and P. Ryan                     *
 *****************************************************************************/

#ifndef CutL2ElectronAccept_h
#define CutL2ElectronAccept_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutL2ElectronAccept : public HistoCut 
{

public:

  // Parameterized Constructor
  CutL2ElectronAccept(EventContainer *EventContainerObj, std::string TriggerNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutL2ElectronAccept();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutL2ElectronAccept_" + _L2ElectronName; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hL2ElectronBefore;      // Histogram Electron Before Trigger Cut
  myTH1F* _hL2ElectronAfter;      // Histogram Electron Before Trigger Cut

  // Tmp
  Int_t tmpCounter;

  // Cut parameters
  std::string _L2ElectronName;         // Electron trigger name
  
};


#endif
