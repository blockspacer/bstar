/******************************************************************************
 * CutElectronTriggerMatch.hpp                                                           *
 *                                                                            *
 * Cuts on electron N (can cut on tight or veto Electrons)               *
 * Must pass tight or veto to constructor                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutElectronTriggerMatch class                              *
 *    CutElectronTriggerMatch()                     -- Parameterized Constructor         *
 *    ~CutElectronTriggerMatch()                    -- Destructor                        *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutElectronTriggerMatch"                *
 *                                                                            *
 * Private Data Members of CutElectronTriggerMatch class                                 *
 *    myTH1F* _hElectronTriggerMatchBefore;    -- Hist mu N before cut         *
 *    myTH1F* _hElectronTriggerMatchAfter;     -- Hist mu N of jets after cut  *
 *                                                                            *
 *    Int_t _ElectronTriggerMatchMin;          -- Minimum Electron N           *
 *    Int_t _ElectronTriggerMatchMax;          -- Maximum Electron N           *
 *                                                                            *
 * History                                                                    *
 *      15 Jan 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutElectronTriggerMatch_h
#define CutElectronTriggerMatch_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutElectronTriggerMatch : public HistoCut 
{

public:

  // Parameterized Constructor
  CutElectronTriggerMatch(EventContainer *EventContainerObj, TString electronType);
  
  // Destructor
  ~CutElectronTriggerMatch();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutElectronTriggerMatch"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hElectronTriggerMatchBefore;   // Histogram Electron N before the cut
  myTH1F* _hElectronTriggerMatchAfter;    // Histogram Electron N of jets after the cut

  // Cut parameters
  Bool_t _ElectronTriggerMatch;        
 
  // Electron Type (veto, tight)
  TString electronType;

};


#endif







