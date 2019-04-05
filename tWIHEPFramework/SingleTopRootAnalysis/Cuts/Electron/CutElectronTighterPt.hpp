/******************************************************************************
 * CutElectronTighterPt.hpp                                                    *
 *                                                                            *
 * Cuts on electron Number (can cut on tight or veto Electrons)               *
 * Must pass tight or veto to constructor                                     *
 * With special/extra cut for leading/subLeading electron                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutElectronTighterPt class                       *
 *    CutElectronTighterPt()                     -- Parameterized Constructor  *
 *    ~CutElectronTighterPt()                    -- Destructor                 *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutElectronTighterPt"         *
 *                                                                            *
 * Private Data Members of CutElectronTighterPt class                          *
 *    myTH1F* _hElectronNumberBefore;    -- Hist mu Number before cut         *
 *    myTH1F* _hElectronNumberAfter;     -- Hist mu Number of jets after cut  *
 *                                                                            *
 *    Int_t _ElectronNumberMin;          -- Minimum Electron Number           *
 *    Int_t _ElectronNumberMax;          -- Maximum Electron Number           *
 *    Int_t _ElectronLeadingPt;          -- Leading Electron Pt Cut           *
 *    Int_t _ElectronSubLeadingPt;       -- SubLeading Electron Pt Cut        *
 *                                                                            *
 * History                                                                    *
 *      15 May 2011 - Created by Huaqiao ZHANG @ CERN                         *
 *****************************************************************************/

#ifndef CutElectronTighterPt_h
#define CutElectronTighterPt_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutElectronTighterPt : public HistoCut 
{

public:

  // Parameterized Constructor
  CutElectronTighterPt(EventContainer *EventContainerObj, TString electronType);
  
  // Destructor
  ~CutElectronTighterPt();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutElectronTighterPt"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hElectronNumberBefore;   // Histogram Electron Number before the cut
  myTH1F* _hElectronNumberAfter;    // Histogram Electron Number of jets after the cut

  // Cut parameters
  Int_t _ElectronNumberMin;         // Minimum Electron Number
  Int_t _ElectronNumberMax;         // Maximum Electron Number
  Int_t _ElectronLeadingPt;         // Leading Electron Pt
  Int_t _ElectronSubLeadingPt;      // SubLeading Electron Pt

  // Electron Type (veto, tight)
  TString electronType;

};


#endif







