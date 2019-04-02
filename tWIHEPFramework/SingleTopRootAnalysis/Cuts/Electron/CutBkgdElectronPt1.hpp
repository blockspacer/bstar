/******************************************************************************
 * CutBkgdElectronPt1.hpp                                                         *
 *                                                                            *
 * Cuts on electron Pt1 (can cut on tight or veto Electrons)                  *
 * Must pass tight or veto to constructor                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutBkgdElectronPt1 class                            *
 *    CutBkgdElectronPt1()                     -- Parameterized Constructor       *
 *    ~CutBkgdElectronPt1()                    -- Destructor                      *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutBkgdElectronPt1"              *
 *                                                                            *
 * Private Data Members of CutBkgdElectronPt1 class                               *
 *    myTH1F* _hBkgdElectronPt1Before;    -- Hist mu Pt1 before cut               *
 *    myTH1F* _hBkgdElectronPt1After;     -- Hist mu Pt1 of jets after cut        *
 *                                                                            *
 *    Int_t _BkgdElectronPt1Min;          -- Minimum Electron Pt1                 *
 *    Int_t _BkgdElectronPt1Max;          -- Maximum Electron Pt1                 *
 *                                                                            *
 * History                                                                    *
 *      06 May 2011 -Created by J. Holzbauer                                  *
 *****************************************************************************/

#ifndef CutBkgdElectronPt1_h
#define CutBkgdElectronPt1_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutBkgdElectronPt1 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutBkgdElectronPt1(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutBkgdElectronPt1();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutBkgdElectronPt1"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hBkgdElectronPt1Before;   // Histogram Electron Pt1 before the cut
  myTH1F* _hBkgdElectronPt1After;    // Histogram Electron Pt1 of jets after the cut

  // Cut parameters
  Double_t _BkgdElectronPt1Min;         // Minimum Electron Pt1
  Double_t _BkgdElectronPt1Max;         // Maximum Electron Pt1
  Double_t BkgdElectronPt1;

  
};


#endif







