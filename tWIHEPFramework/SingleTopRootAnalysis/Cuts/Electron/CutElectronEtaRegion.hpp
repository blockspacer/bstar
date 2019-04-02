/******************************************************************************
 * CutElectronEtaRegion.hpp                                                   *
 *                                                                            *
 * Cuts event if electron lies in a certain eta region                        *
 * Cuts the range in both positive and negativ eta                            *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutElectronEtaRegion class                      *
 *    CutElectronEtaRegion()            -- Parameterized Constructor          *
 *    ~CutElectronEtaRegion()           -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutElectronEtaRegion"     *
 *                                                                            *
 * Private Data Members of CutElectronEtaRegion class                         *
 *    myTH1F* _hElectronEtaRegionBefore; -- Hist e Eta before cut             *
 *    myTH1F* _hElectronEtaRegionAfter;  -- Hist e Eta of jets after cut      *
 *                                                                            *
 *    Double_t _ElectronEtaRegionMin;    -- Minimum Electron Eta              *
 *    Double_t _ElectronEtaRegionMax;    -- Maximum Electron Eta              *
 *                                                                            *
 * History                                                                    *
 *      21 Aug 2007 - Created by P. Ryan                                      *
 *      18 Jan 2008 - Use of tight/veto electrons requires that eta region    *
 *                    cut uses all electrons                                  *
 *****************************************************************************/

#ifndef CutElectronEtaRegion_h
#define CutElectronEtaRegion_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutElectronEtaRegion : public HistoCut 
{

public:

  // Parameterized Constructor
  CutElectronEtaRegion(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutElectronEtaRegion();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutElectronEtaRegion"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hElectronEtaRegionBefore;    // Histogram Electron 1 Eta before the cut
  myTH1F* _hElectronEtaRegionAfter;     // Histogram Electron 1 Eta of jets after the cut

  // Cut parameters
  Int_t _ElectronNumberMin;             // Minimum Electron Number
  Int_t _ElectronNumberMax;             // Maximum Electron Number

  Double_t _ElectronEtaRegionMin;       // Minimum Electron 1 Eta
  Double_t _ElectronEtaRegionMax;       // Maximum Electron 1 Eta

};


#endif







