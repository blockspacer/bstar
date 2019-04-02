/******************************************************************************
 * CutMuonPhiRegion.hpp                                                       *
 *                                                                            *
 * Cuts event if muon lies in a certain phi region                            *
 * Cuts the range in both positive and negative phi                           *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMuonPhiRegion class                          *
 *    CutMuonPhiRegion()                -- Parameterized Constructor          *
 *    ~CutMuonPhiRegion()               -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutMuonPhiRegion"         *
 *                                                                            *
 * Private Data Members of CutMuonPhiRegion class                             *
 *    myTH1F* _hMuonPhiRegionBefore; -- Hist e Phi before cut                 *
 *    myTH1F* _hMuonPhiRegionAfter;  -- Hist e Phi of jets after cut          *
 *                                                                            *
 *    Double_t _MuonPhiRegionMin;    -- Minimum Muon Phi                      *
 *    Double_t _MuonPhiRegionMax;    -- Maximum Muon Phi                      *
 *                                                                            *
 * History                                                                    *
 *      19 Nov 2008 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutMuonPhiRegion_h
#define CutMuonPhiRegion_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutMuonPhiRegion : public HistoCut 
{

public:

  // Parameterized Constructor
  CutMuonPhiRegion(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutMuonPhiRegion();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutMuonPhiRegion"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hMuonPhiRegionBefore;    // Histogram Muon 1 Phi before the cut
  myTH1F* _hMuonPhiRegionAfter;     // Histogram Muon 1 Phi of jets after the cut

  // Cut parameters
  Int_t _MuonNumberMin;             // Minimum Muon Number
  Int_t _MuonNumberMax;             // Maximum Muon Number

  Double_t _MuonPhiRegionMin;       // Minimum Muon 1 Phi
  Double_t _MuonPhiRegionMax;       // Maximum Muon 1 Phi

};


#endif







