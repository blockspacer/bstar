/******************************************************************************
 * CutMuonTighterPt.hpp                                                        *
 *                                                                            *
 * Cuts on muon Number (can cut on All, UnIsolated, Tight, or Veto Muons)     *
 * Must pass All, UnIsolated, Tight, or Veto to constructor                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutMuonTighterPt class                           *
 *    CutMuonTighterPt()                     -- Parameterized Constructor      *
 *    ~CutMuonTighterPt()                    -- Destructor                     *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutMuonTighterPt"             *
 *                                                                            *
 * Private Data Members of CutMuonTighterPt class                              *
 *    myTH1F* _hMuonNumberBefore;    -- Hist mu Number before cut             *
 *    myTH1F* _hMuonNumberAfter;     -- Hist mu Number of jets after cut      *
 *                                                                            *
 *    Int_t _MuonNumberMin;          -- Minimum Muon Number                   *
 *    Int_t _MuonNumberMax;          -- Maximum Muon Number                   *
 *    Int_t _MuonLeadingPt;          -- Leading Muon Pt Cut                   *
 *    Int_t _MuonSubLeadingPt;       -- SubLeading Muon Pt Cut                *
 *                                                                            *
 * History                                                                    *
 *      15 May 2011 - Created by Huaqiao ZHANG @ CERN                         *
 *****************************************************************************/

#ifndef CutMuonTighterPt_h
#define CutMuonTighterPt_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutMuonTighterPt : public HistoCut 
{

public:

  // Parameterized Constructor
  CutMuonTighterPt(EventContainer *EventContainerObj, TString muonType);
  
  // Destructor
  ~CutMuonTighterPt();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutMuonTighterPt"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hMuonNumberBefore;   // Histogram Muon Number before the cut
  myTH1F* _hMuonNumberAfter;    // Histogram Muon Number of jets after the cut

  // Cut parameters
  Int_t _MuonNumberMin;         // Minimum Muon Number
  Int_t _MuonNumberMax;         // Maximum Muon Number
  Int_t _MuonLeadingPt;         // Leading Muon Pt
  Int_t _MuonSubLeadingPt;      // SubLeading Muon Pt

  // Muon Type (veto, tight, unIsolated)
  TString muonType;

};


#endif







