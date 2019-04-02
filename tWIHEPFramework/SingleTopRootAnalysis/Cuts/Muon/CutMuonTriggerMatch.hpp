/******************************************************************************
 * CutMuonTriggerMatch.hpp                                                           *
 *                                                                            *
 * Cuts on muon N (can cut on tight or veto Muons)               *
 * Must pass tight or veto to constructor                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutMuonTriggerMatch class                              *
 *    CutMuonTriggerMatch()                     -- Parameterized Constructor         *
 *    ~CutMuonTriggerMatch()                    -- Destructor                        *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutMuonTriggerMatch"                *
 *                                                                            *
 * Private Data Members of CutMuonTriggerMatch class                                 *
 *    myTH1F* _hMuonTriggerMatchBefore;    -- Hist mu N before cut         *
 *    myTH1F* _hMuonTriggerMatchAfter;     -- Hist mu N of jets after cut  *
 *                                                                            *
 *    Int_t _MuonTriggerMatchMin;          -- Minimum Muon N           *
 *    Int_t _MuonTriggerMatchMax;          -- Maximum Muon N           *
 *                                                                            *
 * History                                                                    *
 *      15 Jan 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutMuonTriggerMatch_h
#define CutMuonTriggerMatch_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutMuonTriggerMatch : public HistoCut 
{

public:

  // Parameterized Constructor
  CutMuonTriggerMatch(EventContainer *EventContainerObj, TString muonType);
  
  // Destructor
  ~CutMuonTriggerMatch();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutMuonTriggerMatch"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hMuonTriggerMatchBefore;   // Histogram Muon N before the cut
  myTH1F* _hMuonTriggerMatchAfter;    // Histogram Muon N of jets after the cut

  // Cut parameters
  Bool_t _MuonTriggerMatch;        
 
  // Muon Type (veto, tight)
  TString muonType;

};


#endif







