/******************************************************************************
 * CutJetTriggerMatch.hpp                                                           *
 *                                                                            *
 * Cuts on jet N (can cut on tight or veto Jets)               *
 * Must pass tight or veto to constructor                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutJetTriggerMatch class                              *
 *    CutJetTriggerMatch()                     -- Parameterized Constructor         *
 *    ~CutJetTriggerMatch()                    -- Destructor                        *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutJetTriggerMatch"                *
 *                                                                            *
 * Private Data Members of CutJetTriggerMatch class                                 *
 *    myTH1F* _hJetTriggerMatchBefore;    -- Hist mu N before cut         *
 *    myTH1F* _hJetTriggerMatchAfter;     -- Hist mu N of jets after cut  *
 *                                                                            *
 *    Int_t _JetTriggerMatchMin;          -- Minimum Jet N           *
 *    Int_t _JetTriggerMatchMax;          -- Maximum Jet N           *
 *                                                                            *
 * History                                                                    *
 *      15 Jan 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutJetTriggerMatch_h
#define CutJetTriggerMatch_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutJetTriggerMatch : public HistoCut 
{

public:

  // Parameterized Constructor
  CutJetTriggerMatch(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutJetTriggerMatch();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutJetTriggerMatch"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetTriggerMatchBefore;   // Histogram Jet N before the cut
  myTH1F* _hJetTriggerMatchAfter;    // Histogram Jet N of jets after the cut

  // Cut parameters
  Bool_t _JetTriggerMatch;        
 

};


#endif







