/******************************************************************************
 * CutLeptonCharge.hpp                                                             *
 *                                                                            *
 * Cuts on lepton Number (can cut on unIsolated, tight, or veto Leptons)      *
 * Must pass unisolated, tight, or veto to constructor                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutLeptonCharge class                                *
 *    CutLeptonCharge()                     -- Parameterized Constructor           *
 *    ~CutLeptonCharge()                    -- Destructor                          *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutLeptonCharge"                  *
 *                                                                            *
 * Private Data Members of CutLeptonCharge class                                   *
 *    myTH1F* _hLeptonMassumberBefore;    -- Hist lepton Number before cut       *
 *    myTH1F* _hLeptonMassumberAfter;     -- Hist lepton Number of jets after cut*
 *                                                                            *
 *    Int_t _LeptonNumberMin;          -- Minimum Lepton Number               *
 *    Int_t _LeptonNumberMax;          -- Maximum Lepton Number               *
 *                                                                            *
 * History                                                                    *
 *      24 Mar 2009 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutLeptonCharge_h
#define CutLeptonCharge_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutLeptonCharge : public HistoCut 
{

public:

  // Parameterized Constructor
  CutLeptonCharge(EventContainer *EventContainerObj, TString leptonType);
  
  // Destructor
  ~CutLeptonCharge();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutLeptonCharge"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hLeptonChargeBefore;   // Histogram Lepton Number before the cut
  myTH1F* _hLeptonChargeAfter;    // Histogram Lepton Number of jets after the cut

  //Do we want same or opposite sign?
  bool _LeptonSameSign;	  // true to accept only within the region, false to reject region

  // Lepton Type (veto, tight, unIsolated)
  TString leptonType;

};


#endif







