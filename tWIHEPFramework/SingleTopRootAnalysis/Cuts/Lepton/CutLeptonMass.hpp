/******************************************************************************
 * CutLeptonMass.hpp                                                             *
 *                                                                            *
 * Cuts on lepton Number (can cut on unIsolated, tight, or veto Leptons)      *
 * Must pass unisolated, tight, or veto to constructor                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutLeptonMass class                                *
 *    CutLeptonMass()                     -- Parameterized Constructor           *
 *    ~CutLeptonMass()                    -- Destructor                          *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutLeptonMass"                  *
 *                                                                            *
 * Private Data Members of CutLeptonMass class                                   *
 *    myTH1F* _hLeptonMassumberBefore;    -- Hist lepton Number before cut       *
 *    myTH1F* _hLeptonMassumberAfter;     -- Hist lepton Number of jets after cut*
 *                                                                            *
 *    Int_t _LeptonNumberMin;          -- Minimum Lepton Number               *
 *    Int_t _LeptonNumberMax;          -- Maximum Lepton Number               *
 *                                                                            *
 * History                                                                    *
 *      24 Mar 2009 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutLeptonMass_h
#define CutLeptonMass_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutLeptonMass : public HistoCut 
{

public:

  // Parameterized Constructor
  CutLeptonMass(EventContainer *EventContainerObj, TString leptonType, TString massCutPostfix = "");
  
  // Destructor
  ~CutLeptonMass();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutLeptonMass"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hLeptonMassBefore;   // Histogram Lepton Number before the cut
  myTH1F* _hLeptonMassAfter;    // Histogram Lepton Number of jets after the cut

  // Cut parameters
  Float_t _LeptonMassMin;         // Minimum Dilepton mass
  Float_t _LeptonMassMax;         // Maximum Dilepton mass
  
  bool _LeptonMassReject;	  // true to accept only within the region, false to reject region

  // Lepton Type (veto, tight, unIsolated)
  TString leptonType;
  TString massCutPostfix;

};


#endif







