/******************************************************************************
 * CutLeptonOppositeCharge.hpp                                                *
 *                                                                            *
 * Selects events with 2 leptons having opposite charge                       *
 * Requires that the event has EXACTLY two leptons                            *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutLeptonOppositeCharge class                   *
 *    CutLeptonOppositeCharge()   -- Parameterized Constructor                *
 *    ~CutLeptonOppositeCharge()  -- Destructor                               *
 *    BookHistogram()             -- Book histograms                          *
 *    Apply()                     -- Apply cuts and fill histograms           *
 *    GetCutName()                -- Returns "CutLeptonOppositeCharge"        *
 *                                                                            *
 * Private Data Members of CutJetN                                            *
 *    Int_t _muonCharge                -- muon charge                         *
 *    Int_t _electronCharge            -- electron charge                     *
 *    myTH1F*   _hOppositeChargeBefore -- Opposite charge before cut          *
 *    myTH1F*   _hOppositeChargeAfter  -- Opposite charge after cut           *
 *                                                                            *
 * History                                                                    *
 *      24 Mar 2009 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutLeptonOppositeCharge_h
#define CutLeptonOppositeCharge_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutLeptonOppositeCharge : public HistoCut
{

public:

  // Parameterized Constructor
  CutLeptonOppositeCharge(EventContainer *obj);
  
  // Destructor
  ~CutLeptonOppositeCharge();

  // Book Histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Get the name describing the cut
  inline std::string GetCutName() { return "CutLeptonOppositeCharge"; }


private:

  // Cut parameters
  Int_t _muonCharge;         // Charge of 1st veto muon
  Int_t _electronCharge;     // Charge of 1st veto electron

  // Histograms
  myTH1F*   _hOppositeChargeBefore;  // Opposite charge before cut
  myTH1F*   _hOppositeChargeAfter;   // Opposite charge after cut


};


#endif







