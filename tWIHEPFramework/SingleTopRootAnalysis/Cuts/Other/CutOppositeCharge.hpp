/******************************************************************************
 * CutOppositeCharge.hpp                                                      *
 *                                                                            *
 * Selects events with muon and electron having opposite charge               *
 * Compares 1st veto muon 1st veto electron                                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutOppositeCharge class                         *
 *    CutOppositeCharge()                       -- Parameterized Constructor  *
 *    ~CutOppositeCharge()                      -- Destructor                 *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutOppositeCharge"        *
 *                                                                            *
 * Private Data Members of CutJetN                                            *
 *    Int_t _muonCharge                -- muon charge                         *
 *    Int_t _electronCharge            -- electron charge                     *
 *    myTH1F*   _hOppositeChargeBefore -- Opposite charge before cut          *
 *    myTH1F*   _hOppositeChargeAfter  -- Opposite charge after cut           *
 *                                                                            *
 * History                                                                    *
 *      17 Feb 2008 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutOppositeCharge_h
#define CutOppositeCharge_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutOppositeCharge : public HistoCut
{

public:

  // Parameterized Constructor
  CutOppositeCharge(EventContainer *obj);
  
  // Destructor
  ~CutOppositeCharge();

  // Book Histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Get the name describing the cut
  inline std::string GetCutName() { return "CutOppositeCharge"; }


private:

  // Cut parameters
  Int_t _muonCharge;         // Charge of 1st veto muon
  Int_t _electronCharge;     // Charge of 1st veto electron

  // Histograms
  myTH1F*   _hOppositeChargeBefore;  // Opposite charge before cut
  myTH1F*   _hOppositeChargeAfter;   // Opposite charge after cut


};


#endif







