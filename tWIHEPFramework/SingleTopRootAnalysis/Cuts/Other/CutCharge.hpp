/******************************************************************************
 * CutCharge.hpp                                                              *
 *                                                                            *
 * Cuts on charge and flavor of lepton                                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutCharge class                                 *
 *    CutCharge()                       -- Parameterized Constructor          *
 *    ~CutCharge()                      -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutCharge"                *
 *                                                                            *
 * Private Data Members of CutJetN                                           *
 *    Int_t _leptonCharge                -- lepton charge cut value           *
 *    Int_t _leptonFlavor                -- lepton flavor cut value           *
 *    myTH1F* _hLeptonCharge             -- Hist of lepton charge             *
 *    myTH1F* _hLeptonFlavor             -- Hist of lepton flavor             *
 *                                                                            *
 * History                                                                    *
 *      28 Apr 2004 - Created by R. Schwienhorst                              *
 *      14 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#ifndef CutCharge_h
#define CutCharge_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutCharge : public HistoCut
{

public:

  // Parameterized Constructor
  CutCharge(EventContainer *obj);
  
  // Destructor
  ~CutCharge();

  // Book Histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Get the name describing the cut
  inline std::string GetCutName() { return "CutCharge"; }


private:

  // Cut parameters
  Int_t _leptonCharge;  // Sign of the lepton charge
  Int_t _leptonFlavor;  // Flavor of the lepton

  // Histograms
  myTH1F*   _hLeptonChargeBefore;  // Lepton charge before cut
  myTH1F*   _hLeptonFlavorBefore;  // Lepton flavor before cut
  myTH1F*   _hLeptonChargeAfter;   // Lepton charge after cut
  myTH1F*   _hLeptonFlavorAfter;   // Lepton flavor after cut

};


#endif







