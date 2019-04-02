/******************************************************************************
 * CutBadPileupEvent.hpp                                                              *
 *                                                                            *
 * Cuts on charge and flavor of lepton                                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutBadPileupEvent class                                 *
 *    CutBadPileupEvent()                       -- Parameterized Constructor          *
 *    ~CutBadPileupEvent()                      -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutBadPileupEvent"                *
 *                                                                            *
 * Private Data Members of CutJetN                                           *
 *    Int_t _leptonCharge                -- lepton charge cut value           *
 *    Int_t _leptonFlavor                -- lepton flavor cut value           *
 *    myTH1F* _hLeptonCharge             -- Hist of lepton charge             *
 *    myTH1F* _hLeptonFlavor             -- Hist of lepton flavor             *
 *                                                                            *
 * History                                                                    *
 *      09 Jun 2010 - Created by J. Holzbauer                                *
 *****************************************************************************/

#ifndef CutBadPileupEvent_h
#define CutBadPileupEvent_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

class CutBadPileupEvent : public HistoCut
{

public:

  // Parameterized Constructor
  CutBadPileupEvent(EventContainer *obj);
  
  // Destructor
  ~CutBadPileupEvent();

  // Book Histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Get the name describing the cut
  inline std::string GetCutName() { return "CutBadPileupEvent"; }


private:

  int _ENint;
  // Histograms
  myTH1F*   _hBadPileupEventBefore;  // Lepton charge before cut
  myTH1F*   _hBadPileupEventAfter;   // Lepton charge after cut


};


#endif







