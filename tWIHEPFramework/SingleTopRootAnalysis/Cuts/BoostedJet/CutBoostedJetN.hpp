/******************************************************************************
 * CutBoostedJetN.hpp                                                                *
 *                                                                            *
 * Cuts on number of jets in an event                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutBoostedJetN class                                   *
 *    CutBoostedJetN()                         -- Parameterized Constructor          *
 *    ~CutBoostedJetN()                        -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutBoostedJetN"                  *
 *                                                                            *
 * Private Data Members of CutBoostedJetN                                            *
 *    myTH1F* _hBoostedJetNumberBefore         -- Hist of Num of jets before cuts    *
 *    myTH1F* _hBoostedJetNumberAfter          -- Hist of Num of jets before cuts    *
 *    Int_t _BoostedJetNumberMin;              -- Minimum number of jets to require  *
 *    Int_t _BoostedJetNumberMax;              -- Maximum number of jets allowed     *
 *                                                                            *
 * History                                                                    *
 *       7 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      10 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      11 Jan 2007 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#ifndef CutBoostedJetN_h
#define CutBoostedJetN_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutBoostedJetN : public HistoCut 
{

public:

  // Parameterized Constructor
  CutBoostedJetN(EventContainer *EventContainerObj, Int_t nBoostedJetsDefault = -1);
  
  // Destructor
  ~CutBoostedJetN();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutBoostedJetN"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hBoostedJetNumberBefore;    // Histogram number of jets before the cut
  myTH1F* _hBoostedJetNumberAfter;     // Histogram number of jets after the cut

  // Cut parameters
  Int_t _BoostedJetNumberMin;           // Minimum number of jets to require
  Int_t _BoostedJetNumberMax;           // Maximum number of jets allowed

};


#endif







