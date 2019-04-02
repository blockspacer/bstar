/******************************************************************************
 * CutBadJet.hpp                                                         *
 *                                                                            *
 * Apply a default set of cuts to all classes deriving from BaseCuts          *
 *                                                                            *
 * Derived from BaseCut                                                       *
 *                                                                            *
 * Public Member Functions of CutBadJet class                            *
 *    CutBadJet()                        -- Parameterized Constructor    *
 *    ~CutBadJet()                       -- Destructor                   *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutBadJet"           *
 *                                                                            *

 * History                                                                    *
 *      16 June 2011 - Created by J. Holzbauer                                     *
 ******************************************************************************/

#ifndef CutBadJet_h
#define CutBadJet_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>
class CutBadJet : public HistoCut
{

public:
  
  // Parameterized constructor
  CutBadJet(EventContainer  *EventContainerObj);
  
  // Destructor
  ~CutBadJet();

  // Book histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Return name of Cut
  inline std::string GetCutName() { return "CutBadJet"; };
  

private:

  // Histograms
  myTH1F* _hBadJetBefore;     // Histogram Before GRL Cut
  myTH1F* _hBadJetAfter;      // Histogram After GRL Cut
  
  // Tmp
  Int_t tmpCounter;
    
  
  Int_t BadJetNumber;

};


#endif







