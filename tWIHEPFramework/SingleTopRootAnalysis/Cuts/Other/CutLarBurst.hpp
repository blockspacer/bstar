/******************************************************************************
 * CutLarBurst.hpp                                                         *
 *                                                                            *
 * Apply a default set of cuts to all classes deriving from BaseCuts          *
 *                                                                            *
 * Derived from BaseCut                                                       *
 *                                                                            *
 * Public Member Functions of CutLarBurst class                            *
 *    CutLarBurst()                        -- Parameterized Constructor    *
 *    ~CutLarBurst()                       -- Destructor                   *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutLarBurst"           *
 *                                                                            *

 * History                                                                    *
 *      16 June 2011 - Created by J. Holzbauer                                     *
 ******************************************************************************/

#ifndef CutLarBurst_h
#define CutLarBurst_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>
class CutLarBurst : public HistoCut
{

public:
  
  // Parameterized constructor
  CutLarBurst(EventContainer  *EventContainerObj);
  
  // Destructor
  ~CutLarBurst();

  // Book histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Return name of Cut
  inline std::string GetCutName() { return "CutLarBurst"; };
  

private:

  // Histograms
  myTH1F* _hLarBurstBefore;     // Histogram Before GRL Cut
  myTH1F* _hLarBurstAfter;      // Histogram After GRL Cut
  
  // Tmp
  Int_t tmpCounter;
  
  // Cut parameters
  TString _GoodRunsListName; 
  
  
  Int_t LarBurstNumber;

};


#endif







