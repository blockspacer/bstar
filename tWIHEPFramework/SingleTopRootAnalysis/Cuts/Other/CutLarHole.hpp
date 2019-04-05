/******************************************************************************
 * CutLarHole.hpp                                                         *
 *                                                                            *
 * Apply a default set of cuts to all classes deriving from BaseCuts          *
 *                                                                            *
 * Derived from BaseCut                                                       *
 *                                                                            *
 * Public Member Functions of CutLarHole class                            *
 *    CutLarHole()                        -- Parameterized Constructor    *
 *    ~CutLarHole()                       -- Destructor                   *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutLarHole"           *
 *                                                                            *

 * History                                                                    *
 *      16 June 2011 - Created by J. Holzbauer                                     *
 ******************************************************************************/

#ifndef CutLarHole_h
#define CutLarHole_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>
class CutLarHole : public HistoCut
{

public:
  
  // Parameterized constructor
  CutLarHole(EventContainer  *EventContainerObj);
  
  // Destructor
  ~CutLarHole();

  // Book histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Return name of Cut
  inline std::string GetCutName() { return "CutLarHole"; };
  

private:

  // Histograms
  myTH1F* _hLarHoleBefore;     // Histogram Before GRL Cut
  myTH1F* _hLarHoleAfter;      // Histogram After GRL Cut
  
  // Tmp
  Int_t tmpCounter;
  
  // Cut parameters
  TString _GoodRunsListName; 
  
  
  Int_t LarHoleNumber;

};


#endif







