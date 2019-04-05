/******************************************************************************
 * DumpEventNumber.hpp                                                         *
 *                                                                            *
 * Apply a default set of cuts to all classes deriving from BaseCuts          *
 *                                                                            *
 * Derived from BaseCut                                                       *
 *                                                                            *
 * Public Member Functions of DumpEventNumber class                            *
 *    DumpEventNumber()                        -- Parameterized Constructor    *
 *    ~DumpEventNumber()                       -- Destructor                   *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "DumpEventNumber"           *
 *                                                                            *

 * History                                                                    *
 *      16 June 2011 - Created by J. Holzbauer                                     *
 ******************************************************************************/

#ifndef DumpEventNumber_h
#define DumpEventNumber_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>
#include <fstream>
class DumpEventNumber : public HistoCut
{

public:
  
  // Parameterized constructor
  DumpEventNumber(EventContainer  *EventContainerObj, TString fileNameToPrint);
  
  // Destructor
  ~DumpEventNumber();

  // Book histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Return name of Cut
  inline std::string GetCutName() { return "DumpEventNumber"; };
  

private:

  TString _particle;
  TString _particleType;

  ofstream _fileToPrint;
};


#endif
