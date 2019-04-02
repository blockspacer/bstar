/******************************************************************************
 * DumpEventInfo.hpp                                                         *
 *                                                                            *
 * Apply a default set of cuts to all classes deriving from BaseCuts          *
 *                                                                            *
 * Derived from BaseCut                                                       *
 *                                                                            *
 * Public Member Functions of DumpEventInfo class                            *
 *    DumpEventInfo()                        -- Parameterized Constructor    *
 *    ~DumpEventInfo()                       -- Destructor                   *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "DumpEventInfo"           *
 *                                                                            *

 * History                                                                    *
 *      16 June 2011 - Created by J. Holzbauer                                     *
 ******************************************************************************/

#ifndef DumpEventInfo_h
#define DumpEventInfo_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>
class DumpEventInfo : public HistoCut
{

public:
  
  // Parameterized constructor
  DumpEventInfo(EventContainer  *EventContainerObj);
  
  // Destructor
  ~DumpEventInfo();

  // Book histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Return name of Cut
  inline std::string GetCutName() { return "DumpEventInfo"; };
  


};


#endif
