/******************************************************************************
 * CutTauN.hpp                                                                *
 *                                                                            *
 * Cuts on electron Number                                                    *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTauN class                                   *
 *    CutTauN()                     -- Parameterized Constructor              *
 *    ~CutTauN()                    -- Destructor                             *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutTauN"                     *
 *                                                                            *
 * Private Data Members of CutTauN class                                      *
 *    myTH1F* _hTauNumberBefore;    -- Hist e Number before cut               *
 *    myTH1F* _hTauNumberAfter;     -- Hist e Number of jets after cut        *
 *                                                                            *
 *    Int_t _TauNumberMin;          -- Minimum Electron Number                *
 *    Int_t _TauNumberMax;          -- Maximum Electron Number                *
 *                                                                            *
 * History                                                                    *
 *      22 Aug 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutTauN_h
#define CutTauN_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutTauN : public HistoCut 
{

public:

  // Parameterized Constructor
  CutTauN(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutTauN();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutTauN"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hTauNumberBefore;   // Histogram Tau Number before the cut
  myTH1F* _hTauNumberAfter;    // Histogram Tau Number of jets after the cut

  // Cut parameters
  Int_t _TauNumberMin;         // Minimum Tau Number
  Int_t _TauNumberMax;         // Maximum Tau Number

};


#endif







