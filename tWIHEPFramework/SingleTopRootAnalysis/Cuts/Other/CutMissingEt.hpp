/******************************************************************************
 * CutMissingEt.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMissingEt class                              *
 *    CutMissingEt()                    -- Parameterized Constructor          *
 *    ~CutMissingEt()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutMissingEt"             *
 *                                                                            *
 * Private Data Members of CutMissingEt                                       *
 *    myTH1F* _hMissingEtBefore         -- Hist of MissingEt before cuts      *
 *    myTH1F* _hMissingEtAfter          -- Hist of MissingEt after cuts       *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutMissingEt_h
#define CutMissingEt_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <vector>
#include <sstream>

class CutMissingEt : public HistoCut 
{

public:

  // Parameterized Constructor
  CutMissingEt(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutMissingEt();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutMissingEt"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hMissingEtBefore;    // Histogram Missing Et
  myTH1F* _hMissingEtAfter;     // Histogram Missing Et

  // Cut parameters
  Double_t _missingEtMin;           // Minimum Missing Et
  Double_t _missingEtMax;           // Maximum Missing Et

};


#endif







