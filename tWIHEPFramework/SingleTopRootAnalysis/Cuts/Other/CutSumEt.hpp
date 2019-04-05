/******************************************************************************
 * CutSumEt.hpp                                                           *
 *                                                                            *
 * Cuts on Sum Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutSumEt class                              *
 *    CutSumEt()                    -- Parameterized Constructor          *
 *    ~CutSumEt()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutSumEt"             *
 *                                                                            *
 * Private Data Members of CutSumEt                                       *
 *    myTH1F* _hSumEtBefore         -- Hist of SumEt before cuts      *
 *    myTH1F* _hSumEtAfter          -- Hist of SumEt after cuts       *
 *    Int_t _sumEtMin;              -- Minimum Sum Et                 *
 *    Int_t _sumEtMin;              -- Minimum Sum Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutSumEt_h
#define CutSumEt_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutSumEt : public HistoCut 
{

public:

  // Parameterized Constructor
  CutSumEt(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutSumEt();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutSumEt"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hSumEtBefore;    // Histogram Sum Et
  myTH1F* _hSumEtAfter;     // Histogram Sum Et

  // Cut parameters
  Double_t _sumEtMin;           // Minimum Sum Et
  Double_t _sumEtMax;           // Maximum Sum Et

};


#endif







