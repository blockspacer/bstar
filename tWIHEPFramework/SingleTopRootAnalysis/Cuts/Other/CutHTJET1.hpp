/******************************************************************************
 * CutHTJET1.hpp                                                           *
 *                                                                            *
 * Cuts on Sum Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutHTJET1 class                              *
 *    CutHTJET1()                    -- Parameterized Constructor          *
 *    ~CutHTJET1()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutHTJET1"             *
 *                                                                            *
 * Private Data Members of CutHTJET1                                       *
 *    myTH1F* _hHTJET1Before         -- Hist of HTJET1 before cuts      *
 *    myTH1F* _hHTJET1After          -- Hist of HTJET1 after cuts       *
 *    Int_t _htJET1Min;              -- Minimum Sum Et                 *
 *    Int_t _htJET1Min;              -- Minimum Sum Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutHTJET1_h
#define CutHTJET1_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutHTJET1 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutHTJET1(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutHTJET1();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutHTJET1"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hHTJET1Before;    // Histogram Sum Et
  myTH1F* _hHTJET1After;     // Histogram Sum Et

  // Cut parameters
  Double_t _htJET1Min;           // Minimum Sum Et
  Double_t _htJET1Max;           // Maximum Sum Et

};


#endif







