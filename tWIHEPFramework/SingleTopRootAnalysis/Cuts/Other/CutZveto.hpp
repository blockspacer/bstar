/******************************************************************************
 * CutZveto.hpp                                                               *
 *                                                                            *
 * Cuts on Sum Et                                                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutZveto class                                  *
 *    CutZveto()                    -- Parameterized Constructor              *
 *    ~CutZveto()                   -- Destructor                             *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutZveto"                 *
 *                                                                            *
 * Private Data Members of CutZveto                                           *
 *    myTH1F* _hZvetoBefore         -- Hist of Zveto before cuts              *
 *    myTH1F* _hZvetoAfter          -- Hist of Zveto after cuts               *
 *    Double_t _ZvetoMin;              -- Minimum Zveto window                *
 *    Double_t _ZvetoMin;              -- Minimum Zveto window                *
 *                                                                            *
 * History                                                                    *
 *      03 Dec 2010 - Created by H. Zhang                                     *
 *****************************************************************************/

#ifndef CutZveto_h
#define CutZveto_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutZveto : public HistoCut 
{

public:

  // Parameterized Constructor
  CutZveto(EventContainer *EventContainerObj,TString leptonTypePassed);
  
  // Destructor
  ~CutZveto();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutZveto"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hZvetoBefore;    // Histogram Sum Et
  myTH1F* _hZvetoAfter;     // Histogram Sum Et

  // Cut parameters
  Double_t _ZvetoMin;           // Minimum Sum Et
  Double_t _ZvetoMax;           // Maximum Sum Et

  TString leptonType; // Type of electron being plotted
};


#endif







