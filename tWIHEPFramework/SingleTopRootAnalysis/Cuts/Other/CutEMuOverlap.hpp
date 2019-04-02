/******************************************************************************
 * CutEMuOverlap.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutEMuOverlap class                              *
 *    CutEMuOverlap()                    -- Parameterized Constructor          *
 *    ~CutEMuOverlap()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutEMuOverlap"             *
 *                                                                            *
 * Private Data Members of CutEMuOverlap                                       *
 *    myTH1F* _hEMuOverlapBefore         -- Hist of EMuOverlap before cuts      *
 *    myTH1F* _hEMuOverlapAfter          -- Hist of EMuOverlap after cuts       *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutEMuOverlap_h
#define CutEMuOverlap_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutEMuOverlap : public HistoCut 
{

public:

  // Parameterized Constructor
  CutEMuOverlap(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutEMuOverlap();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutEMuOverlap"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hEMuOverlapBefore;    // Histogram Missing Et
  myTH1F* _hEMuOverlapAfter;     // Histogram Missing Et

};


#endif







