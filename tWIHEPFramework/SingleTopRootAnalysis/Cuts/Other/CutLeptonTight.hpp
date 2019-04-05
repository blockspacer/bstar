/******************************************************************************
 * CutLeptonTight.cpp                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutLeptonTight class                            *
 *    CutLeptonTight()                     -- Parameterized Constructor       *
 *    ~CutLeptonTight()                    -- Destructor                      *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutLeptonTight"              *
 *                                                                            *
 * Private Data Members of CutLeptonTight class                               *
 *    myTH1F* _hLeptonTightBefore;    -- Hist before cut                      *
 *    myTH1F* _hLeptonTightAfter;     -- Hist after cut                       *
 *                                                                            *
 * This is meant for use with trees like JetsMedium which require EXACTLY one *
 * loose lepton.  This class removes events in that tree in which the lepton  *
 * does not also satisft the tight requirement, turning the tree into one that*
 * requires exactly one tight lepton instead.                                 *
 *                                                                            *
 * History                                                                    *
 *      Nov 10 2010:  Created by J. Holzbauer                                 *
 *****************************************************************************/

#ifndef CutLeptonTight_h
#define CutLeptonTight_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutLeptonTight : public HistoCut 
{

public:

  // Parameterized Constructor
  CutLeptonTight(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutLeptonTight();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutLeptonTight"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hLeptonTightBefore;    // Histogram Missing Et
  myTH1F* _hLeptonTightAfter;     // Histogram Missing Et

 };


#endif







