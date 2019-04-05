/******************************************************************************
 * CutMCJetN.hpp                                                              *
 *                                                                            *
 * Cuts on number of MC jets in an event                                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMCJetN class                                 *
 *    CutMCJetN()                       -- Parameterized Constructor          *
 *    ~CutMCJetN()                      -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutMCJetN"                *
 *                                                                            *
 * Private Data Members of CutMCJetN                                          *
 *    myTH1F* _hJetNumberBefore         -- Hist of Num of jets before cuts    *
 *    myTH1F* _hJetNumberAfter          -- Hist of Num of jets before cuts    *
 *    Int_t _JetNumberMin;              -- Minimum number of jets to require  *
 *    Int_t _JetNumberMax;              -- Maximum number of jets allowed     *
 *                                                                            *
 * History                                                                    *
 *       26 Mar 2008 - Created by P. Ryan                                     *
 *****************************************************************************/

#ifndef CutMCJetN_h
#define CutMCJetN_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutMCJetN : public HistoCut 
{

public:

  // Parameterized Constructor
  CutMCJetN(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutMCJetN();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutMCJetN"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetNumberBefore;    // Histogram number of jets before the cut
  myTH1F* _hJetNumberAfter;     // Histogram number of jets after the cut

  // Cut parameters
  Int_t _JetNumberMin;           // Minimum number of jets to require
  Int_t _JetNumberMax;           // Maximum number of jets allowed

};


#endif







