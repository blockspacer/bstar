/******************************************************************************
 * CutTaggedJetN.hpp                                                          *
 *                                                                            *
 * Cuts on number of tagged jets in an event                                  *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTaggedJetN class                             *
 *    CutTaggedJetN()                   -- Parameterized Constructor          *
 *    ~CutTaggedJetN()                  -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutTaggedJetN"            *
 *                                                                            *
 * Private Data Members of CutTaggedJetN                                      *
 *    myTH1F* _hJetNumberBefore         -- Hist of Num of jets before cuts    *
 *    myTH1F* _hJetNumberAfter          -- Hist of Num of jets before cuts    *
 *    Int_t _JetNumberMin;              -- Minimum number of jets to require  *
 *    Int_t _JetNumberMax;              -- Maximum number of jets allowed     *
 *                                                                            *
 * History                                                                    *
 *      23 Feb 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutTaggedJetN_h
#define CutTaggedJetN_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutTaggedJetN : public HistoCut 
{

public:

  // Parameterized Constructor
  CutTaggedJetN(EventContainer *EventContainerObj, Int_t nbJetsDefault);
  
  // Destructor
  ~CutTaggedJetN();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutTaggedJetN"; };

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







