/*****************************************************************************
 * CutTaggedJetEta1.hpp                                                             *
 *                                                                            *
 * Cuts on Eta of 1st TaggedJet                                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTaggedJetEta1 class                                *
 *    CutTaggedJetEta1()                      -- Parameterized Constructor          *
 *    ~CutTaggedJetEta1()                     -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutTaggedJetEta1"               *
 *                                                                            *
 * Private Data Members of CutTaggedJetEta1                                         *
 *    myTH1F* _hTaggedJetEta1Before            -- Hist of Num of jets before cuts   *
 *    myTH1F* _hTaggedJetEta1After             -- Hist of Num of jets before cuts   *
 *    Int_t _TaggedJetEta1Min;                 -- Minimum TaggedJet 1 Eta                 *
 *    Int_t _TaggedJetEta1Min;                 -- Minimum TaggedJet 2 Eta                 *
 *                                                                            *
 * History                                                                    *
 *      28 Aug 2007- created by J. Holzbauer based on CutJetEta1 by P. Ryan        *
 *****************************************************************************/

#ifndef CutTaggedJetEta1_h
#define CutTaggedJetEta1_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutTaggedJetEta1 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutTaggedJetEta1(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutTaggedJetEta1();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutTaggedJetEta1"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hTaggedJetEta1Before;    // Histogram Eta of jet 1 before the cut
  myTH1F* _hTaggedJetEta1After;     // Histogram Eta of jet 1 after the cut

  // Cut parameters
  Double_t _TaggedJetEta1Min;           // Minimum Eta of jet 1
  Double_t _TaggedJetEta1Max;           // Maximum Eta of jet 1

};


#endif






