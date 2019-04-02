/******************************************************************************
 * CutJetEta1.hpp                                                             *
 *                                                                            *
 * Cuts on Eta of 1st Jet                                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetEta1 class                                *
 *    CutJetEta1()                      -- Parameterized Constructor          *
 *    ~CutJetEta1()                     -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetEta1"               *
 *                                                                            *
 * Private Data Members of CutJetEta1                                         *
 *    myTH1F* _hJetEta1Before            -- Hist of Num of jets before cuts   *
 *    myTH1F* _hJetEta1After             -- Hist of Num of jets before cuts   *
 *    Int_t _JetEta1Min;                 -- Minimum Jet 1 Eta                 *
 *    Int_t _JetEta1Min;                 -- Minimum Jet 2 Eta                 *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by P. Ryan                                      *
 *      11 Jan 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#ifndef CutJetEta1_h
#define CutJetEta1_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutJetEta1 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutJetEta1(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutJetEta1();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutJetEta1"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetEta1Before;    // Histogram Eta of jet 1 before the cut
  myTH1F* _hJetEta1After;     // Histogram Eta of jet 1 after the cut

  // Cut parameters
  Double_t _JetEta1Min;           // Minimum Eta of jet 1
  Double_t _JetEta1Max;           // Maximum Eta of jet 1

};


#endif







