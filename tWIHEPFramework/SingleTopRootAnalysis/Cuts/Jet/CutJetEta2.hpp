/******************************************************************************
 * CutJetEta2.hpp                                                             *
 *                                                                            *
 * Cuts on Eta of 2st Jet                                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetEta2 class                                *
 *    CutJetEta2()                      -- Parameterized Constructor          *
 *    ~CutJetEta2()                     -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetEta2"               *
 *                                                                            *
 * Private Data Members of CutJetEta2                                         *
 *    myTH1F* _hJetEta2Before            -- Hist of Num of jets before cuts   *
 *    myTH1F* _hJetEta2After             -- Hist of Num of jets before cuts   *
 *    Int_t _JetEta2Min;                 -- Minimum Jet 2 Eta                 *
 *    Int_t _JetEta2Min;                 -- Minimum Jet 2 Eta                 *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by P. Ryan                                      *
 *      11 Jan 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#ifndef CutJetEta2_h
#define CutJetEta2_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutJetEta2 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutJetEta2(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutJetEta2();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutJetEta2"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetEta2Before;    // Histogram Eta of jet 2 before the cut
  myTH1F* _hJetEta2After;     // Histogram Eta of jet 2 after the cut

  // Cut parameters
  Double_t _JetEta2Min;           // Minimum Eta of jet 2
  Double_t _JetEta2Max;           // Maximum Eta of jet 2

};


#endif







