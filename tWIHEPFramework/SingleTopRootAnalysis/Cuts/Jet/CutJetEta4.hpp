/******************************************************************************
 * CutJetEta4.hpp                                                             *
 *                                                                            *
 * Cuts on Eta of 2st Jet                                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetEta4 class                                *
 *    CutJetEta4()                      -- Parameterized Constructor          *
 *    ~CutJetEta4()                     -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetEta4"               *
 *                                                                            *
 * Private Data Members of CutJetEta4                                         *
 *    myTH1F* _hJetEta4Before            -- Hist of Num of jets before cuts   *
 *    myTH1F* _hJetEta4After             -- Hist of Num of jets before cuts   *
 *    Int_t _JetEta4Min;                 -- Minimum Jet 4 Eta                 *
 *    Int_t _JetEta4Min;                 -- Minimum Jet 4 Eta                 *
 *                                                                            *
 * History                                                                    *
 *      24 Nov 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutJetEta4_h
#define CutJetEta4_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutJetEta4 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutJetEta4(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutJetEta4();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutJetEta4"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetEta4Before;    // Histogram Eta of jet 4 before the cut
  myTH1F* _hJetEta4After;     // Histogram Eta of jet 4 after the cut

  // Cut parameters
  Double_t _JetEta4Min;           // Minimum Eta of jet 4
  Double_t _JetEta4Max;           // Maximum Eta of jet 4

};


#endif







