/******************************************************************************
 * CutJetEta3.hpp                                                             *
 *                                                                            *
 * Cuts on Eta of 3st Jet                                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetEta3 class                                *
 *    CutJetEta3()                      -- Parameterized Constructor          *
 *    ~CutJetEta3()                     -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetEta3"               *
 *                                                                            *
 * Private Data Members of CutJetEta3                                         *
 *    myTH1F* _hJetEta3Before            -- Hist of Num of jets before cuts   *
 *    myTH1F* _hJetEta3After             -- Hist of Num of jets before cuts   *
 *    Int_t _JetEta3Min;                 -- Minimum Jet 3 Eta                 *
 *    Int_t _JetEta3Min;                 -- Minimum Jet 3 Eta                 *
 *                                                                            *
 * History                                                                    *
 *      24 Nov 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutJetEta3_h
#define CutJetEta3_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutJetEta3 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutJetEta3(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutJetEta3();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutJetEta3"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetEta3Before;    // Histogram Eta of jet 3 before the cut
  myTH1F* _hJetEta3After;     // Histogram Eta of jet 3 after the cut

  // Cut parameters
  Double_t _JetEta3Min;           // Minimum Eta of jet 3
  Double_t _JetEta3Max;           // Maximum Eta of jet 3

};


#endif







