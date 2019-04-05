/******************************************************************************
 * CutJetAll.hpp                                                              *
 *                                                                            *
 * Selects events with the correct number of jets                             *
 * Takes into account the pT and eta requirements used in jet definition      *
 *                                                                            *
 * Differs from an event selection that requires, for example,                *
 *     2 <= N_Jets <= 4                                                       *
 *     pT_1 > 30, pT_2 > 30, pT_3 > 15, pT_4 > 15                             *
 *     -2.5 < Eta < 2.5                                                       *
 * where the highest pT cut is applied only to the 1st jet in a pT ordered    *
 * list.  It could be possible that this jet fails the eta requirement but    *
 * the second jet in the list satisfies the requirement of the highest pT jet *
 *                                                                            *
 * Also, cutting on the 3rd and 4th jet pT (in a linear event selection like  *
 * this program performs, the existance of a 3rd and 4th jet is implicitly    *
 * required.                                                                  *
 *                                                                            *
 * This Cut Class counts the number of valid jets based upon the pT and Eta   *
 * requirements.  Jets in the list failing to meet the pT and Eta             *
 * requirements are not counted as jets.                                      *
 *                                                                            *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutJetAll class                                 *
 *    CutJetAll()                       -- Parameterized Constructor          *
 *    ~CutJetAll()                      -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetAll"                *
 *                                                                            *
 * Private Data Members of CutJetAll                                          *
 *    myTH1F* _hJetAllBefore            -- Hist of Num of jets before cuts    *
 *    myTH1F* _hJetAllAfter             -- Hist of Num of jets before cuts    *
 *    Double_t _JetAllMin               -- Minimum Jet 1 Pt                   *
 *    Double_t _JetAllMax               -- Maximum Jet 1 Pt                   *
 *                                                                            *
 * History                                                                    *
 *      22 Nov 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutJetAll_h
#define CutJetAll_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutJetAll : public HistoCut 
{

public:

  // Parameterized Constructor
  CutJetAll(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutJetAll();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutJetAll"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hJetListNumberBefore;    // Histogram of number of jets in array before cut
  myTH1F* _hJetListNumberAfter;     // Histogram of number of jets in array after cut

  myTH1F* _hJetNumberBefore;        // Number of jets before cut - total
  myTH1F* _hJetNumberAfter;         // Number of jets after cut  - total
  myTH1F* _hJetNumberValidBefore;   // Number of jets before cut - pass validity requirement
  myTH1F* _hJetNumberValidAfter;    // Number of jets after cut  - pass validity requirement

  myTH1F* _hJetPtBefore;        // Pt of jets before cut - total
  myTH1F* _hJetPtAfter;         // Pt of jets after cut  - total
  myTH1F* _hJetPtValidBefore;   // Pt of jets before cut - pass validity requirement
  myTH1F* _hJetPtValidAfter;    // Pt of jets after cut  - pass validity requirement

  myTH1F* _hJetEtaBefore;        // Eta of jets before cut - total
  myTH1F* _hJetEtaAfter;         // Eta of jets after cut  - total
  myTH1F* _hJetEtaValidBefore;   // Eta of jets before cut - pass validity requirement
  myTH1F* _hJetEtaValidAfter;    // Eta of jets after cut  - pass validity requirement

  myTH1F* _hJetIndexBefore;        // Index of jets before cut - total
  myTH1F* _hJetIndexAfter;         // Index of jets after cut  - total
  myTH1F* _hJetIndexValidBefore;   // Index of jets before cut - pass validity requirement
  myTH1F* _hJetIndexValidAfter;    // Index of jets after cut  - pass validity requirement


  // Cut parameters
  Int_t _JetAllNumberMin;    // Minimum Jet Number (from parameter file)
  Int_t _JetAllNumberMax;    // Maximum Jet Number (from parameter file)

  Double_t _JetPtMinDefault;  // Default jet pT min
  Double_t _JetPtMaxDefault;  // Default jet pT max
  Double_t _JetEtaMinDefault; // Default jet Eta min
  Double_t _JetEtaMaxDefault; // Default jet Eta max

  //TMP
  Double_t _MuonPt1Min;       // Minimum Pt of muon 1

  Int_t JetNumberMin;       // Minimum Jet Number used in looping
  Int_t JetNumberMax;       // Maximum Jet Number used in looping

  std::vector<Double_t> _JetPtMinVector;   // Jet Pt Cuts for all jets being considered
  std::vector<Double_t> _JetPtMaxVector;   // Jet Pt Cuts for all jets being considered
  std::vector<Double_t> _JetEtaMinVector;  // Jet Eta Cuts for all jets being considered
  std::vector<Double_t> _JetEtaMaxVector;  // Jet Eta Cuts for all jets being considered


};


#endif







