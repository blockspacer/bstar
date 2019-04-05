/******************************************************************************
 * HistogrammingMCTop.hpp                                                     *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMCTop()              -- Parameterized Constructor          *
 *    ~HistogrammingMCTop()             -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMCTop"       *
 *                                                                            *
 * Private Data Members of HistogrammingMCTop class                           *
 * - lots of histograms of top propreties                                     *
 *                                                                            *
 * History                                                                    *
 *      15 May 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/

#ifndef HistogrammingMCTop_h
#define HistogrammingMCTop_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMCTop : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMCTop(EventContainer *obj);
  
  // Destructor
  ~HistogrammingMCTop();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMCTop"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of tops
  myTH1F* _hPtObj1;   // Top 1 PT
  myTH1F* _hPtObj2;   // Top 2 PT
  myTH1F* _hEtaObj1;  // Top 1 Eta
  myTH1F* _hEtaObj2;  // Top 2 Eta
  myTH1F* _hPhiObj1;  // Top 1 Phi
  myTH1F* _hPhiObj2;  // Top 2 Phi

  myTH1F *_hPt;       // Pt of all tops
  myTH1F *_hEta;      // Eta of all tops
  myTH1F *_hPhi;      // Phi of all tops
  myTH1F *_hCharge;   // Charge of all tops 
  myTH1F *_hTopPdgId;  //PDG ID of all tops

};


#endif




