/*****************************************************************************
 * HistogrammingMCMuon.hpp                                                    *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMCMuon()              -- Parameterized Constructor         *
 *    ~HistogrammingMCMuon()             -- Destructor                        *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMCMuon"      *
 *                                                                            *
 * Private Data Members of HistogrammingMCMuon class                          *
 * - lots of histograms of muon propreties                                    *
 *                                                                            *
 * History                                                                    *
 *      07 Feb 2007 - Created by J. Holzbauer to plot MCMuon Histograms       *
 *      25 Apr 2007 - Modified by J. Holzbauer for DeltaR, PDG ID plots       *
 *****************************************************************************/

#ifndef HistogrammingMCMuon_h
#define HistogrammingMCMuon_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMCMuon : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMCMuon(EventContainer *obj);
  
  // Destructor
  ~HistogrammingMCMuon();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMCMuon"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of muons
  myTH1F* _hPtObj1;   // Muon 1 PT
  myTH1F* _hPtObj2;   // Muon 2 PT
  myTH1F* _hEtaObj1;  // Muon 1 Eta
  myTH1F* _hEtaObj2;  // Muon 2 Eta
  myTH1F* _hPhiObj1;  // Muon 1 Phi
  myTH1F* _hPhiObj2;  // Muon 2 Phi

  myTH1F *_hPt;       // Pt of all muons
  myTH1F *_hEta;      // Eta of all muons
  myTH1F *_hPhi;      // Phi of all muons
  myTH1F *_hCharge;   // Charge of all muons 
  myTH1F *_hMuPdgId;  //PDG ID of all muons
  myTH1F *_hdeltaRtree;  //DeltaR between MC and Reco Muons, from tree
  myTH2F *_hdeltaEta;  //DeltaR between MC and Reco Muons, from tree

  myTH1F* _hdeltaRTagged;
 
};


#endif






