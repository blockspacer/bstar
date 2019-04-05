/******************************************************************************
 * HistogrammingZmumu.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Must pass a paramter describing the type of muon you want to plot          *
 *   Choices are All, Veto, Tight                                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingZmm()           -- Parameterized Constructor               *
 *    ~HistogrammingZmm()          -- Destructor                              *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingZmm"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of mueon    properties                               *
 *                                                                            *
 * History                                                                    *
 *      10 Dec 2010 - Created by J. Koll for Wt channel Drell-Yan studies    *
 *****************************************************************************/

#ifndef HistogrammingZmumu_h
#define HistogrammingZmumu_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"



class HistogrammingZmumu : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingZmumu(EventContainer *obj, TString muonTypePassed, bool doPileupWeight=false);
  
  // Destructor
  ~HistogrammingZmumu();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingZmumu"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();

  double xy_to_phi(double,double);
 
  double dPhi(double, double);
  
private:
  bool _doPileupWeight;
  // Histograms declarations 
  myTH1F* _hNObj;            // Histogram of number of muons
  
  myTH1F* _hPtObj1;          // Muon 1 PT
  myTH1F* _hPtBlowupObj1;    // Muon 1 PT, region 0-100GeV.
  myTH1F* _hEtaObj1;         // Muon 1 Eta
  myTH1F* _hPhiObj1;         // Muon 1 Phi
  myTH2F* _hPtEtaObj1;       // Muon 1 Pt-Eta
  myTH2F* _hPtPhiObj1;       // Muon 1 Pt-Eta
  myTH2F* _hEtaPhiObj1;      // Muon 1 Pt-Eta

  myTH1F* _hPtObj2;          // Muon 2 PT
  myTH1F* _hEtaObj2;         // Muon 2 Eta
  myTH1F* _hPhiObj2;         // Muon 2 Phi
  myTH2F* _hPtEtaObj2;       // Muon 2 Pt-Eta
  myTH2F* _hPtPhiObj2;       // Muon 2 Pt-Eta
  myTH2F* _hEtaPhiObj2;      // Muon 2 Pt-Eta

  myTH1F *_hPt;              // Pt of all muons
  myTH1F *_hEta;             // Eta of all muons
  myTH1F *_hPhi;             // Phi of all muons
  myTH2F* _hPtEta;           // Muon all Pt-Eta
  myTH2F* _hPtPhi;           // Muon all Pt-Eta
  myTH2F* _hEtaPhi;          // Muon all Pt-Eta

  myTH1F *_hMZmumu_OS;         // invairant mass of two opposite sign muons
  myTH1F *_hMZmumu_0jet_OS;    // invairant mass of two opposite sign muons
  myTH1F *_hMZmumu_1jet_OS;    // invairant mass of two opposite sign muons
  myTH1F *_hMZmumu_2jet_OS;    // invairant mass of two opposite sign muons
  myTH1F *_hMZmumu_3jet_OS;    // invairant mass of two opposite sign muons
  myTH1F *_hMZmumu_4jet_OS;    // invairant mass of two opposite sign muons
  myTH1F *_hMZmumu_1plusjet_OS;    // invairant mass of two opposite sign muons
  myTH1F *_hMZmumu_2plusjet_OS;    // invairant mass of two opposite sign muons
  myTH1F *_hMZmumu_3plusjet_OS;    // invairant mass of two opposite sign muons
  myTH1F *_hMZmumu_SS;         // invairant mass of two same sign muons
  myTH1F *_hMZmumu;            // invairant mass of any two muons
  myTH1F *_hMET;             // MET of this events
  myTH2F* _hMmumuMET;          // Mmumu vs MET
  myTH2F* _hMmumu_OSMET;       // Opposite sign Mmumu vs MET
  myTH2F* _hMmumu_0jet_OSMET;  // Opposite sign Mmumu vs MET
  myTH2F* _hMmumu_1jet_OSMET;  // Opposite sign Mmumu vs MET
  myTH2F* _hMmumu_2jet_OSMET;  // Opposite sign Mmumu vs MET
  myTH2F* _hMmumu_3jet_OSMET;  // Opposite sign Mmumu vs MET
  myTH2F* _hMmumu_4jet_OSMET;  // Opposite sign Mmumu vs MET
  myTH2F* _hMmumu_1plusjet_OSMET;  // Opposite sign Mmumu vs MET
  myTH2F* _hMmumu_2plusjet_OSMET;  // Opposite sign Mmumu vs MET
  myTH2F* _hMmumu_3plusjet_OSMET;  // Opposite sign Mmumu vs MET
  myTH2F* _hMmumu_SSMET;       // Same sign Mmumu vs MET
  myTH2F* _hMmumudPhiMETJ1;   // Mmumu vs dPhi MET-leading Jet
  myTH2F* _hMmumudPhiMETE1;   // Mmumu vs dPhi MET-leading muon
  myTH2F* _hMmumudPhiMETZt;   // Mmumu vs dPhi MET-Z

  myTH1F *_hmissingPhi;      // missing Phi
  myTH1F *_hCharge;          // Charge of all muons 
  myTH1F *_hEoverP;          // E over P of all muons 

  myTH1F *_hJetDeltaPhi;
  myTH1F *_hJetDeltaR;

  TString muonType; // Type of muon being plotted

};


#endif







