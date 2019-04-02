/******************************************************************************
 * HistogrammingElectron.hpp                                                  *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Must pass a paramter describing the type of muon you want to plot          *
 *   Choices are All, Veto, Tight                                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingElectron()           -- Parameterized Constructor          *
 *    ~HistogrammingElectron()          -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingElectron"    *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of electron properties                               *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      18 Jan 2007 - P. Ryan - modified for all/tight/veto/                  *
 *****************************************************************************/

#ifndef HistogrammingElectron_h
#define HistogrammingElectron_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"



class HistogrammingElectron : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingElectron(EventContainer *obj, TString electronTypePassed);
  
  // Destructor
  ~HistogrammingElectron();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingElectron"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;            // Histogram of number of electrons
  
  myTH1F* _hPtObj1;          // Electron 1 PT
  myTH1F* _hPtBlowupObj1;    // Electron 1 PT, region 0-100GeV.
  myTH1F* _hEtaObj1;         // Electron 1 Eta
  myTH1F* _hPhiObj1;         // Electron 1 Phi
  myTH2F* _hPtEtaObj1;       // Electron 1 Pt-Eta
  myTH2F* _hPtPhiObj1;       // Electron 1 Pt-Eta
  myTH2F* _hEtaPhiObj1;      // Electron 1 Pt-Eta

  myTH1F* _hPtObj2;          // Electron 2 PT
  myTH1F* _hEtaObj2;         // Electron 2 Eta
  myTH1F* _hPhiObj2;         // Electron 2 Phi
  myTH2F* _hPtEtaObj2;       // Electron 2 Pt-Eta
  myTH2F* _hPtPhiObj2;       // Electron 2 Pt-Eta
  myTH2F* _hEtaPhiObj2;      // Electron 2 Pt-Eta

  myTH1F *_hPt;              // Pt of all electrons
  myTH1F *_hEta;             // Eta of all electrons
  myTH1F *_hPhi;             // Phi of all electrons
  myTH2F* _hPtEta;           // Electron all Pt-Eta
  myTH2F* _hPtPhi;           // Electron all Pt-Eta
  myTH2F* _hEtaPhi;          // Electron all Pt-Eta

  myTH1F *_hCharge;          // Charge of all electrons 
  myTH1F *_hEoverP;          // E over P of all electrons 
  myTH1F *_hEtCone30OverPt;     // calorimeter ET in a cone 30 over Pt for selected electrons 
  myTH1F *_hPtCone30OverPt;     // inner dectector track sum Pt in a cone 30 over Pt for selected electrons 

 myTH1F *_hJetDeltaPhi;
 myTH1F *_hJetDeltaR;

  TString electronType; // Type of electron being plotted

};


#endif







