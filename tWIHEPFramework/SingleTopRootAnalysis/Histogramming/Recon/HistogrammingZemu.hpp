/******************************************************************************
 * HistogrammingZemu.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Must pass a paramter describing the type of muon you want to plot          *
 *   Choices are All, Veto, Tight                                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingZemu()           -- Parameterized Constructor               *
 *    ~HistogrammingZemu()          -- Destructor                              *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingZemu"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of electron properties                               *
 *                                                                            *
 * History                                                                    *
 *      21 Nov 2010 - Created by H. Zhang for Wt channel Drell-Yan studies    *
 *****************************************************************************/

#ifndef HistogrammingZemu_h
#define HistogrammingZemu_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include "TLorentzVector.h"
#include "TMath.h"

class HistogrammingZemu : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingZemu(EventContainer *obj, TString electronTypePassed, TString muonTypePassed, bool doPileupWeight=false);
  
  // Destructor
  ~HistogrammingZemu();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingZemu"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();


  
private:

  double xy_to_phi(double,double);

  double dPhi(double, double);

  void soulveZtautau(TLorentzVector lep1, TLorentzVector lep2, double METX, double METY, TLorentzVector & nv1, TLorentzVector & nv2);

  bool _doPileupWeight;
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

  myTH1F *_hPtu;              // Pt of all electrons
  myTH1F *_hEtau;             // Eta of all electrons
  myTH1F *_hPhiu;             // Phi of all electrons
  myTH2F* _hPtEtau;           // Electron all Pt-Eta
  myTH2F* _hPtPhiu;           // Electron all Pt-Eta
  myTH2F* _hEtaPhiu;          // Electron all Pt-Eta

  myTH1F *_hMZemu_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZemu_0jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZemu_1jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZemu_2jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZemu_3jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZemu_4jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZemu_1plusjet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZemu_2plusjet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZemu_3plusjet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZemu_SS;         // invairant mass of two same sign electrons
  myTH1F *_hMZemu;            // invairant mass of any two electrons
  myTH1F *_hM_METcoLinearemu_OS;         // invairant mass of any two leptons plus coLinear solved neutrinos
  myTH1F *_hMET;             // MET of this events
  myTH2F* _hMemuMET;          // Memu vs MET
  myTH2F* _hM_METcoLinearemu_MET_OS;          // M_METcoLinearemu vs MET
  myTH2F* _hMemu_OSMET;       // Opposite sign Memu vs MET
  myTH2F* _hMemu_0jet_OSMET;  // Opposite sign Memu vs MET
  myTH2F* _hMemu_1jet_OSMET;  // Opposite sign Memu vs MET
  myTH2F* _hMemu_2jet_OSMET;  // Opposite sign Memu vs MET
  myTH2F* _hMemu_3jet_OSMET;  // Opposite sign Memu vs MET
  myTH2F* _hMemu_4jet_OSMET;  // Opposite sign Memu vs MET
  myTH2F* _hMemu_1plusjet_OSMET;  // Opposite sign Memu vs MET
  myTH2F* _hMemu_2plusjet_OSMET;  // Opposite sign Memu vs MET
  myTH2F* _hMemu_3plusjet_OSMET;  // Opposite sign Memu vs MET
  myTH2F* _hMemu_SSMET;       // Same sign Memu vs MET
  myTH2F* _hMemu_OSdPhiMETJ1;   // Memu_OS vs dPhi MET-leading Jet
  myTH2F* _hMemu_OSdPhiMETE1;   // Memu_OS vs dPhi MET-leading elec
  myTH2F* _hMemu_OSdPhiMETZt;   // Memu_OS vs dPhi MET-Z
  myTH2F* _hMemu_OSdPhiZtJet;   // Memu_OS vs dPhi Z-Jet

  myTH2F* _hMemu_OSMET_projZt;   // Memu vs MET project along Z
  myTH2F* _hMemu_OSMET_projE1;   // Memu vs MET project along leading Electron
  myTH2F* _hMemu_OSMET_projJ1;   // Memu vs MET project along leading Jet

  myTH2F* _hMETdPhiMETZt;   // MET vs dPhi MET-Z
  myTH2F* _hMETdPhiMETEl;   // MET vs dPhi MET-leading Electron
  myTH2F* _hMETdPhiMETJ1;   // MET vs dPhi MET-leading Jet
  myTH2F* _hMETdPhiZtJet;   // MET vs dPhi Z-Jet
  
  myTH1F *_hmissingPhi;      // missing Phi
  myTH1F *_hCharge;          // Charge of all electrons 
  myTH1F *_hChargeu;          // Charge of all electrons 
  myTH1F *_hEoverP;          // E over P of all electrons 

  myTH1F *_hJetDeltaPhi;
  myTH1F *_hJetDeltaR;

  TString electronType; // Type of electron being plotted
  TString muonType; // Type of electron being plotted

};


#endif







