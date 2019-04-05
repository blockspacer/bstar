/******************************************************************************
 * HistogrammingZee.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Must pass a paramter describing the type of muon you want to plot          *
 *   Choices are All, Veto, Tight                                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingZee()           -- Parameterized Constructor               *
 *    ~HistogrammingZee()          -- Destructor                              *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingZee"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of electron properties                               *
 *                                                                            *
 * History                                                                    *
 *      21 Nov 2010 - Created by H. Zhang for Wt channel Drell-Yan studies    *
 *****************************************************************************/

#ifndef HistogrammingZee_h
#define HistogrammingZee_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"



class HistogrammingZee : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingZee(EventContainer *obj, TString electronTypePassed, bool doPileupWeight=false);
  
  // Destructor
  ~HistogrammingZee();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingZee"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
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

  myTH1F *_hMZee_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZee_0jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZee_1jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZee_2jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZee_3jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZee_4jet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZee_1plusjet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZee_2plusjet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZee_3plusjet_OS;         // invairant mass of two opposite sign electrons
  myTH1F *_hMZee_SS;         // invairant mass of two same sign electrons
  myTH1F *_hMZee;            // invairant mass of any two electrons
  myTH1F *_hMET;             // MET of this events
  myTH2F* _hMeeMET;          // Mee vs MET
  myTH2F* _hMee_OSMET;       // Opposite sign Mee vs MET
  myTH2F* _hMee_0jet_OSMET;       // Opposite sign Mee vs MET
  myTH2F* _hMee_1jet_OSMET;       // Opposite sign Mee vs MET
  myTH2F* _hMee_2jet_OSMET;       // Opposite sign Mee vs MET
  myTH2F* _hMee_3jet_OSMET;       // Opposite sign Mee vs MET
  myTH2F* _hMee_4jet_OSMET;       // Opposite sign Mee vs MET
  myTH2F* _hMee_1plusjet_OSMET;       // Opposite sign Mee vs MET
  myTH2F* _hMee_2plusjet_OSMET;       // Opposite sign Mee vs MET
  myTH2F* _hMee_3plusjet_OSMET;       // Opposite sign Mee vs MET
  myTH2F* _hMee_SSMET;       // Same sign Mee vs MET
  myTH2F* _hMee_OSdPhiMETJ1;   // Mee_OS vs dPhi MET-leading Jet
  myTH2F* _hMee_OSdPhiMETE1;   // Mee_OS vs dPhi MET-leading elec
  myTH2F* _hMee_OSdPhiMETZt;   // Mee_OS vs dPhi MET-Z
  myTH2F* _hMee_OSdPhiZtJet;   // Mee_OS vs dPhi Z-Jet

  myTH2F* _hMee_OSMET_projZt;   // Mee vs MET project along Z
  myTH2F* _hMee_OSMET_projE1;   // Mee vs MET project along leading Electron
  myTH2F* _hMee_OSMET_projJ1;   // Mee vs MET project along leading Jet

  myTH2F* _hMETdPhiMETZt;   // MET vs dPhi MET-Z
  myTH2F* _hMETdPhiMETEl;   // MET vs dPhi MET-leading Electron
  myTH2F* _hMETdPhiMETJ1;   // MET vs dPhi MET-leading Jet
  myTH2F* _hMETdPhiZtJet;   // MET vs dPhi Z-Jet
  
  myTH1F *_hmissingPhi;      // missing Phi
  myTH1F *_hCharge;          // Charge of all electrons 
  myTH1F *_hEoverP;          // E over P of all electrons 

  myTH1F *_hJetDeltaPhi;
  myTH1F *_hJetDeltaR;

  TString electronType; // Type of electron being plotted

};


#endif







