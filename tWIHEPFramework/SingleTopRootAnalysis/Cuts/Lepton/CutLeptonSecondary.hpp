/******************************************************************************
 * CutLeptonSecondary.hpp                                                     *
 *                                                                            *
 * Cut on Secondary Lepton                                                    *
 * Rejects lepton (electorn or muon) in central region                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutLeptonSecondary class                        *
 *    CutLeptonSecondary()              -- Parameterized Constructor          *
 *    ~CutLeptonSecondary()             -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutLeptonSecondary"       *
 *                                                                            *
 * Private Data Members of CutLeptonSecondary class                           *
 *    myTH1F* _hElectronPt1Before;      // Hist Ele 1 p_T before cut          *
 *    myTH1F* _hElectronPt1After;       // Hist Ele1 p_T of jets after cut    *
 *                                                                            *
 *    myTH1F* _hMuonPt1Before;          // Hist Muon 1 p_T before cut         *
 *    myTH1F* _hMuonPt1After;           // Hist Muon 1 p_T of jets after cut  *
 *                                                                            *
 *    Int_t _ElectronNumberMin;         // Minimum Electron Number            *
 *    Int_t _ElectronNumberMax;         // Maximum Electron Number            *
 *                                                                            *
 *    Int_t _MuonNumberMin;             // Minimum Muon Number                *
 *    Int_t _MuonNumberMax;             // Maximum Muon Number                *
 *                                                                            *
 *    Double_t _ElectronPt1Min;         // Minimum Electron 1 p_T             *
 *    Double_t _ElectronPt1Max;         // Maximum Electron 1 p_T             *
 *                                                                            *
 *    Double_t _MuonPt1Min;             // Minimum Muon 1 p_T                 *
 *    Double_t _MuonPt1Max;             // Maximum Muon 1 p_T                 *
 *                                                                            *
 * History                                                                    *
 *      26 Feb 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutLeptonSecondary_h
#define CutLeptonSecondary_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutLeptonSecondary : public HistoCut 
{

public:

  // Parameterized Constructor
  CutLeptonSecondary(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutLeptonSecondary();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutLeptonSecondary"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hElectronEtaBefore;      // Histogram electron eta before the cut
  myTH1F* _hElectronEtaAfter;       // Histogram electron eta after the cut

  myTH1F* _hMuonEtaBefore;          // Histogram muon eta before the cut
  myTH1F* _hMuonEtaAfter;           // Histogram muon eta after the cut

  myTH1F* _hElectronNBefore;        // Histogram electron number before the cut
  myTH1F* _hElectronNAfter;         // Histogram electron number after the cut

  myTH1F* _hMuonNBefore;            // Histogram muon number before the cut
  myTH1F* _hMuonNAfter;             // Histogram muon number after the cut

  myTH1F* _hElectronPtBefore;       // Histogram electron pT before the cut
  myTH1F* _hElectronPtAfter;        // Histogram electron pT after the cut

  myTH1F* _hMuonPtBefore;       // Histogram electron pT before the cut
  myTH1F* _hMuonPtAfter;        // Histogram electron pT after the cut

  myTH1F* _hEventsBefore;           // Histogram events before the cut
  myTH1F* _hEventsAfter;            // Histogram events after the cut

  // Cut parameters
  Double_t _ElectronPt;             // Secondary Electron pT
  Double_t _MuonPt;                 // Secondary Muon pT
  Double_t _CutElectron;            // 1 if high pt electron sample, 0 if high pt muon sample

  Double_t _ElectronEtaMin;         // Minimum Electron Eta for rejection range
  Double_t _ElectronEtaMax;         // Maximum Electron Eta for rejection range

  Double_t _MuonEtaMin;             // Minimum Muon Eta for rejection range
  Double_t _MuonEtaMax;             // Maximum Muon Eta for rejection range
  
  Int_t _ElectronNMax;              // Max number of electrons allowed in event
  Int_t _MuonNMax;                  // Max number of muons allowed in event


};


#endif







