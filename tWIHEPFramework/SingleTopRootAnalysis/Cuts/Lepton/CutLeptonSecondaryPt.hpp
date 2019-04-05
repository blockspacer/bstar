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
 *      09 Aug - Created by J. Holzbauer based on CutLeptonSecondary by P.Ryan*
 *               Cuts all secondary leptons at same pt, and only pt cut       *
 *               Lepton = electron or isolated muon                           *
 *****************************************************************************/

#ifndef CutLeptonSecondaryPt_h
#define CutLeptonSecondaryPt_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutLeptonSecondaryPt : public HistoCut 
{

public:

  // Parameterized Constructor
  CutLeptonSecondaryPt(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutLeptonSecondaryPt();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutLeptonSecondaryPt"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hElectronPtBefore;      // Histogram electron pt before the cut
  myTH1F* _hElectronPtAfter;       // Histogram electron pt after the cut

  myTH1F* _hIsoMuonPtBefore;          // Histogram muon pt before the cut
  myTH1F* _hIsoMuonPtAfter;           // Histogram muon pt after the cut

  // Cut parameters

  Double_t _ElectronPtMin;         // Minimum Electron Pt for rejection range
  Double_t _ElectronPtMax;         // Maximum Electron Pt for rejection range

  Double_t _IsoMuonPtMin;             // Minimum IsoMuon Pt for rejection range
  Double_t _IsoMuonPtMax;             // Maximum IsoMuon Pt for rejection range

};


#endif







