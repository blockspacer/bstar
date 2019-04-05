/******************************************************************************
 * CutLeptonPt1.hpp                                                           *
 *                                                                            *
 * Cuts on lepton pT                                                          *
 * Cut on Electron pT *OR* Muon pT                                            *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutLeptonPt1 class                              *
 *    CutLeptonPt1()                    -- Parameterized Constructor          *
 *    ~CutLeptonPt1()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutLeptonPt1"             *
 *                                                                            *
 * Private Data Members of CutLeptonPt1 class                                 *
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
 *      23 Feb 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutLeptonPt1_h
#define CutLeptonPt1_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutLeptonPt1 : public HistoCut 
{

public:

  // Parameterized Constructor
  CutLeptonPt1(EventContainer *EventContainerObj, TString leptonType);
  
  // Destructor
  ~CutLeptonPt1();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutLeptonPt1"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  //  myTH1F* _hElectronPt1Before;      // Histogram Electron 1 p_T before the cut
  //myTH1F* _hElectronPt1After;       // Histogram Electron 1 p_T of jets after the cut

  //myTH1F* _hMuonPt1Before;          // Histogram Muon 1 p_T before the cut
  //  myTH1F* _hMuonPt1After;           // Histogram Muon 1 p_T of jets after the cut

  // The type of lepton we are using
  TString _leptonType;

  myTH1F* _hLeadingLeptonPtBefore;  // Histogram of leading pt before cut
  myTH1F* _hLeadingLeptonPtAfter;   // Histogram of leading pt after cut

  // Cut parameters
  Double_t _LeadingMuonPtCut;
  Double_t _LeadingElectronPtCut;

  /*Int_t _ElectronNumberMin;         // Minimum Electron Number
  Int_t _ElectronNumberMax;         // Maximum Electron Number

  Int_t _MuonNumberMin;             // Minimum Muon Number
  Int_t _MuonNumberMax;             // Maximum Muon Number

  Double_t _ElectronPt1Min;         // Minimum Electron 1 p_T
  Double_t _ElectronPt1Max;         // Maximum Electron 1 p_T

  Double_t _MuonPt1Min;             // Minimum Muon 1 p_T
  Double_t _MuonPt1Max;             // Maximum Muon 1 p_T
  */

};


#endif







