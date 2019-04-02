/******************************************************************************
 * HistogrammingTau.cpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTau()              -- Parameterized Constructor            *
 *    ~HistogrammingTau()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTau"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of tau properties                                    *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingTau.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingTau::HistogrammingTau(EventContainer *obj)          *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingTau::HistogrammingTau(EventContainer *obj)
{
  SetEventContainer(obj);
} //HistogrammingTau()

/******************************************************************************
 * HistogrammingTau::~HistogrammingTau()                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingTau::~HistogrammingTau()
{
  
} //HistogrammingTau

/******************************************************************************
 * void HistogrammingTau::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingTau::BookHistogram(){
  
  // Histogram of Number of Objects
  _hNObj = DeclareTH1F("TauN","Number of Taus",10,0.,10.);
  _hNObj->SetXAxisTitle("N");
  _hNObj->SetYAxisTitle("Events");
  // Tau 1 Pt
  _hPtObj1 = DeclareTH1F("Tau1Pt","Tau 1 P_{T}",200,0.,400.);
  _hPtObj1->SetXAxisTitle("P_{T} (tau 1) [GeV]");
  _hPtObj1->SetYAxisTitle("Events");
  // Tau 1 Eta
  _hEtaObj1 = DeclareTH1F("Tau1Eta","Tau 1 Pseudorapidity",100,-5.,5.);
  _hEtaObj1->SetXAxisTitle("#eta (tau 1) ");
  _hEtaObj1->SetYAxisTitle("Events");
  // Tau 1 Phi
  _hPhiObj1 = DeclareTH1F("Tau1Phi","Tau 1 phi",100,-1.*TMath::TwoPi(),TMath::TwoPi());
  _hPhiObj1->SetXAxisTitle("#phi (tau 1)");
  _hPhiObj1->SetYAxisTitle("Events");
  // Tau 2 Pt
  _hPtObj2 = DeclareTH1F("Tau2Pt","Tau 2 P_{T}",200,0.,400.);
  _hPtObj2->SetXAxisTitle("P_{T} (tau 2) [GeV]");
  _hPtObj2->SetYAxisTitle("Events");
  // Tau 2 Eta
  _hEtaObj2 = DeclareTH1F("Tau2Eta","Tau 2 Pseudorapidity",100,-5.,5.);
  _hEtaObj2->SetXAxisTitle("#eta (tau 2) ");
  _hEtaObj2->SetYAxisTitle("Events");
  // Tau 2 Phi
  _hPhiObj2 = DeclareTH1F("Tau2Phi","Tau 2 phi",100,-1.*TMath::TwoPi(),TMath::TwoPi());
  _hPhiObj2->SetXAxisTitle("#phi (tau 2)");
  _hPhiObj2->SetYAxisTitle("Events");

  // any Tau Pt
  _hPt = DeclareTH1F("TauPt","any Tau P_{T}",200,0.,400.);
  _hPt->SetXAxisTitle("P_{T} (tau) [GeV]");
  _hPt->SetYAxisTitle("Events");
  // any Tau Eta
  _hEta = DeclareTH1F("TauEta","any Tau Pseudorapidity",100,-5.,5.);
  _hEta->SetXAxisTitle("#eta (tau) ");
  _hEta->SetYAxisTitle("Events");
  // any Tau Phi
  _hPhi = DeclareTH1F("TauPhi","Any tau phi",100,-1.*TMath::TwoPi(),TMath::TwoPi());
  _hPhi->SetXAxisTitle("#phi (tau)");
  _hPhi->SetYAxisTitle("Events");
  // any Tau Charge
  _hCharge = DeclareTH1F("TauCharge","Any tau charge",3,-1.5,1.5);
  _hCharge->SetXAxisTitle("charge (tau)");
  _hCharge->SetYAxisTitle("Events");
  // any Tau number of tracks
  _hNumTracks = DeclareTH1F("TauNumTracks","Any tau number of tracks",10,0.,10.);
  _hNumTracks->SetXAxisTitle("number of tracks");
  _hNumTracks->SetYAxisTitle("Events");
  
  //cout<<"end of HistogrammingTau::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingTau::Apply()                                         *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingTau::Apply()
{
  //cout<<"Begin of HistogrammingTau::Apply()"<<endl;
  
  // Get Event Tree
  EventContainer *evc = GetEventContainer();
  
  // Fill Histograms
  _hNObj -> Fill(evc->taus.size());

  if(evc->taus.size()>0) {
    _hPtObj1  -> Fill(evc->taus[0].Pt());
    _hEtaObj1 -> Fill(evc->taus[0].Eta());
    _hPhiObj1 -> Fill(evc->taus[0].Phi());
  }
  if(evc->taus.size()>1) {
    _hPtObj2  -> Fill(evc->taus[1].Pt());
    _hEtaObj2 -> Fill(evc->taus[1].Eta());
    _hPhiObj2 -> Fill(evc->taus[1].Phi());
  }

  // now loop over all taus and histogram their properties
  for(int io=0;io<evc->taus.size();io++) {
    _hPt  -> Fill(evc->taus[io].Pt());
    _hEta -> Fill(evc->taus[io].Eta());
    _hPhi -> Fill(evc->taus[io].Phi());
    _hCharge->Fill(evc->taus[io].GetCharge());
    _hNumTracks->Fill(evc->taus[io].GetNumTracks());
  }    

  //cout<<"End of HistogrammingTau::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












