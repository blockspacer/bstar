/******************************************************************************
 * Histogrammingneutrino.cpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingNeutrino()              -- Parameterized Constructor            *
 *    ~HistogrammingNeutrino()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingNeutrino"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 * - lots of histograms of w properties                                       *
 *                                                                            *
 * History                                                                    *
 *      15 May 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingNeutrino.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingNeutrino::HistogrammingNeutrino(EventContainer *obj)                    *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingNeutrino::HistogrammingNeutrino(EventContainer *obj)
{
  SetEventContainer(obj);
} //HistogrammingNeutrino()

/******************************************************************************
 * HistogrammingNeutrino::~HistogrammingNeutrino()                                      *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingNeutrino::~HistogrammingNeutrino()
{
  
} //HistogrammingNeutrino

/******************************************************************************
 * void HistogrammingNeutrino::BookHistogram()                                     *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingNeutrino::BookHistogram(){
  
  // Histogram of Number of Objects
  _hNObj = DeclareTH1F("neutrinoN","Number of neutrinos",10,0.,10.);
  _hNObj->SetXAxisTitle("N");
  _hNObj->SetYAxisTitle("Events");
  // neutrino 1 Pt
  _hPtObj1 = DeclareTH1F("neutrino1Pt","neutrino 1 P_{T}",200,0.,400.);
  _hPtObj1->SetXAxisTitle("P_{T} (neutrino 1) [GeV]");
  _hPtObj1->SetYAxisTitle("Events");
  // neutrino 1 Eta
  _hEtaObj1 = DeclareTH1F("neutrino1Eta","neutrino 1 Pseudorapidity",100,-5.,5.);
  _hEtaObj1->SetXAxisTitle("#eta (neutrino 1) ");
  _hEtaObj1->SetYAxisTitle("Events");
  // neutrino 1 Phi
  _hPhiObj1 = DeclareTH1F("neutrino1Phi","neutrino 1 phi",100,-1.*TMath::TwoPi(),TMath::TwoPi());
  _hPhiObj1->SetXAxisTitle("#phi (neutrino 1)");
  _hPhiObj1->SetYAxisTitle("Events");
  // neutrino 2 Pt
  _hPtObj2 = DeclareTH1F("neutrino2Pt","neutrino 2 P_{T}",200,0.,400.);
  _hPtObj2->SetXAxisTitle("P_{T} (neutrino 2) [GeV]");
  _hPtObj2->SetYAxisTitle("Events");
  // neutrino 2 Eta
  _hEtaObj2 = DeclareTH1F("neutrino2Eta","neutrino 2 Pseudorapidity",100,-5.,5.);
  _hEtaObj2->SetXAxisTitle("#eta (neutrino 2) ");
  _hEtaObj2->SetYAxisTitle("Events");
  // neutrino 2 Phi
  _hPhiObj2 = DeclareTH1F("neutrino2Phi","neutrino 2 phi",100,-1.*TMath::TwoPi(),TMath::TwoPi());
  _hPhiObj2->SetXAxisTitle("#phi (neutrino 2)");
  _hPhiObj2->SetYAxisTitle("Events");

  // any neutrino Pt
  _hPt = DeclareTH1F("neutrinoPt","any neutrino P_{T}",200,0.,400.);
  _hPt->SetXAxisTitle("P_{T} (neutrino) [GeV]");
  _hPt->SetYAxisTitle("Events");
 // any neutrino Pz
  _hPz = DeclareTH1F("neutrinoPz","any neutrino P_{z}",200,0.,400.);
  _hPz->SetXAxisTitle("P_{z} (neutrino) [GeV]");
  _hPz->SetYAxisTitle("Events");
  // any neutrino Eta
  _hEta = DeclareTH1F("neutrinoEta","any neutrino Pseudorapidity",100,-5.,5.);
  _hEta->SetXAxisTitle("#eta (neutrino) ");
  _hEta->SetYAxisTitle("Events");
  // any neutrino Phi
  _hPhi = DeclareTH1F("neutrinoPhi","Any neutrino phi",100,-1.*TMath::TwoPi(),TMath::TwoPi());
  _hPhi->SetXAxisTitle("#phi (neutrino)");
  _hPhi->SetYAxisTitle("Events");
  // any neutrino Charge
  _hCharge = DeclareTH1F("neutrinoCharge","Any neutrino charge",3,-1.5,1.5);
  _hCharge->SetXAxisTitle("charge (neutrino)");
  _hCharge->SetYAxisTitle("Events");
// any neutrino E
  _hE = DeclareTH1F("neutrinoE","any neutrino E",200,0.,400.);
  _hE->SetXAxisTitle("E (neutrino) [GeV]");
  _hE->SetYAxisTitle("Events");
// any neutrino M
  _hM = DeclareTH1F("neutrinoM","any neutrino M",200,0.,400.);
  _hM->SetXAxisTitle("M (neutrino) [GeV]");
  _hM->SetYAxisTitle("Events");

  
  //cout<<"end of HistogrammingNeutrino::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingNeutrino::Apply()                                           *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingNeutrino::Apply()
{
  //cout<<"Begin of HistogrammingNeutrino::Apply()"<<endl;
  
  // Get Event Tree
  EventContainer *evc = GetEventContainer();
  
  // Fill Histograms
  _hNObj -> Fill(evc->neutrinos.size());

  if(evc->neutrinos.size()>0) {
    _hPtObj1  -> Fill(evc->neutrinos[0].Pt());
    _hEtaObj1 -> Fill(evc->neutrinos[0].Eta());
    _hPhiObj1 -> Fill(evc->neutrinos[0].Phi());
  }
  if(evc->neutrinos.size()>1) {
    _hPtObj2  -> Fill(evc->neutrinos[1].Pt());
    _hEtaObj2 -> Fill(evc->neutrinos[1].Eta());
    _hPhiObj2 -> Fill(evc->neutrinos[1].Phi());
  }

  // now loop over all neutrinos and histogram their properties
  for(int io=0;io<evc->neutrinos.size();io++) {
    _hPt  -> Fill(evc->neutrinos[io].Pt());
    _hPz  -> Fill(evc->neutrinos[io].Pz());
    _hEta -> Fill(evc->neutrinos[io].Eta());
    _hPhi -> Fill(evc->neutrinos[io].Phi());
    _hCharge->Fill(evc->neutrinos[io].GetCharge());
    _hE  -> Fill(evc->neutrinos[io].E());
    _hM  -> Fill(evc->neutrinos[io].M());

   }    

  //cout<<"End of HistogrammingNeutrino::Apply()"<<endl;
  return kTRUE;  
  
} //Apply

