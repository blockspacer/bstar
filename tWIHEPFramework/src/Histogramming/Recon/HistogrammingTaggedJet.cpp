/******************************************************************************
 * HistogrammingTaggedJet.cpp                                                 *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTaggedJet()          -- Parameterized Constructor          *
 *    ~HistogrammingTaggedJet()         -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTaggedJet"   *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of jet properties                                    *
 *                                                                            *
 * History                                                                    *
 *      21 Feb 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingTaggedJet.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingTaggedJet::HistogrammingTaggedJet(EventContainer *obj)        *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingTaggedJet::HistogrammingTaggedJet(EventContainer *obj)
{
  SetEventContainer(obj);
} //HistogrammingTaggedJet()

/******************************************************************************
 * HistogrammingTaggedJet::~HistogrammingTaggedJet()                          *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingTaggedJet::~HistogrammingTaggedJet()
{
  
} //HistogrammingTaggedJet

/******************************************************************************
 * void HistogrammingTaggedJet::BookHistogram()                               *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingTaggedJet::BookHistogram(){
  
  // Histogram of Number of Objects
  _hNObj = DeclareTH1F("JetN","Number of Jets",10,0.,10.);
  _hNObj->SetXAxisTitle("N^{Jet}");
  _hNObj->SetYAxisTitle("Events");

  // Jet 1 Pt
  _hPtObj1 = DeclareTH1F("Jet1Pt","Jet 1 p_{T}",100,0.,400.);
  _hPtObj1->SetXAxisTitle("p_{T}^{Jet 1} [GeV]");
  _hPtObj1->SetYAxisTitle("Events");
  // Jet 1 Eta
  _hEtaObj1 = DeclareTH1F("Jet1Eta","Jet 1 Pseudorapidity",50,-5.,5.);
  _hEtaObj1->SetXAxisTitle("#eta^{Jet 1}");
  _hEtaObj1->SetYAxisTitle("Events");
  // Jet 1 Phi
  _hPhiObj1 = DeclareTH1F("Jet1Phi","Jet 1 phi",50,-1.*TMath::Pi(),TMath::Pi());
  _hPhiObj1->SetXAxisTitle("#phi^{Jet 1}");
  _hPhiObj1->SetYAxisTitle("Events");

  // Jet 2 Pt
  _hPtObj2 = DeclareTH1F("Jet2Pt","Jet 2 p_{T}",100,0.,400.);
  _hPtObj2->SetXAxisTitle("p_{T}^{Jet 2} [GeV]");
  _hPtObj2->SetYAxisTitle("Events");
  // Jet 2 Eta
  _hEtaObj2 = DeclareTH1F("Jet2Eta","Jet 2 Pseudorapidity",50,-5.,5.);
  _hEtaObj2->SetXAxisTitle("#eta^{Jet 2} ");
  _hEtaObj2->SetYAxisTitle("Events");
  // Jet 2 Phi
  _hPhiObj2 = DeclareTH1F("Jet2Phi","Jet 2 phi",50,-1.*TMath::Pi(),TMath::Pi());
  _hPhiObj2->SetXAxisTitle("#phi^{Jet 2}");
  _hPhiObj2->SetYAxisTitle("Events");

  // Jet 3 Pt
  _hPtObj3 = DeclareTH1F("Jet3Pt","Jet 3 p_{T}",100,0.,400.);
  _hPtObj3->SetXAxisTitle("p_{T}^{Jet 3} [GeV]");
  _hPtObj3->SetYAxisTitle("Events");
  // Jet 3 Eta
  _hEtaObj3 = DeclareTH1F("Jet3Eta","Jet 3 Pseudorapidity",50,-5.,5.);
  _hEtaObj3->SetXAxisTitle("#eta^{Jet 3}");
  _hEtaObj3->SetYAxisTitle("Events");
  // Jet 3 Phi
  _hPhiObj3 = DeclareTH1F("Jet3Phi","Jet 3 phi",50,-1.*TMath::Pi(),TMath::Pi());
  _hPhiObj3->SetXAxisTitle("#phi^{Jet 3}");
  _hPhiObj3->SetYAxisTitle("Events");
  
  // Jet 4 Pt
  _hPtObj4 = DeclareTH1F("Jet4Pt","Jet 4 p_{T}",100,0.,400.);
  _hPtObj4->SetXAxisTitle("p_{T}^{Jet 4} [GeV]");
  _hPtObj4->SetYAxisTitle("Events");
  // Jet 4 Eta
  _hEtaObj4 = DeclareTH1F("Jet4Eta","Jet 4 Pseudorapidity",50,-5.,5.);
  _hEtaObj4->SetXAxisTitle("#eta^{Jet 4}");
  _hEtaObj4->SetYAxisTitle("Events");
  // Jet 4 Phi
  _hPhiObj4 = DeclareTH1F("Jet4Phi","Jet 4 phi",50,-1.*TMath::Pi(),TMath::Pi());
  _hPhiObj4->SetXAxisTitle("#phi^{Jet 4}");
  _hPhiObj4->SetYAxisTitle("Events");

  // Jet 5 Pt
  _hPtObj5 = DeclareTH1F("Jet5Pt","Jet 5 p_{T}",100,0.,400.);
  _hPtObj5->SetXAxisTitle("p_{T}^{Jet 5} [GeV]");
  _hPtObj5->SetYAxisTitle("Events");
  // Jet 5 Eta
  _hEtaObj5 = DeclareTH1F("Jet5Eta","Jet 5 Pseudorapidity",50,-5.,5.);
  _hEtaObj5->SetXAxisTitle("#eta^{Jet 5}");
  _hEtaObj5->SetYAxisTitle("Events");
  // Jet 5 Phi
  _hPhiObj5 = DeclareTH1F("Jet5Phi","Jet 5 phi",50,-1.*TMath::Pi(),TMath::Pi());
  _hPhiObj5->SetXAxisTitle("#phi^{Jet 5}");
  _hPhiObj5->SetYAxisTitle("Events");

  // Jet 6 Pt
  _hPtObj6 = DeclareTH1F("Jet6Pt","Jet 6 p_{T}",100,0.,400.);
  _hPtObj6->SetXAxisTitle("p_{T}^{Jet 6} [GeV]");
  _hPtObj6->SetYAxisTitle("Events");
  // Jet 6 Eta
  _hEtaObj6 = DeclareTH1F("Jet6Eta","Jet 6 Pseudorapidity",50,-5.,5.);
  _hEtaObj6->SetXAxisTitle("#eta^{Jet 6}");
  _hEtaObj6->SetYAxisTitle("Events");
  // Jet 6 Phi
  _hPhiObj6 = DeclareTH1F("Jet6Phi","Jet 6 phi",50,-1.*TMath::Pi(),TMath::Pi());
  _hPhiObj6->SetXAxisTitle("#phi^{Jet 6}");
  _hPhiObj6->SetYAxisTitle("Events");

  // Jet 7 Pt
  _hPtObj7 = DeclareTH1F("Jet7Pt","Jet 7 p_{T}",100,0.,400.);
  _hPtObj7->SetXAxisTitle("p_{T}^{Jet 7} [GeV]");
  _hPtObj7->SetYAxisTitle("Events");
  // Jet 7 Eta
  _hEtaObj7 = DeclareTH1F("Jet7Eta","Jet 7 Pseudorapidity",50,-5.,5.);
  _hEtaObj7->SetXAxisTitle("#eta^{Jet 7}");
  _hEtaObj7->SetYAxisTitle("Events");
  // Jet 7 Phi
  _hPhiObj7 = DeclareTH1F("Jet7Phi","Jet 7 phi",50,-1.*TMath::Pi(),TMath::Pi());
  _hPhiObj7->SetXAxisTitle("#phi^{Jet 7}");
  _hPhiObj7->SetYAxisTitle("Events");

  // any Jet Pt
  _hPt = DeclareTH1F("JetPt","All Jet p_{T}",100,0.,400.);
  _hPt->SetXAxisTitle("p_{T}^{Jet} [GeV]");
  _hPt->SetYAxisTitle("Events");
  // any Jet Eta
  _hEta = DeclareTH1F("JetEta","All Jet Pseudorapidity",50,-5.,5.);
  _hEta->SetXAxisTitle("#eta^{Jet}");
  _hEta->SetYAxisTitle("Events");
  // any Jet Phi
  _hPhi = DeclareTH1F("JetPhi","All Jet phi",50,-1.*TMath::Pi(),TMath::Pi());
  _hPhi->SetXAxisTitle("#phi^{Jet}");
  _hPhi->SetYAxisTitle("Events");
  // any Jet Charge
  _hCharge = DeclareTH1F("JetCharge","All jet charge",3,-1.5,1.5);
  _hCharge->SetXAxisTitle("Charge^{jet}");
  _hCharge->SetYAxisTitle("Events");
  // any Jet E over P
  _hEtEM0 = DeclareTH1F("JetEtEM0","All Jet E_{T} in EM layer 0",100,0.,100.);
  _hEtEM0->SetXAxisTitle("E_{T}^{jet} [GeV]");
  _hEtEM0->SetYAxisTitle("Events");
  
  //cout<<"end of HistogrammingTaggedJet::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingTaggedJet::Apply()                                     *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingTaggedJet::Apply()
{
  //cout<<"Begin of HistogrammingTaggedJet::Apply()"<<endl;
  
  // Get Event Tree
  EventContainer *evc = GetEventContainer();
  
  // Fill Histograms
  _hNObj -> Fill(evc->taggedJets.size());

  // Jet 1
  if(evc->taggedJets.size()>0) {
    _hPtObj1  -> Fill(evc->taggedJets[0].Pt());
    _hEtaObj1 -> Fill(evc->taggedJets[0].Eta());
    _hPhiObj1 -> Fill(evc->taggedJets[0].Phi());
  } //if

  // Jet 2
  if(evc->taggedJets.size()>1) {
    _hPtObj2  -> Fill(evc->taggedJets[1].Pt());
    _hEtaObj2 -> Fill(evc->taggedJets[1].Eta());
    _hPhiObj2 -> Fill(evc->taggedJets[1].Phi());
  } //if
  
  // Jet 3
  if(evc->taggedJets.size()>2) {
    _hPtObj3  -> Fill(evc->taggedJets[2].Pt());
    _hEtaObj3 -> Fill(evc->taggedJets[2].Eta());
    _hPhiObj3 -> Fill(evc->taggedJets[2].Phi());
  } //if

  // Jet 4
  if(evc->taggedJets.size()>3) {
    _hPtObj4  -> Fill(evc->taggedJets[3].Pt());
    _hEtaObj4 -> Fill(evc->taggedJets[3].Eta());
    _hPhiObj4 -> Fill(evc->taggedJets[3].Phi());
  } //if

  // Jet 5
  if(evc->taggedJets.size()>4) {
    _hPtObj5  -> Fill(evc->taggedJets[4].Pt());
    _hEtaObj5 -> Fill(evc->taggedJets[4].Eta());
    _hPhiObj5 -> Fill(evc->taggedJets[4].Phi());
  } //if

  // Jet 6
  if(evc->taggedJets.size()>5) {
    _hPtObj6  -> Fill(evc->taggedJets[5].Pt());
    _hEtaObj6 -> Fill(evc->taggedJets[5].Eta());
    _hPhiObj6 -> Fill(evc->taggedJets[5].Phi());
  } //if

  // Jet 7
  if(evc->taggedJets.size()>6) {
    _hPtObj7  -> Fill(evc->taggedJets[6].Pt());
    _hEtaObj7 -> Fill(evc->taggedJets[6].Eta());
    _hPhiObj7 -> Fill(evc->taggedJets[6].Phi());
  } //if

  // Loop over all jets and histogram their properties
  for(int io = 0; io < evc -> taggedJets.size(); io++) {
    _hPt     -> Fill(evc -> taggedJets[io].Pt());
    _hEta    -> Fill(evc -> taggedJets[io].Eta());
    _hPhi    -> Fill(evc -> taggedJets[io].Phi());
    _hCharge -> Fill(evc -> taggedJets[io].GetCharge());
    //_hEtEM0  -> Fill(evc -> taggedJets[io].GetEtEM0());
  } //for 

  //cout<<"End of HistogrammingTaggedJet::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












