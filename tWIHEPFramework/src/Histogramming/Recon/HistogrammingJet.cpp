/******************************************************************************
 * HistogrammingJet.cpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingJet()              -- Parameterized Constructor            *
 *    ~HistogrammingJet()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingJet"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of jet properties                                    *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingJet.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingJet::HistogrammingJet(EventContainer *obj)                    *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingJet::HistogrammingJet(EventContainer *obj)
{
  SetEventContainer(obj);
} //HistogrammingJet()

/******************************************************************************
 * HistogrammingJet::~HistogrammingJet()                                      *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingJet::~HistogrammingJet()
{

} //HistogrammingJet

/******************************************************************************
 * void HistogrammingJet::BookHistogram()                                     *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingJet::BookHistogram(){
  
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
  // Jet 1 Tag
  _hTagObj1 = DeclareTH1F("Jet1Tag","Jet 1 tag",50,-1.,1.);
  _hTagObj1->SetXAxisTitle("Tag^{Jet 1}");
  _hTagObj1->SetYAxisTitle("Events");
  

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
  // Jet 2 Tag
  _hTagObj2 = DeclareTH1F("Jet2Tag","Jet 2 tag",50,-1.,1.);
  _hTagObj2->SetXAxisTitle("Tag^{Jet 2}");
  _hTagObj2->SetYAxisTitle("Events");

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
  // Jet 3 Tag
  _hTagObj3 = DeclareTH1F("Jet3Tag","Jet 3 tag",50,-1.,1.);
  _hTagObj3->SetXAxisTitle("Tag^{Jet 3}");
  _hTagObj3->SetYAxisTitle("Events");
  
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
  // Jet 4 Tag
  _hTagObj4 = DeclareTH1F("Jet4Tag","Jet 4 tag",50,-1.,1.);
  _hTagObj4->SetXAxisTitle("Tag^{Jet 4}");
  _hTagObj4->SetYAxisTitle("Events");

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
  // Jet 5 Tag
  _hTagObj5 = DeclareTH1F("Jet5Tag","Jet 5 tag",50,-1.,1.);
  _hTagObj5->SetXAxisTitle("Tag^{Jet 5}");
  _hTagObj5->SetYAxisTitle("Events");

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
  // Jet 6 Tag
  _hTagObj6 = DeclareTH1F("Jet6Tag","Jet 6 tag",50,-1.,1.);
  _hTagObj6->SetXAxisTitle("Tag^{Jet 6}");
  _hTagObj6->SetYAxisTitle("Events");

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
  // Jet 7 Tag
  _hTagObj7 = DeclareTH1F("Jet7Tag","Jet 7 tag",50,-1.,1.);
  _hTagObj7->SetXAxisTitle("Tag^{Jet 7}");
  _hTagObj7->SetYAxisTitle("Events");

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
  // any Jet E over P
  _hTruthLabel = DeclareTH1F("JetTruthLabel","All Jet Truth ID",100,-1050.,100.);
  _hTruthLabel->SetXAxisTitle("Truth Label");
  _hTruthLabel->SetYAxisTitle("Events");
  
  //cout<<"end of HistogrammingJet::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingJet::Apply()                                         *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingJet::Apply()
{
  //cout<<"Begin of HistogrammingJet::Apply()"<<endl;
  // Get Event Tree
  EventContainer *evc = GetEventContainer();

  // Fill Histograms
  _hNObj -> Fill(evc->jets.size());

  // Jet 1
  if(evc->jets.size()>0) {
    _hPtObj1  -> Fill(evc->jets[0].Pt());
    _hEtaObj1 -> Fill(evc->jets[0].Eta());
    _hPhiObj1 -> Fill(evc->jets[0].Phi());
    _hTagObj1 -> Fill(evc->jets[0].bDiscriminator());
  } //if

  // Jet 2
  if(evc->jets.size()>1) {
    _hPtObj2  -> Fill(evc->jets[1].Pt());
    _hEtaObj2 -> Fill(evc->jets[1].Eta());
    _hPhiObj2 -> Fill(evc->jets[1].Phi());
    _hTagObj2 -> Fill(evc->jets[1].bDiscriminator());
  } //if
  
  // Jet 3
  if(evc->jets.size()>2) {
    _hPtObj3  -> Fill(evc->jets[2].Pt());
    _hEtaObj3 -> Fill(evc->jets[2].Eta());
    _hPhiObj3 -> Fill(evc->jets[2].Phi());
    _hTagObj3 -> Fill(evc->jets[2].bDiscriminator());
  } //if

  // Jet 4
  if(evc->jets.size()>3) {
    _hPtObj4  -> Fill(evc->jets[3].Pt());
    _hEtaObj4 -> Fill(evc->jets[3].Eta());
    _hPhiObj4 -> Fill(evc->jets[3].Phi());
    _hTagObj4 -> Fill(evc->jets[3].bDiscriminator());
  } //if

  // Jet 5
  if(evc->jets.size()>4) {
    _hPtObj5  -> Fill(evc->jets[4].Pt());
    _hEtaObj5 -> Fill(evc->jets[4].Eta());
    _hPhiObj5 -> Fill(evc->jets[4].Phi());
    _hTagObj5 -> Fill(evc->jets[4].bDiscriminator());
  } //if

  // Jet 6
  if(evc->jets.size()>5) {
    _hPtObj6  -> Fill(evc->jets[5].Pt());
    _hEtaObj6 -> Fill(evc->jets[5].Eta());
    _hPhiObj6 -> Fill(evc->jets[5].Phi());
    _hTagObj6 -> Fill(evc->jets[5].bDiscriminator());
  } //if

  // Jet 7
  if(evc->jets.size()>6) {
    _hPtObj7  -> Fill(evc->jets[6].Pt());
    _hEtaObj7 -> Fill(evc->jets[6].Eta());
    _hPhiObj7 -> Fill(evc->jets[6].Phi());
    _hTagObj7 -> Fill(evc->jets[6].bDiscriminator());
  } //if

  // Loop over all jets and histogram their properties
  for(int io = 0; io < evc -> jets.size(); io++) {
    _hPt     -> Fill(evc -> jets[io].Pt());
    _hEta    -> Fill(evc -> jets[io].Eta());
    _hPhi    -> Fill(evc -> jets[io].Phi());
    _hCharge -> Fill(evc -> jets[io].GetCharge());
    //_hEtEM0  -> Fill(evc -> jets[io].GetEtEM0());
    //_hTruthLabel  -> Fill(evc -> jets[io].GetPdgId());
  } //for 

  //cout<<"End of HistogrammingJet::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












