/******************************************************************************
 * HistogrammingMuon.cpp                                                      *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMuon()              -- Parameterized Constructor           *
 *    ~HistogrammingMuon()             -- Destructor                          *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMuon"        *
 *                                                                            *
 * Private Data Members of this class                                         *
 * - lots of histograms of muon properties                                    *
 *                                                                            *
 * NOTE: Detector histograms commented out- uncommented if wish to use        *
 *                                                                            *
 * History                                                                    *
 *      17 June 2015 - Created by Huaqiao ZHANG                               *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingMuon.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingMuon::HistogrammingMuon(EventContainer *obj)                  *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingMuon::HistogrammingMuon(EventContainer *obj, TString muonTypePassed)
{
  // Check muonType parameter
  if( muonTypePassed.CompareTo("All")   && muonTypePassed.CompareTo("UnIsolated") && 
      muonTypePassed.CompareTo("Tight") && muonTypePassed.CompareTo("Veto") &&
      muonTypePassed.CompareTo("Isolated") && muonTypePassed.CompareTo("PtEtaCut") ){
    std::cout << "ERROR " << "<HistogrammingMuon::HistogrammingMuon()> " 
              << "Must pass Tight, Veto, Isolated, UnIsolated, PtEtaCut, or All to constructor" << std::endl;
    exit(8);
  } //if
  muonType = muonTypePassed;

  SetEventContainer(obj);
} //HistogrammingMuon()

/******************************************************************************
 * HistogrammingMuon::~HistogrammingMuon()                                    *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingMuon::~HistogrammingMuon()
{
  
} //HistogrammingMuon

/******************************************************************************
 * void HistogrammingMuon::BookHistogram()                                    *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingMuon::BookHistogram(){
  
  Int_t RPCEtaBin    = 15;
  Double_t RPCEtaMin = 0.0;
  Double_t RPCEtaMax = 15.0;

  Int_t RPCPhiBin    = 15;
  Double_t RPCPhiMin = 0.0;
  Double_t RPCPhiMax = 15.0;

  Int_t TGCEtaBin    = 15;
  Double_t TGCEtaMin = 0.0;
  Double_t TGCEtaMax = 15.0;

  Int_t TGCPhiBin    = 15;
  Double_t TGCPhiMin = 0.0;
  Double_t TGCPhiMax = 15.0;

  Int_t ptBin    = 25;
  Int_t ptBin2D  = 10;
  Double_t ptMin = 0.0;
  Double_t ptMax = 100.0;

  Int_t relIsoBin = 30;
  Double_t relIsoMin = 0.0;
  Double_t relIsoMax = 2.0;

  //  Int_t etaBin    = 60;
  Int_t etaBin    = 24;
  Int_t etaBin2D  = 10;
  Double_t etaMin = -3.0;
  Double_t etaMax = 3.0;

  //  Int_t phiBin    = 25;
  Int_t phiBin    = 16;
  Int_t phiBin2D  = 10;
  Double_t phiMin = -1.0 * TMath::Pi();
  Double_t phiMax = TMath::Pi();
    

  ////////////////////////////////////////////////////////////////
  // Number of Objects
  ////////////////////////////////////////////////////////////////
  _hNObj = DeclareTH1F("MuN","Number of Muons",10,0.,10.);
  _hNObj->SetXAxisTitle("N");
  _hNObj->SetYAxisTitle("Events");


  ////////////////////////////////////////////////////////////////
  // Kinematic Quantities - Individual
  ////////////////////////////////////////////////////////////////
  // 1 Pt
  _hPtObj1 = DeclareTH1F("Mu1Pt", "Muon 1 P_{T}", ptBin, ptMin, ptMax);
  _hPtObj1->SetXAxisTitle("P_{T} (muon 1) [GeV]");
  _hPtObj1->SetYAxisTitle("Events");
  // 1 Eta
  _hEtaObj1 = DeclareTH1F("Mu1Eta", "Muon 1 Pseudorapidity", etaBin, etaMin, etaMax);
  _hEtaObj1->SetXAxisTitle("#eta (muon 1) ");
  _hEtaObj1->SetYAxisTitle("Events");
  // 1 Phi
  _hPhiObj1 = DeclareTH1F("Mu1Phi", "Muon 1 phi", phiBin, phiMin, phiMax);
  _hPhiObj1->SetXAxisTitle("#phi (muon 1)");
  _hPhiObj1->SetYAxisTitle("Events");
  // 1 relIso
  _hRelIsoObj1 = DeclareTH1F("Mu1RelIso", "Muon 1 relIso", relIsoBin, relIsoMin, relIsoMax);
  _hRelIsoObj1->SetXAxisTitle("relIso (muon 1)");
  _hRelIsoObj1->SetYAxisTitle("Events");
  // 1 Pt-Eta
  _hPtEtaObj1 = DeclareTH2F("Mu1PtEta", "Muon 1 pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEtaObj1->SetXAxisTitle("P_{T} [GeV]");
  _hPtEtaObj1->SetYAxisTitle("#eta");
  // 1 Pt-Phi
  _hPtPhiObj1 = DeclareTH2F("Mu1PtPhi", "Muon 1 pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhiObj1->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhiObj1->SetYAxisTitle("#phi");
  // 1 Eta-Phi
  _hEtaPhiObj1 = DeclareTH2F("Mu1EtaPhi", "Muon 1 eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhiObj1->SetXAxisTitle("#eta [GeV]");
  _hEtaPhiObj1->SetYAxisTitle("#phi");

  // 2 Pt
  _hPtObj2 = DeclareTH1F("Mu2Pt", "Muon 2 P_{T}", ptBin, ptMin, ptMax);
  _hPtObj2->SetXAxisTitle("P_{T} (muon 2) [GeV]");
  _hPtObj2->SetYAxisTitle("Events");
  // 2 Eta
  _hEtaObj2 = DeclareTH1F("Mu2Eta", "Muon 2 Pseudorapidity", etaBin, etaMin, etaMax);
  _hEtaObj2->SetXAxisTitle("#eta (muon 2) ");
  _hEtaObj2->SetYAxisTitle("Events");
  // 2 Phi
  _hPhiObj2 = DeclareTH1F("Mu2Phi", "Muon 2 phi", phiBin, phiMin, phiMax);
  _hPhiObj2->SetXAxisTitle("#phi (muon 2)");
  _hPhiObj2->SetYAxisTitle("Events");
  // 2 relIso
  _hRelIsoObj2 = DeclareTH1F("Mu2RelIso", "Muon 2 relIso", relIsoBin, relIsoMin, relIsoMax);
  _hRelIsoObj2->SetXAxisTitle("relIso (muon 2)");
  _hRelIsoObj2->SetYAxisTitle("Events");
  // 2 Pt-Eta
  _hPtEtaObj2 = DeclareTH2F("Mu2PtEta", "Muon 2 pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEtaObj2->SetXAxisTitle("P_{T} [GeV]");
  _hPtEtaObj2->SetYAxisTitle("#eta");
  // 2 Pt-Phi
  _hPtPhiObj2 = DeclareTH2F("Mu2PtPhi", "Muon 2 pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhiObj2->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhiObj2->SetYAxisTitle("#phi");
  // 2 Eta-Phi
  _hEtaPhiObj2 = DeclareTH2F("Mu2EtaPhi", "Muon 2 eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhiObj2->SetXAxisTitle("#eta");
  _hEtaPhiObj2->SetYAxisTitle("#phi");

// All Muon jet Phi
  _hJetDeltaPhi = DeclareTH1F("MuJetPhi", "Leading Muon Jet Delta phi", 160, phiMin, phiMax);
  _hJetDeltaPhi->SetXAxisTitle("#phi (Mu)");
  _hJetDeltaPhi->SetYAxisTitle("Events");


// All Muon jet deltar
  _hJetDeltaR = DeclareTH1F("MuJetR", "Leading Muon Jet Delta R", 160, 0, 6);
  _hJetDeltaR->SetXAxisTitle("#deltaR (Mu)");
  _hJetDeltaR->SetYAxisTitle("Events");

  ////////////////////////////////////////////////////////////////
  // Kinematic Qunatities - All
  ////////////////////////////////////////////////////////////////
  // All Muon Pt
  _hPt = DeclareTH1F("MuPt", "Any Muon P_{T}", ptBin, ptMin, ptMax);
  _hPt->SetXAxisTitle("P_{T} (muon) [GeV]");
  _hPt->SetYAxisTitle("Events");
  // All Muon Eta
  _hEta = DeclareTH1F("MuEta", "Any Muon Pseudorapidity", etaBin, etaMin, etaMax);
  _hEta->SetXAxisTitle("#eta (muon) ");
  _hEta->SetYAxisTitle("Events");
  // All Muon Phi
  _hPhi = DeclareTH1F("MuPhi", "Any Muon phi", phiBin, phiMin, phiMax);
  _hPhi->SetXAxisTitle("#phi (muon)");
  _hPhi->SetYAxisTitle("Events");
  // All relIso
  _hRelIso = DeclareTH1F("MuRelIso", "Any Muon relIso", relIsoBin, relIsoMin, relIsoMax);
  _hRelIso->SetXAxisTitle("relIso (muon)");
  _hRelIso->SetYAxisTitle("Events");
  // All Pt-Eta
  _hPtEta = DeclareTH2F("MuPtEta", "Any Muon pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEta->SetXAxisTitle("P_{T} [GeV]");
  _hPtEta->SetYAxisTitle("#eta");
  // All Pt-Phi
  _hPtPhi = DeclareTH2F("MuPtPhi", "Any Muon pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhi->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhi->SetYAxisTitle("#phi");
  // All Eta-Phi
  _hEtaPhi = DeclareTH2F("MuEtaPhi", "Any Muon eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhi->SetXAxisTitle("#eta");
  _hEtaPhi->SetYAxisTitle("#phi");

  // All Muon Charge
  _hCharge = DeclareTH1F("MuCharge", "Any muon charge", 3, -1.5, 1.5);
  _hCharge->SetXAxisTitle("charge (muon)");
  _hCharge->SetYAxisTitle("Events");

  ////////////////////////////////////////////////////////////////
  //cout<<"end of HistogrammingMuon::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingMuon::Apply()                                          *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingMuon::Apply()
{
  //cout<<"Begin of HistogrammingMuon::Apply()"<<endl;

  // Get Event Tree
  EventContainer *evc = GetEventContainer();

  // Muons class depends on 
  std::vector<Muon> muonVector;
  if(     "All"        == muonType) muonVector.assign(evc -> muons.begin(),           evc -> muons.end());
  else if("Veto"       == muonType) muonVector.assign(evc -> vetoMuons.begin(),       evc -> vetoMuons.end());
  else if("PtEtaCut"   == muonType) muonVector.assign(evc -> ptetaMuons.begin(),      evc -> ptetaMuons.end());
  else if("Tight"      == muonType) muonVector.assign(evc -> tightMuons.begin(),      evc -> tightMuons.end());
  else if("Isolated"   == muonType) muonVector.assign(evc -> isolatedMuons.begin(),   evc -> isolatedMuons.end());
  else if("UnIsolated" == muonType) muonVector.assign(evc -> unIsolatedMuons.begin(), evc -> unIsolatedMuons.end());
  else{
    std::cout << "ERROR " << "<HistogramminMuon::Apply()> " 
	      << "muonType must be All, Tight, Veto, Isolated, or UnIsolated, PtEtaCut" << std::endl;
    exit(8);
  } //else
    
  // Fill Histograms
  _hNObj -> Fill(muonVector.size());

  if(muonVector.size()>0) {

    _hPtObj1  -> Fill(muonVector[0].Pt());
    _hEtaObj1 -> Fill(muonVector[0].Eta());
    _hPhiObj1 -> Fill(muonVector[0].Phi());
    _hRelIsoObj1 -> Fill(muonVector[0].relIsoR04());
    
    _hPtEtaObj1  -> Fill(muonVector[0].Pt(),  muonVector[0].Eta());
    _hPtPhiObj1  -> Fill(muonVector[0].Pt(),  muonVector[0].Phi());
    _hEtaPhiObj1 -> Fill(muonVector[0].Eta(), muonVector[0].Phi());

  } // if size < 0

  if(muonVector.size()>1) {
    _hPtObj2  -> Fill(muonVector[1].Pt());
    _hEtaObj2 -> Fill(muonVector[1].Eta());
    _hPhiObj2 -> Fill(muonVector[1].Phi());
    _hRelIsoObj2 -> Fill(muonVector[1].relIsoR04());

    _hPtEtaObj2  -> Fill(muonVector[1].Pt(),  muonVector[1].Eta());
    _hPtPhiObj2  -> Fill(muonVector[1].Pt(),  muonVector[1].Phi());
    _hEtaPhiObj2 -> Fill(muonVector[1].Eta(), muonVector[1].Phi());

  } // if size < 1

  // now loop over all vetoMuons and histogram their properties
  for(int io=0;io<muonVector.size();io++) {
    _hPt  -> Fill(muonVector[io].Pt());
    _hEta -> Fill(muonVector[io].Eta());
    _hPhi -> Fill(muonVector[io].Phi());
    _hRelIso -> Fill(muonVector[io].relIsoR04());

    _hPtEta  -> Fill(muonVector[io].Pt(),  muonVector[io].Eta());
    _hPtPhi  -> Fill(muonVector[io].Pt(),  muonVector[io].Phi());
    _hEtaPhi -> Fill(muonVector[io].Eta(), muonVector[io].Phi());

    _hCharge->Fill(muonVector[io].GetCharge());


    double minDphi = 9999.; double minDr = 9999.;
    for(int ijet=0; ijet<evc->jets.size(); ijet++){
	if(fabs(evc->jets[0].DeltaPhi(muonVector[0])) < minDphi) minDphi = fabs(evc->jets[0].DeltaPhi(muonVector[0]));
	if(fabs(evc->jets[0].DeltaR(muonVector[0])) < minDr) minDr = fabs(evc->jets[0].DeltaR(muonVector[0]));
    }
    _hJetDeltaPhi-> Fill(minDphi);
    _hJetDeltaR  -> Fill(minDr);

  } //for    

  
  //cout<<"End of HistogrammingMuon::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












