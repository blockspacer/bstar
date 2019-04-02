/******************************************************************************
 * HistogrammingZmumu.hpp                                                     *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Must pass a paramter describing the type of muon you want to plot          *
 *   Choices are All, Veto, Tight                                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingZmumu()           -- Parameterized Constructor             *
 *    ~HistogrammingZmumu()          -- Destructor                            *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingZmumu"       *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of muon properties                                   *
 *                                                                            *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      18 Jan 2007 - P. Ryan - modified for all/tight/veto/                  *
 *      06 Feb 2008 - P. Ryan - add 2d histograms                             *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingZmumu.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingZmumu::HistogrammingZmumu(EventContainer *obj)          *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/

HistogrammingZmumu::HistogrammingZmumu(EventContainer *obj, TString muonTypePassed, bool doPileupWeight)
{
  // Check muonType parameter
  if( muonTypePassed.CompareTo("All") && muonTypePassed.CompareTo("Isolated") && muonTypePassed.CompareTo("UnIsolated") &&
      muonTypePassed.CompareTo("Tight") && muonTypePassed.CompareTo("Veto") ){
    std::cout << "ERROR " << "<HistogrammingZmumu::HistogrammingZmumu()> " 
              << "Must pass Tight, Veto, or UnIsolated to constructor" << std::endl;
    exit(8);
  } //if
  muonType = muonTypePassed; 
  _doPileupWeight = doPileupWeight;

  SetEventContainer(obj);
} //HistogrammingZmumu()

/******************************************************************************
 * HistogrammingZmumu::~HistogrammingZmumu()                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingZmumu::~HistogrammingZmumu()
{
  
} //HistogrammingZmumu

/******************************************************************************
 * void HistogrammingZmumu::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/double HistogrammingZmumu::xy_to_phi(double x, double y){
   if(y>0) return  acos(x/sqrt(x*x+y*y));
   if(y<0) return -acos(x/sqrt(x*x+y*y));
}
double HistogrammingZmumu::dPhi(double phi1,double phi2){
        double dphi = fabs(phi1 - phi2);
        if(dphi>TMath::Pi()) dphi = dphi - 2*TMath::Pi();
        if(dphi<-TMath::Pi()) dphi = dphi + 2*TMath::Pi();
	return dphi;
}
void HistogrammingZmumu::BookHistogram(){
  
  Int_t ptBin    = 25;
  Int_t ptBin2D  = 10;
  Double_t ptMin = 0.0;
  Double_t ptMax = 100.0;

  Int_t etaBin    = 24;
  Int_t etaBin2D  = 10;
  Double_t etaMin = -3.0;
  Double_t etaMax = 3.0;

  Int_t phiBin    = 16;
  Int_t phiBin2D  = 10;
  Double_t phiMin = -1.0 * TMath::Pi();
  Double_t phiMax = TMath::Pi();

  ////////////////////////////////////////////////////////////////
  // Number of Objects
  ////////////////////////////////////////////////////////////////
  _hNObj = DeclareTH1F("MuonN","Number of Muons",10,0.,10.);
  _hNObj->SetXAxisTitle("N^{#mu}");
  _hNObj->SetYAxisTitle("Events");

  ////////////////////////////////////////////////////////////////
  // Kinematic Quantities - Individual
  ////////////////////////////////////////////////////////////////
  // 1 Pt
  _hPtObj1 = DeclareTH1F("Muon1Pt", "Muon 1 P_{T}", ptBin, ptMin, ptMax);
  _hPtObj1->SetXAxisTitle("P_{T} (Muon 1) [GeV]");
  _hPtObj1->SetYAxisTitle("Events");
  // 1 Pt - Blown up
  _hPtBlowupObj1 = DeclareTH1F("Muon1PtBlowup","Muon 1 P_{T}",200,0.,100.);
  _hPtBlowupObj1->SetXAxisTitle("p_{T} (Muon 1) [GeV]");
  _hPtBlowupObj1->SetYAxisTitle("Events");
  // 1 Eta
  _hEtaObj1 = DeclareTH1F("Muon1Eta", "Muon 1 Pseudorapidity", etaBin, etaMin, etaMax);
  _hEtaObj1->SetXAxisTitle("#eta (Muon 1) ");
  _hEtaObj1->SetYAxisTitle("Events");
  // 1 Phi
  _hPhiObj1 = DeclareTH1F("Muon1Phi", "Muon 1 phi", phiBin, phiMin, phiMax);
  _hPhiObj1->SetXAxisTitle("#phi (Muon 1)");
  _hPhiObj1->SetYAxisTitle("Events");
  // 1 Pt-Eta
  _hPtEtaObj1 = DeclareTH2F("Muon1PtEta", "Muon 1 pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEtaObj1->SetXAxisTitle("P_{T} [GeV]");
  _hPtEtaObj1->SetYAxisTitle("#eta");
  // 1 Pt-Phi
  _hPtPhiObj1 = DeclareTH2F("Muon1PtPhi", "Muon 1 pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhiObj1->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhiObj1->SetYAxisTitle("#phi");
  // 1 Eta-Phi
  _hEtaPhiObj1 = DeclareTH2F("Muon1EtaPhi", "Muon 1 eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhiObj1->SetXAxisTitle("#eta [GeV]");
  _hEtaPhiObj1->SetYAxisTitle("#phi");

  // 2 Pt
  _hPtObj2 = DeclareTH1F("Muon2Pt", "Muon 2 P_{T}", ptBin, ptMin, ptMax);
  _hPtObj2->SetXAxisTitle("P_{T} (Muon 2) [GeV]");
  _hPtObj2->SetYAxisTitle("Events");
  // 2 Eta
  _hEtaObj2 = DeclareTH1F("Muon2Eta", "Muon 2 Pseudorapidity", etaBin, etaMin, etaMax);
  _hEtaObj2->SetXAxisTitle("#eta (Muon 2) ");
  _hEtaObj2->SetYAxisTitle("Events");
  // 2 Phi
  _hPhiObj2 = DeclareTH1F("Muon2Phi", "Muon 2 phi", phiBin, phiMin, phiMax);
  _hPhiObj2->SetXAxisTitle("#phi (Muon 2)");
  _hPhiObj2->SetYAxisTitle("Events");
  // 2 Pt-Eta
  _hPtEtaObj2 = DeclareTH2F("Muon2PtEta", "Muon 2 pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEtaObj2->SetXAxisTitle("P_{T} [GeV]");
  _hPtEtaObj2->SetYAxisTitle("#eta");
  // 2 Pt-Phi
  _hPtPhiObj2 = DeclareTH2F("Muon2PtPhi", "Muon 2 pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhiObj2->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhiObj2->SetYAxisTitle("#phi");
  // 2 Eta-Phi
  _hEtaPhiObj2 = DeclareTH2F("Muon2EtaPhi", "Muon 2 eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhiObj2->SetXAxisTitle("#eta [GeV]");
  _hEtaPhiObj2->SetYAxisTitle("#phi");


// All muon jet Phi
  _hJetDeltaPhi = DeclareTH1F("MuonJetPhi", "Leading Muon Jet Delta phi", 160, phiMin, phiMax);
  _hJetDeltaPhi->SetXAxisTitle("#phi (Muon)");
  _hJetDeltaPhi->SetYAxisTitle("Events");

// All muon jet deltar
  _hJetDeltaR = DeclareTH1F("MuonJetR", "Leading Muon Jet Delta R", 160, 0, 6);
  _hJetDeltaR->SetXAxisTitle("#deltaR (Muon)");
  _hJetDeltaR->SetYAxisTitle("Events");


  ////////////////////////////////////////////////////////////////
  // Kinematic Qunatities - All
  // All muon Pt
  ////////////////////////////////////////////////////////////////
  _hPt = DeclareTH1F("MuonPt", "Any Muon P_{T}", ptBin, ptMin, ptMax);
  _hPt->SetXAxisTitle("P_{T} (Muon) [GeV]");
  _hPt->SetYAxisTitle("Events");
  // All muon Eta
  _hEta = DeclareTH1F("MuonEta", "Any Muon Pseudorapidity", etaBin, etaMin, etaMax);
  _hEta->SetXAxisTitle("#eta (Muon) ");
  _hEta->SetYAxisTitle("Events");
  // All muon Phi
  _hPhi = DeclareTH1F("MuonPhi", "Any Muon phi", phiBin, phiMin, phiMax);
  _hPhi->SetXAxisTitle("#phi (Muon)");
  _hPhi->SetYAxisTitle("Events");
  // All Pt-Eta
  _hPtEta = DeclareTH2F("MuonPtEta", "Any Muon pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEta->SetXAxisTitle("P_{T} [GeV]");
  _hPtEta->SetYAxisTitle("#eta");
  // All Pt-Phi
  _hPtPhi = DeclareTH2F("MuonPtPhi", "Any Muon pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhi->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhi->SetYAxisTitle("#phi");
  // All Eta-Phi
  _hEtaPhi = DeclareTH2F("MuonEtaPhi", "Any Muon eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhi->SetXAxisTitle("#eta");
  _hEtaPhi->SetYAxisTitle("#phi");

  // any muon Charge
  _hCharge = DeclareTH1F("MuonCharge","Any Muon charge",3,-1.5,1.5);
  _hCharge->SetXAxisTitle("Charge^{Muon}");
  _hCharge->SetYAxisTitle("Events");
  // any muon E over P
  _hEoverP = DeclareTH1F("MuonEoverP","Any Muon Energy over Momentum",50,0.,5.);
  _hEoverP->SetXAxisTitle("E/P^{Muon}");
  _hEoverP->SetYAxisTitle("Events");
  // any Zmumu invariant mass
  _hMZmumu = DeclareTH1F("MZ#mu#mu","Any Z#mu#mu invariant mass ",208,-8.,200.);
  _hMZmumu->SetXAxisTitle("M_{#mu#mu}");
  _hMZmumu->SetYAxisTitle("Events");
  // any Zmumu opposite sign invariant mass
  _hMZmumu_OS = DeclareTH1F("MZ#mu#mu_OS","Any Z#mu#mu Opposite Sign invariant mass ",208,-8.,200.);
  _hMZmumu_OS->SetXAxisTitle("M_{#mu#mu}^{OS}");
  _hMZmumu_OS->SetYAxisTitle("Events");
  // any Zmumu opposite sign invariant mass
  _hMZmumu_0jet_OS = DeclareTH1F("MZmumui_0jet_OS","Any Z#mu#mu Opposite Sign invariant mass in zero jet bin",208,-8.,200.);
  _hMZmumu_0jet_OS->SetXAxisTitle("M_{#mu#mu}^{OS}");
  _hMZmumu_0jet_OS->SetYAxisTitle("Events");
  // any Zmumu opposite sign invariant mass
  _hMZmumu_1jet_OS = DeclareTH1F("MZmumui_1jet_OS","Any Z#mu#mu Opposite Sign invariant mass in one jet bin",208,-8.,200.);
  _hMZmumu_1jet_OS->SetXAxisTitle("M_{#mu#mu}^{OS}");
  _hMZmumu_1jet_OS->SetYAxisTitle("Events");
  // any Zmumu opposite sign invariant mass
  _hMZmumu_2jet_OS = DeclareTH1F("MZmumui_2jet_OS","Any Z#mu#mu Opposite Sign invariant mass in two jets bin",208,-8.,200.);
  _hMZmumu_2jet_OS->SetXAxisTitle("M_{#mu#mu}^{OS}");
  _hMZmumu_2jet_OS->SetYAxisTitle("Events");
  // any Zmumu opposite sign invariant mass
  _hMZmumu_3jet_OS = DeclareTH1F("MZmumui_3jet_OS","Any Z#mu#mu Opposite Sign invariant mass in three jets bin",208,-8.,200.);
  _hMZmumu_3jet_OS->SetXAxisTitle("M_{#mu#mu}^{OS}");
  _hMZmumu_3jet_OS->SetYAxisTitle("Events");
  // any Zmumu opposite sign invariant mass
  _hMZmumu_4jet_OS = DeclareTH1F("MZmumui_4jet_OS","Any Z#mu#mu Opposite Sign invariant mass in four jets bin",208,-8.,200.);
  _hMZmumu_4jet_OS->SetXAxisTitle("M_{#mu#mu}^{OS}");
  _hMZmumu_4jet_OS->SetYAxisTitle("Events");
  // any Zmumu opposite sign invariant mass
  _hMZmumu_1plusjet_OS = DeclareTH1F("MZmumui_1plusjet_OS","Any Z#mu#mu Opposite Sign invariant mass in greater or equal one jets bin",208,-8.,200.);
  _hMZmumu_1plusjet_OS->SetXAxisTitle("M_{#mu#mu}^{OS}");
  _hMZmumu_1plusjet_OS->SetYAxisTitle("Events");
  // any Zmumu opposite sign invariant mass
  _hMZmumu_2plusjet_OS = DeclareTH1F("MZmumui_2plusjet_OS","Any Z#mu#mu Opposite Sign invariant mass in greater or equal two jets bin",208,-8.,200.);
  _hMZmumu_2plusjet_OS->SetXAxisTitle("M_{#mu#mu}^{OS}");
  _hMZmumu_2plusjet_OS->SetYAxisTitle("Events");
  // any Zmumu opposite sign invariant mass
  _hMZmumu_3plusjet_OS = DeclareTH1F("MZmumui_3plusjet_OS","Any Z#mu#mu Opposite Sign invariant mass in greater or equal three jets bin",208,-8.,200.);
  _hMZmumu_3plusjet_OS->SetXAxisTitle("M_{#mu#mu}^{OS}");
  _hMZmumu_3plusjet_OS->SetYAxisTitle("Events");
  // any Zmumu same sign invariant mass
  _hMZmumu_SS = DeclareTH1F("MZ#mu#mu_SS","Any Z#mu#mu same Sign invariant mass ",208,-8.,200.);
  _hMZmumu_SS->SetXAxisTitle("M_{#mu#mu}^{SS}");
  _hMZmumu_SS->SetYAxisTitle("Events");
  // missing Phi
  _hmissingPhi = DeclareTH1F("missingPhi"," missing phi ",80,-4.,4.);
  _hmissingPhi->SetXAxisTitle("#phi^{miss}");
  _hmissingPhi->SetYAxisTitle("Events");

  // Histogram of MET
  _hMET = DeclareTH1F("MET","Missing transverse energy",200,0.,200.);
  _hMET->SetXAxisTitle("MET [GeV]");
  _hMET->SetYAxisTitle("Events");
  // Mmumu vs MET
  _hMmumuMET = DeclareTH2F("M#mu#muMET", "M#mu#mu vs MET", 208, -8., 200, 200, 0., 200.);
  _hMmumuMET->SetXAxisTitle("M#mu#mu [GeV]");
  _hMmumuMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_OSMET = DeclareTH2F("M#mu#mu_OSMET", "M#mu#mu_OS vs MET", 208, -8., 200, 200, 0., 200.);
  _hMmumu_OSMET->SetXAxisTitle("M_{#mu#mu}^{OS} [GeV]");
  _hMmumu_OSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_0jet_OSMET = DeclareTH2F("Mmumu_0jet_OSMET", "M#mu#mu_OS vs MET zero jet bin", 208, -8., 200, 200, 0., 200.);
  _hMmumu_0jet_OSMET->SetXAxisTitle("M_{#mu#mu}^{OS} [GeV]");
  _hMmumu_0jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_1jet_OSMET = DeclareTH2F("Mmumu_1jet_OSMET", "M#mu#mu_OS vs MET one jet bin", 208, -8., 200, 200, 0., 200.);
  _hMmumu_1jet_OSMET->SetXAxisTitle("M_{#mu#mu}^{OS} [GeV]");
  _hMmumu_1jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_2jet_OSMET = DeclareTH2F("Mmumu_2jet_OSMET", "M#mu#mu_OS vs MET two jets bin", 208, -8., 200, 200, 0., 200.);
  _hMmumu_2jet_OSMET->SetXAxisTitle("M_{#mu#mu}^{OS} [GeV]");
  _hMmumu_2jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_3jet_OSMET = DeclareTH2F("Mmumu_3jet_OSMET", "M#mu#mu_OS vs MET three jets bin", 208, -8., 200, 200, 0., 200.);
  _hMmumu_3jet_OSMET->SetXAxisTitle("M_{#mu#mu}^{OS} [GeV]");
  _hMmumu_3jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_4jet_OSMET = DeclareTH2F("Mmumu_4jet_OSMET", "M#mu#mu_OS vs MET four jets bin", 208, -8., 200, 200, 0., 200.);
  _hMmumu_4jet_OSMET->SetXAxisTitle("M_{#mu#mu}^{OS} [GeV]");
  _hMmumu_4jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_1plusjet_OSMET = DeclareTH2F("Mmumu_1plusjet_OSMET", "M#mu#mu_OS vs MET greater equal one jets bin", 208, -8., 200, 200, 0., 200.);
  _hMmumu_1plusjet_OSMET->SetXAxisTitle("M_{#mu#mu}^{OS} [GeV]");
  _hMmumu_1plusjet_OSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_2plusjet_OSMET = DeclareTH2F("Mmumu_2plusjet_OSMET", "M#mu#mu_OS vs MET greater equal two jets bin", 208, -8., 200, 200, 0., 200.);
  _hMmumu_2plusjet_OSMET->SetXAxisTitle("M_{#mu#mu}^{OS} [GeV]");
  _hMmumu_2plusjet_OSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_3plusjet_OSMET = DeclareTH2F("Mmumu_3plusjet_OSMET", "M#mu#mu_OS vs MET greater equal three jets bin", 208, -8., 200, 200, 0., 200.);
  _hMmumu_3plusjet_OSMET->SetXAxisTitle("M_{#mu#mu}^{OS} [GeV]");
  _hMmumu_3plusjet_OSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu_OS vs MET
  _hMmumu_SSMET = DeclareTH2F("M#mu#mu_SSMET", "M#mu#mu_SS vs MET", 208, -8., 200, 200, 0., 200.);
  _hMmumu_SSMET->SetXAxisTitle("M_{#mu#mu}^{SS} [GeV]");
  _hMmumu_SSMET->SetYAxisTitle("MET [GeV]");
  // Mmumu vs dphi MET-Leading Jet
  _hMmumudPhiMETJ1 = DeclareTH2F("M#mu#mudphiMET_Jet", "M#mu#mu vs dPhi MET-Jet", 208, -8., 200, 32, 0, TMath::Pi());
  _hMmumudPhiMETJ1->SetXAxisTitle("M#mu#mu [GeV]");
  _hMmumudPhiMETJ1->SetYAxisTitle("#Delta #Phi (MET,Jet)");
  // Mmumu vs dphi MET-Leading Muon
  _hMmumudPhiMETE1 = DeclareTH2F("M#mu#mudphiMET_Muon", "M#mu#mu vs dPhi MET-Muon", 208, -8., 200, 32, 0, TMath::Pi());
  _hMmumudPhiMETE1->SetXAxisTitle("M#mu#mu [GeV]");
  _hMmumudPhiMETE1->SetYAxisTitle("#Delta #Phi (MET,Muon)");
  // Mmumu vs dphi MET-Z
  _hMmumudPhiMETZt = DeclareTH2F("M#mu#mudphiMET_Z", "M#mu#mu vs dPhi MET-Z", 208, -8., 200, 32, 0, TMath::Pi());
  _hMmumudPhiMETZt->SetXAxisTitle("M#mu#mu [GeV]");
  _hMmumudPhiMETZt->SetYAxisTitle("#Delta #Phi (MET,ZPt)");
  
  //cout<<"end of HistogrammingZmumu::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingZmumu::Apply()                                         *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingZmumu::Apply()
{
  //cout<<"Begin of HistogrammingZmumu::Apply()"<<endl;
  
  // Get Event Tree
  EventContainer *evc = GetEventContainer();
  
  // Muons class depends on 
  std::vector<Muon> muonVector;
  if(     "All"        == muonType) muonVector.assign(evc -> muons.begin(),           evc -> muons.end());
  else if("Veto"       == muonType) muonVector.assign(evc -> vetoMuons.begin(),       evc -> vetoMuons.end());
  else if("Tight"      == muonType) muonVector.assign(evc -> tightMuons.begin(),      evc -> tightMuons.end());
  else if("Isolated"   == muonType) muonVector.assign(evc -> isolatedMuons.begin(),   evc -> isolatedMuons.end());
  else if("UnIsolated" == muonType) muonVector.assign(evc -> unIsolatedMuons.begin(), evc -> unIsolatedMuons.end());
  else{
    std::cout << "ERROR " << "<HistogramminZmumu::Apply()> " 
	      << "muonType must be All, Tight, Isolated, UnIsolated, or Veto" << std::endl;
    exit(8);
  } //else

  double  _PileupWeight = 1;
  if(_doPileupWeight) _PileupWeight = fabs(evc->GetEventPileupWeight()) * evc -> GetEventWeight();
  else _PileupWeight = evc -> GetEventWeight();
  // Fill Histograms
  _hNObj -> Fill(muonVector.size(), _PileupWeight);


  if(muonVector.size()>0) {

    // TMP
    //if( evc -> L1Zmumus[0].GetDeltaRZmumu() < 0.2 )   // DeltaR Recon
    //if( evc -> L1Zmumus[0].GetDeltaRMCZmumu() < 0.2 )  // DeltaR Truth
    //if( evc -> L1Zmumus[0].GetDeltaRL2Zmumu() < 0.2 )  // DeltaR L2
    //if( evc -> L1Zmumus[0].GetDeltaREvtFilterZmumu() < 0.2 )  // DeltaR EF
    
    //if( evc -> L2Zmumus[0].GetDeltaRZmumu() < 0.2 )   // DeltaR Recon
    //if( evc -> L2Zmumus[0].GetDeltaRMCZmumu() < 0.2 )  // DeltaR Truth
    //if( evc -> L2Zmumus[0].GetDeltaRL1Zmumu() < 0.2 )  // DeltaR L1
    //if( evc -> L2Zmumus[0].GetDeltaREvtFilterZmumu() < 0.2 )  // DeltaR EF
    
    //if( evc -> EFZmumus[0].GetDeltaRZmumu() < 0.2 )   // DeltaR Recon
    //if( evc -> EFZmumus[0].GetDeltaRMCZmumu() < 0.2 )  // DeltaR Truth
    //if( evc -> EFZmumus[0].GetDeltaRL2Zmumu() < 0.2 )  // DeltaR L2
    //if( evc -> EFZmumus[0].GetDeltaRL1Zmumu() < 0.2 )  // DeltaR L1
    //{


    _hPtObj1  -> Fill(muonVector[0].Pt() , _PileupWeight);
    _hEtaObj1 -> Fill(muonVector[0].Eta(), _PileupWeight);
    _hPhiObj1 -> Fill(muonVector[0].Phi(), _PileupWeight);
    _hPtBlowupObj1  -> Fill(muonVector[0].Pt(), _PileupWeight);
    
    _hPtEtaObj1  -> Fill(muonVector[0].Pt(),  muonVector[0].Eta(), _PileupWeight);
    _hPtPhiObj1  -> Fill(muonVector[0].Pt(),  muonVector[0].Phi(), _PileupWeight);
    _hEtaPhiObj1 -> Fill(muonVector[0].Eta(), muonVector[0].Phi(), _PileupWeight);
    
    if(evc->jets.size()>0){
      _hJetDeltaR-> Fill(fabs(evc->jets[0].DeltaR(muonVector[0])), _PileupWeight);
      _hJetDeltaPhi-> Fill(fabs(evc->jets[0].DeltaPhi(muonVector[0])), _PileupWeight);
    }

    //      } //if DeltaR
  } // if size < 0
  if(muonVector.size()>1) {
    _hPtObj2  -> Fill(muonVector[1].Pt(), _PileupWeight);
    _hEtaObj2 -> Fill(muonVector[1].Eta(), _PileupWeight);
    _hPhiObj2 -> Fill(muonVector[1].Phi(), _PileupWeight);

    _hPtEtaObj2  -> Fill(muonVector[1].Pt(),  muonVector[1].Eta(), _PileupWeight);
    _hPtPhiObj2  -> Fill(muonVector[1].Pt(),  muonVector[1].Phi(), _PileupWeight);
    _hEtaPhiObj2 -> Fill(muonVector[1].Eta(), muonVector[1].Phi(), _PileupWeight);
  } // if size < 1
// Fill MZmumu
  double missingPhi = xy_to_phi(evc->missingEx, evc->missingEy); 
  if(muonVector.size()>1) {
      for(int it1=0; it1<muonVector.size()-1; it1++){
          for(int it2=it1+1; it2<muonVector.size(); it2++){
              double mmumu = (muonVector[it1] + muonVector[it2]).M();
              _hMZmumu -> Fill( mmumu , _PileupWeight);
              _hMmumuMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
              if(muonVector[it1].GetCharge() * muonVector[it2].GetCharge() < 0 ) {
                  _hMZmumu_OS -> Fill( mmumu , _PileupWeight);
                  _hMmumu_OSMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
	          if(evc->jets.size()==0) {
                      _hMZmumu_0jet_OS -> Fill( mmumu , _PileupWeight);
                      _hMmumu_0jet_OSMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
	          } else if(evc->jets.size()==1) {
                      _hMZmumu_1jet_OS -> Fill( mmumu , _PileupWeight);
                      _hMmumu_1jet_OSMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
	          } else if(evc->jets.size()==2) {
                      _hMZmumu_2jet_OS -> Fill( mmumu , _PileupWeight);
                      _hMmumu_2jet_OSMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
	          } else if(evc->jets.size()==3) {
                      _hMZmumu_3jet_OS -> Fill( mmumu , _PileupWeight);
                      _hMmumu_3jet_OSMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
	          } else if(evc->jets.size()==4) {
                      _hMZmumu_4jet_OS -> Fill( mmumu , _PileupWeight);
                      _hMmumu_4jet_OSMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
                  }
	          if(evc->jets.size()>=1) {
                      _hMZmumu_1plusjet_OS -> Fill( mmumu , _PileupWeight);
                      _hMmumu_1plusjet_OSMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
                  }
	          if(evc->jets.size()>=2) {
                      _hMZmumu_2plusjet_OS -> Fill( mmumu , _PileupWeight);
                      _hMmumu_2plusjet_OSMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
                  }
	          if(evc->jets.size()>=3) {
                      _hMZmumu_3plusjet_OS -> Fill( mmumu , _PileupWeight);
                      _hMmumu_3plusjet_OSMET -> Fill( mmumu, evc->missingEt , _PileupWeight);
                  }
              } else {
                  _hMZmumu_SS -> Fill( mmumu, _PileupWeight );
                  _hMmumu_SSMET -> Fill( mmumu, evc->missingEt, _PileupWeight );
	      }
// Fill Mmumu vs Toplogical variable
	   double delta = -100;
           _hmissingPhi ->Fill ( missingPhi, _PileupWeight );
	   if(evc->jets.size()>0) {
	     delta = fabs(dPhi(missingPhi,evc->jets[0].Phi())); 
	     _hMmumudPhiMETJ1->Fill(mmumu, delta, _PileupWeight); 
	   }
	   delta = fabs(dPhi(missingPhi,muonVector[0].Phi()));
	   _hMmumudPhiMETE1->Fill(mmumu, delta, _PileupWeight);
	   delta =  fabs(dPhi(missingPhi,(muonVector[0] + muonVector[1]).Phi())); 
	   _hMmumudPhiMETZt->Fill(mmumu, delta, _PileupWeight); 
          }
      }
  } else {
      _hMZmumu    -> Fill ( -8. , _PileupWeight);
      _hMZmumu_OS -> Fill ( -8. , _PileupWeight);
      _hMZmumu_0jet_OS -> Fill ( -8. , _PileupWeight);
      _hMZmumu_1jet_OS -> Fill ( -8. , _PileupWeight);
      _hMZmumu_2jet_OS -> Fill ( -8. , _PileupWeight);
      _hMZmumu_3jet_OS -> Fill ( -8. , _PileupWeight);
      _hMZmumu_4jet_OS -> Fill ( -8. , _PileupWeight);
      _hMZmumu_1plusjet_OS -> Fill ( -8. , _PileupWeight);
      _hMZmumu_2plusjet_OS -> Fill ( -8. , _PileupWeight);
      _hMZmumu_3plusjet_OS -> Fill ( -8. , _PileupWeight);
      _hMZmumu_SS -> Fill ( -8. , _PileupWeight);
      _hMmumuMET  -> Fill ( -8, evc->missingEt , _PileupWeight);
      _hMmumu_OSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
      _hMmumu_0jet_OSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
      _hMmumu_1jet_OSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
      _hMmumu_2jet_OSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
      _hMmumu_3jet_OSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
      _hMmumu_4jet_OSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
      _hMmumu_1plusjet_OSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
      _hMmumu_2plusjet_OSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
      _hMmumu_3plusjet_OSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
      _hMmumu_SSMET  -> Fill( -8, evc->missingEt , _PileupWeight);
  }
  // now loop over all muons and histogram their properties
  for(int io=0;io<muonVector.size();io++) {
    _hPt  -> Fill(muonVector[io].Pt() , _PileupWeight);
    _hEta -> Fill(muonVector[io].Eta(), _PileupWeight);
    _hPhi -> Fill(muonVector[io].Phi(), _PileupWeight);
    
    _hPtEta  -> Fill(muonVector[io].Pt(),  muonVector[io].Eta(), _PileupWeight);
    _hPtPhi  -> Fill(muonVector[io].Pt(),  muonVector[io].Phi(), _PileupWeight);
    _hEtaPhi -> Fill(muonVector[io].Eta(), muonVector[io].Phi(), _PileupWeight);

    _hCharge->Fill(muonVector[io].GetCharge(), _PileupWeight);
    //    _hEoverP->Fill(muonVector[io].GetEoverP());
  }    

  _hMET -> Fill( evc->missingEt , _PileupWeight);
  //cout<<"End of HistogrammingZmumu::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












