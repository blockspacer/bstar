/******************************************************************************
 * HistogrammingZemu.hpp                                                  *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Must pass a paramter describing the type of electron you want to plot      *
 *   Choices are All, Veto, Tight                                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingZemu()           -- Parameterized Constructor          *
 *    ~HistogrammingZemu()          -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingZemu"    *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of electron properties                               *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      18 Jan 2007 - P. Ryan - modified for all/tight/veto/                  *
 *      06 Feb 2008 - P. Ryan - add 2d histograms                             *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingZemu.hpp"
#include<iostream>

using namespace std;

void HistogrammingZemu::soulveZtautau(TLorentzVector lep1, TLorentzVector lep2, double METX, double METY, TLorentzVector & nv1, TLorentzVector & nv2){
    double phi1 = lep1.Phi();
    double phi2 = lep2.Phi();
    //double ptnv1 =  (METX/cos(phi2)-METY/sin(phi2))/(cos(phi1)/cos(phi2)-sin(phi1)/sin(phi2));
    double ptnv1, ptnv2;
    if(sin(phi2-phi1)!=0){
      ptnv1 =  (METX*sin(phi2)-METY*cos(phi2))/sin(phi2-phi1);
      ptnv2 =  (METX*sin(phi1)-METY*cos(phi1))/sin(phi1-phi2);
    } else {
        cout << " WARNING: No solution for coliner aproximation " << endl; 
        ptnv1 = sqrt(METX*METX+METY*METY)/2.; 
        ptnv2 = sqrt(METX*METX+METY*METY)/2.;
    }
    nv1.SetPtEtaPhiM(ptnv1,lep1.Eta(),lep1.Phi(),0);
    nv2.SetPtEtaPhiM(ptnv2,lep2.Eta(),lep2.Phi(),0);
}

double HistogrammingZemu::xy_to_phi(double x, double y){
   if(y>0) return  acos(x/sqrt(x*x+y*y));
   if(y<0) return -acos(x/sqrt(x*x+y*y));
}
double HistogrammingZemu::dPhi(double phi1,double phi2){
        double dphi = (phi1 - phi2);
        if(dphi>TMath::Pi()) dphi = dphi - 2*TMath::Pi();
        if(dphi<-TMath::Pi()) dphi = dphi + 2*TMath::Pi();
	return dphi;
}

/******************************************************************************
 * HistogrammingZemu::HistogrammingZemu(EventContainer *obj)          *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingZemu::HistogrammingZemu(EventContainer *obj, TString electronTypePassed, TString muonTypePassed,bool doPileupWeight)
{
  // Check muonType parameter
  if( muonTypePassed.CompareTo("All") && muonTypePassed.CompareTo("Isolated") && muonTypePassed.CompareTo("UnIsolated") &&
      muonTypePassed.CompareTo("Tight") && muonTypePassed.CompareTo("Veto") ){
    std::cout << "ERROR " << "<HistogrammingZemu::HistogrammingZemu()> "
              << "Must pass Tight, Veto, or UnIsolated to constructor" << std::endl;
    exit(8);
  } //if
  muonType = muonTypePassed;
  // Check electronType parameter
  if( electronTypePassed.CompareTo("All") && electronTypePassed.CompareTo("Isolated") && electronTypePassed.CompareTo("UnIsolated") &&
      electronTypePassed.CompareTo("Tight") && electronTypePassed.CompareTo("Veto") ){
    std::cout << "ERROR " << "<HistogrammingZemu::HistogrammingZemu()> " 
              << "Must pass Tight, Veto, or UnIsolated to constructor" << std::endl;
    exit(8);
  } //if
  electronType = electronTypePassed;

  _doPileupWeight = doPileupWeight;
  SetEventContainer(obj);
} //HistogrammingZemu()

/******************************************************************************
 * HistogrammingZemu::~HistogrammingZemu()                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingZemu::~HistogrammingZemu()
{
  
} //HistogrammingZemu

/******************************************************************************
 * void HistogrammingZemu::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingZemu::BookHistogram(){
  

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
  _hNObj = DeclareTH1F("EleN","Number of Electrons",10,0.,10.);
  _hNObj->SetXAxisTitle("N^{e}");
  _hNObj->SetYAxisTitle("Events");

  ////////////////////////////////////////////////////////////////
  // Kinematic Quantities - Individual
  ////////////////////////////////////////////////////////////////

// All electron jet Phi
  _hJetDeltaPhi = DeclareTH1F("EleJetPhi", "Leading Electron Jet Delta phi", 160, phiMin, phiMax);
  _hJetDeltaPhi->SetXAxisTitle("#phi (Ele)");
  _hJetDeltaPhi->SetYAxisTitle("Events");

// All electron jet deltar
  _hJetDeltaR = DeclareTH1F("EleJetR", "Leading Electron Jet Delta R", 160, 0, 6);
  _hJetDeltaR->SetXAxisTitle("#deltaR (Ele)");
  _hJetDeltaR->SetYAxisTitle("Events");


  ////////////////////////////////////////////////////////////////
  // Kinematic Qunatities - All
  // All electron Pt
  ////////////////////////////////////////////////////////////////
  _hPt = DeclareTH1F("ElePt", "Electron P_{T}", ptBin, ptMin, ptMax);
  _hPt->SetXAxisTitle("P_{T} (Ele) [GeV]");
  _hPt->SetYAxisTitle("Events");
  // All electron Eta
  _hEta = DeclareTH1F("EleEta", "Electron Pseudorapidity", etaBin, etaMin, etaMax);
  _hEta->SetXAxisTitle("#eta (Ele) ");
  _hEta->SetYAxisTitle("Events");
  // All electron Phi
  _hPhi = DeclareTH1F("ElePhi", "Electron phi", phiBin, phiMin, phiMax);
  _hPhi->SetXAxisTitle("#phi (Ele)");
  _hPhi->SetYAxisTitle("Events");
  // All Pt-Eta
  _hPtEta = DeclareTH2F("ElePtEta", "Electron pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEta->SetXAxisTitle("P_{T} [GeV]");
  _hPtEta->SetYAxisTitle("#eta");
  // All Pt-Phi
  _hPtPhi = DeclareTH2F("ElePtPhi", "Electron pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhi->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhi->SetYAxisTitle("#phi");
  // All Eta-Phi
  _hEtaPhi = DeclareTH2F("EleEtaPhi", "Electron eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhi->SetXAxisTitle("#eta");
  _hEtaPhi->SetYAxisTitle("#phi");

  // electron Charge
  _hCharge = DeclareTH1F("EleCharge","Electron charge",3,-1.5,1.5);
  _hCharge->SetXAxisTitle("Charge^{Ele}");
  _hCharge->SetYAxisTitle("Events");
  // electron E over P
  _hEoverP = DeclareTH1F("EleEoverP","Electron Energy over Momentum",50,0.,5.);
  _hEoverP->SetXAxisTitle("E/P^{Ele}");
  _hEoverP->SetYAxisTitle("Events");


  _hPtu = DeclareTH1F("MuonPt", "Muon P_{T}", ptBin, ptMin, ptMax);
  _hPtu->SetXAxisTitle("P_{T} (Muon) [GeV]");
  _hPtu->SetYAxisTitle("Events");
  // All muon Eta
  _hEtau = DeclareTH1F("MuonEta", "Muon Pseudorapidity", etaBin, etaMin, etaMax);
  _hEtau->SetXAxisTitle("#eta (Muon) ");
  _hEtau->SetYAxisTitle("Events");
  // All muon Phi
  _hPhiu = DeclareTH1F("MuonPhi", "Muon phi", phiBin, phiMin, phiMax);
  _hPhiu->SetXAxisTitle("#phi (Muon)");
  _hPhiu->SetYAxisTitle("Events");
  // All Pt-Eta
  _hPtEtau = DeclareTH2F("MuonPtEta", "Muon pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEtau->SetXAxisTitle("P_{T} [GeV]");
  _hPtEtau->SetYAxisTitle("#eta");
  // All Pt-Phi
  _hPtPhiu = DeclareTH2F("MuonPtPhi", "Muon pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhiu->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhiu->SetYAxisTitle("#phi");
  // All Eta-Phi
  _hEtaPhiu = DeclareTH2F("MuonEtaPhi", "Muon eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhiu->SetXAxisTitle("#eta");
  _hEtaPhiu->SetYAxisTitle("#phi");

  // muon Charge
  _hChargeu = DeclareTH1F("MuonCharge","Muon charge",3,-1.5,1.5);
  _hChargeu->SetXAxisTitle("Charge^{Muon}");
  _hChargeu->SetYAxisTitle("Events");


  // any Zemu invariant mass
  _hMZemu = DeclareTH1F("MZemu","Any Zemu invariant mass ",208,-8.,200.);
  _hMZemu->SetXAxisTitle("M_{e#mu}");
  _hMZemu->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hM_METcoLinearemu_OS = DeclareTH1F("M_METcoLinearemu_OS","Any Zemu Opposite Sign invariant mass + solved co Linear nutrinos",208,-8.,200.);
  _hM_METcoLinearemu_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hM_METcoLinearemu_OS->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hMZemu_OS = DeclareTH1F("MZemu_OS","Any Zemu Opposite Sign invariant mass ",208,-8.,200.);
  _hMZemu_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hMZemu_OS->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hMZemu_0jet_OS = DeclareTH1F("MZemu_0jet_OS","Any Zemu Opposite Sign invariant mass in zero jet bin",208,-8.,200.);
  _hMZemu_0jet_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hMZemu_0jet_OS->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hMZemu_1jet_OS = DeclareTH1F("MZemu_1jet_OS","Any Zemu Opposite Sign invariant mass in one jet bin",208,-8.,200.);
  _hMZemu_1jet_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hMZemu_1jet_OS->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hMZemu_2jet_OS = DeclareTH1F("MZemu_2jet_OS","Any Zemu Opposite Sign invariant mass in two jets bin",208,-8.,200.);
  _hMZemu_2jet_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hMZemu_2jet_OS->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hMZemu_3jet_OS = DeclareTH1F("MZemu_3jet_OS","Any Zemu Opposite Sign invariant mass in three jets bin",208,-8.,200.);
  _hMZemu_3jet_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hMZemu_3jet_OS->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hMZemu_4jet_OS = DeclareTH1F("MZemu_4jet_OS","Any Zemu Opposite Sign invariant mass in four jets bin",208,-8.,200.);
  _hMZemu_4jet_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hMZemu_4jet_OS->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hMZemu_1plusjet_OS = DeclareTH1F("MZemu_1plusjet_OS","Any Zemu Opposite Sign invariant mass in greater equal one jets bin",208,-8.,200.);
  _hMZemu_1plusjet_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hMZemu_1plusjet_OS->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hMZemu_2plusjet_OS = DeclareTH1F("MZemu_2plusjet_OS","Any Zemu Opposite Sign invariant mass in greater equal two jets bin",208,-8.,200.);
  _hMZemu_2plusjet_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hMZemu_2plusjet_OS->SetYAxisTitle("Events");
  // any Zemu opposite sign invariant mass
  _hMZemu_3plusjet_OS = DeclareTH1F("MZemu_3plusjet_OS","Any Zemu Opposite Sign invariant mass in greater equal three jets bin",208,-8.,200.);
  _hMZemu_3plusjet_OS->SetXAxisTitle("M_{e#mu}^{OS}");
  _hMZemu_3plusjet_OS->SetYAxisTitle("Events");
  // any Zemu same sign invariant mass
  _hMZemu_SS = DeclareTH1F("MZemu_SS","Any Zemu same Sign invariant mass ",208,-8.,200.);
  _hMZemu_SS->SetXAxisTitle("M_{e#mu}^{SS}");
  _hMZemu_SS->SetYAxisTitle("Events");
  // missing Phi
  _hmissingPhi = DeclareTH1F("missingPhi"," missing phi ",80,-4.,4.);
  _hmissingPhi->SetXAxisTitle("#phi^{miss}");
  _hmissingPhi->SetYAxisTitle("Events");

  // Histogram of MET
  _hMET = DeclareTH1F("MET","Missing transverse energy",200,0.,200.);
  _hMET->SetXAxisTitle("MET [GeV]");
  _hMET->SetYAxisTitle("Events");
  // Memu vs MET
  _hMemuMET = DeclareTH2F("MemuMET", "Memu vs MET", 208, -8., 200, 200, 0., 200.);
  _hMemuMET->SetXAxisTitle("Memu [GeV]");
  _hMemuMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_OSMET = DeclareTH2F("Memu_OSMET", "Memu_OS vs MET", 208, -8., 200, 200, 0., 200.);
  _hMemu_OSMET->SetXAxisTitle("M_{e#mu}^{OS} [GeV]");
  _hMemu_OSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_0jet_OSMET = DeclareTH2F("Memu_0jet_OSMET", "Memu_OS vs MET in zero jet bin", 208, -8., 200, 200, 0., 200.);
  _hMemu_0jet_OSMET->SetXAxisTitle("M_{e#mu}^{OS} [GeV]");
  _hMemu_0jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_1jet_OSMET = DeclareTH2F("Memu_1jet_OSMET", "Memu_OS vs MET in one jet bin", 208, -8., 200, 200, 0., 200.);
  _hMemu_1jet_OSMET->SetXAxisTitle("M_{e#mu}^{OS} [GeV]");
  _hMemu_1jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_2jet_OSMET = DeclareTH2F("Memu_2jet_OSMET", "Memu_OS vs MET in two jets bin", 208, -8., 200, 200, 0., 200.);
  _hMemu_2jet_OSMET->SetXAxisTitle("M_{e#mu}^{OS} [GeV]");
  _hMemu_2jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_3jet_OSMET = DeclareTH2F("Memu_3jet_OSMET", "Memu_OS vs MET in three jets bin", 208, -8., 200, 200, 0., 200.);
  _hMemu_3jet_OSMET->SetXAxisTitle("M_{e#mu}^{OS} [GeV]");
  _hMemu_3jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_4jet_OSMET = DeclareTH2F("Memu_4jet_OSMET", "Memu_OS vs MET in four jets bin", 208, -8., 200, 200, 0., 200.);
  _hMemu_4jet_OSMET->SetXAxisTitle("M_{e#mu}^{OS} [GeV]");
  _hMemu_4jet_OSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_1plusjet_OSMET = DeclareTH2F("Memu_1plusjet_OSMET", "Memu_OS vs MET in greater equal one jets bin", 208, -8., 200, 200, 0., 200.);
  _hMemu_1plusjet_OSMET->SetXAxisTitle("M_{e#mu}^{OS} [GeV]");
  _hMemu_1plusjet_OSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_2plusjet_OSMET = DeclareTH2F("Memu_2plusjet_OSMET", "Memu_OS vs MET in greater equal three jets bin", 208, -8., 200, 200, 0., 200.);
  _hMemu_2plusjet_OSMET->SetXAxisTitle("M_{e#mu}^{OS} [GeV]");
  _hMemu_2plusjet_OSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_3plusjet_OSMET = DeclareTH2F("Memu_3plusjet_OSMET", "Memu_OS vs MET in greater equal three jets bin", 208, -8., 200, 200, 0., 200.);
  _hMemu_3plusjet_OSMET->SetXAxisTitle("M_{e#mu}^{OS} [GeV]");
  _hMemu_3plusjet_OSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hM_METcoLinearemu_MET_OS = DeclareTH2F("M_METcoLinearemu_MET_OS", "Memu_OS vs MET", 208, -8., 200, 200, 0., 200.);
  _hM_METcoLinearemu_MET_OS->SetXAxisTitle("M_{e#mu+nutrino}^{OS} [GeV]");
  _hM_METcoLinearemu_MET_OS->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs MET
  _hMemu_SSMET = DeclareTH2F("Memu_SSMET", "Memu_SS vs MET", 208, -8., 200, 200, 0., 200.);
  _hMemu_SSMET->SetXAxisTitle("M_{e#mu}^{SS} [GeV]");
  _hMemu_SSMET->SetYAxisTitle("MET [GeV]");
  // Memu_OS vs dphi MET-Leading Jet
  _hMemu_OSdPhiMETJ1 = DeclareTH2F("Memu_OSdphiMET_Jet", "Memu_OS vs dPhi MET-Jet", 208, -8., 200, 64, -TMath::Pi(), TMath::Pi());
  _hMemu_OSdPhiMETJ1->SetXAxisTitle("Memu_OS [GeV]");
  _hMemu_OSdPhiMETJ1->SetYAxisTitle("#Delta #Phi (MET,Jet)");
  // Memu_OS vs dphi MET-Leading Electron
  _hMemu_OSdPhiMETE1 = DeclareTH2F("Memu_OSdphiMET_Electron", "Memu_OS vs dPhi MET-Electron", 208, -8., 200, 64, -TMath::Pi(), TMath::Pi());
  _hMemu_OSdPhiMETE1->SetXAxisTitle("Memu_OS [GeV]");
  _hMemu_OSdPhiMETE1->SetYAxisTitle("#Delta #Phi (MET,Elec)");
  // Memu_OS vs dphi MET-Z
  _hMemu_OSdPhiMETZt = DeclareTH2F("Memu_OSdphiMET_Z", "Memu_OS vs dPhi MET-Z", 208, -8., 200, 64, -TMath::Pi(), TMath::Pi());
  _hMemu_OSdPhiMETZt->SetXAxisTitle("Memu_OS [GeV]");
  _hMemu_OSdPhiMETZt->SetYAxisTitle("#Delta #Phi (MET,ZPt)");
  // Memu_OS vs dphi Z - Jet
  _hMemu_OSdPhiZtJet = DeclareTH2F("Memu_OSdphiZ_Jet", "Memu_OS vs dPhi Z_Jet", 208, -8., 200, 64, -TMath::Pi(), TMath::Pi());
  _hMemu_OSdPhiZtJet->SetXAxisTitle("Memu_OS [GeV]");
  _hMemu_OSdPhiZtJet->SetYAxisTitle("#Delta #Phi (ZPt, Jet)");

  // Memu_OS vs dphi MET-Leading Jet
  _hMemu_OSMET_projJ1 = DeclareTH2F("Memu_OSMET_projJet", "Memu_OS vs MET project along Jet", 208, -8., 200, 200, 0., 200.);
  _hMemu_OSMET_projJ1->SetXAxisTitle("Memu_OS [GeV]");
  _hMemu_OSMET_projJ1->SetYAxisTitle("MET * #cos( #Phi (MET,Jet) )");
  // Memu_OS vs dphi MET-Leading Electron
  _hMemu_OSMET_projE1 = DeclareTH2F("Memu_OSMET_projElectron", "Memu_OS vs MET project along Electron", 208, -8., 200, 200, 0., 200.);
  _hMemu_OSMET_projE1->SetXAxisTitle("Memu_OS [GeV]");
  _hMemu_OSMET_projE1->SetYAxisTitle("MET * #cos( #Phi (MET,Elec) )");
  // Memu_OS vs dphi MET-Z
  _hMemu_OSMET_projZt = DeclareTH2F("Memu_OSMET_projZ", "Memu_OS vs MET project along Z", 208, -8., 200, 200, 0., 200.);
  _hMemu_OSMET_projZt->SetXAxisTitle("Memu_OS [GeV]");
  _hMemu_OSMET_projZt->SetYAxisTitle("MET * #cos( #Phi (MET,ZPt) )");

  // MET vs dphi MET-Leading Jet
  _hMETdPhiMETJ1 = DeclareTH2F("METdphiMET_Jet", "MET vs dPhi MET-Jet", 208, -8., 200, 64, -TMath::Pi(), TMath::Pi());
  _hMETdPhiMETJ1->SetXAxisTitle("MET [GeV]");
  _hMETdPhiMETJ1->SetYAxisTitle("#Delta #Phi (MET,Jet)");
  // MET vs dphi MET-Leading Electron
  _hMETdPhiMETEl = DeclareTH2F("METdphiMET_Electron", "MET vs dPhi MET-Electron", 208, -8., 200, 64, -TMath::Pi(), TMath::Pi());
  _hMETdPhiMETEl->SetXAxisTitle("MET [GeV]");
  _hMETdPhiMETEl->SetYAxisTitle("#Delta #Phi (MET,Elec)");
  // MET vs dphi MET-Z
  _hMETdPhiMETZt = DeclareTH2F("METdphiMET_Z", "MET vs dPhi MET-Z", 208, -8., 200, 64, -TMath::Pi(), TMath::Pi());
  _hMETdPhiMETZt->SetXAxisTitle("MET [GeV]");
  _hMETdPhiMETZt->SetYAxisTitle("#Delta #Phi (MET,ZPt)");
  // MET vs dphi Z - Jet
  _hMETdPhiZtJet = DeclareTH2F("METdphiZ_Jet", "MET vs dPhi Z_Jet", 208, -8., 200, 64, -TMath::Pi(), TMath::Pi());
  _hMETdPhiZtJet->SetXAxisTitle("MET [GeV]");
  _hMETdPhiZtJet->SetYAxisTitle("#Delta #Phi (ZPt, Jet)");
  
  //cout<<"end of HistogrammingZemu::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingZemu::Apply()                                         *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingZemu::Apply()
{
  //cout<<"Begin of HistogrammingZemu::Apply()"<<endl;
  
  // Get Event Tree
  EventContainer *evc = GetEventContainer();
  
  // Electrons class depends on 
  std::vector<Electron> electronVector;
  if(     "All"        == electronType) electronVector.assign(evc -> electrons.begin(),           evc -> electrons.end());
  else if("Veto"       == electronType) electronVector.assign(evc -> vetoElectrons.begin(),       evc -> vetoElectrons.end());
  else if("Tight"      == electronType) electronVector.assign(evc -> tightElectrons.begin(),      evc -> tightElectrons.end());
  else if("Isolated"   == electronType) electronVector.assign(evc -> isolatedElectrons.begin(),   evc -> isolatedElectrons.end());
  else if("UnIsolated" == electronType) electronVector.assign(evc -> unIsolatedElectrons.begin(), evc -> unIsolatedElectrons.end());
  else{
    std::cout << "ERROR " << "<HistogramminZemu::Apply()> " 
	      << "electronType must be All, Tight, Isolated, UnIsolated, or Veto" << std::endl;
    exit(8);
  } //else

  // muons class depends on 
  std::vector<Muon> muonVector;
  if(     "All"        == muonType) muonVector.assign(evc -> muons.begin(),           evc -> muons.end());
  else if("Veto"       == muonType) muonVector.assign(evc -> vetoMuons.begin(),       evc -> vetoMuons.end());
  else if("Tight"      == muonType) muonVector.assign(evc -> tightMuons.begin(),      evc -> tightMuons.end());
  else if("Isolated"   == muonType) muonVector.assign(evc -> isolatedMuons.begin(),   evc -> isolatedMuons.end());
  else if("UnIsolated" == muonType) muonVector.assign(evc -> unIsolatedMuons.begin(), evc -> unIsolatedMuons.end());
  else{
    std::cout << "ERROR " << "<HistogramminZemu::Apply()> " 
	      << "muonType must be All, Tight, Isolated, UnIsolated, or Veto" << std::endl;
    exit(8);
  } //else


  // Fill Histograms
  double  _PileupWeight = 1;
  if(_doPileupWeight) _PileupWeight = fabs(evc->GetEventPileupWeight()) * evc -> GetEventWeight();
  else _PileupWeight = evc -> GetEventWeight();
  _hNObj -> Fill(electronVector.size() + muonVector.size(),_PileupWeight);

  // now loop over all electrons and histogram their properties
  //for(int io=0;io<electronVector.size();io++) {
  if(electronVector.size() > 0)
    {
    _hPt  -> Fill(electronVector[0].Pt(),_PileupWeight);
    _hEta -> Fill(electronVector[0].Eta(),_PileupWeight);
    _hPhi -> Fill(electronVector[0].Phi(),_PileupWeight);
    
    _hPtEta  -> Fill(electronVector[0].Pt(),  electronVector[0].Eta(),_PileupWeight);
    _hPtPhi  -> Fill(electronVector[0].Pt(),  electronVector[0].Phi(),_PileupWeight);
    _hEtaPhi -> Fill(electronVector[0].Eta(), electronVector[0].Phi(),_PileupWeight);

    _hCharge->Fill(electronVector[0].GetCharge(),_PileupWeight);
   }    

  if(muonVector.size() > 0)
    {
    _hPtu  -> Fill(muonVector[0].Pt(),_PileupWeight);
    _hEtau -> Fill(muonVector[0].Eta(),_PileupWeight);
    _hPhiu -> Fill(muonVector[0].Phi(),_PileupWeight);
    
    _hPtEtau  -> Fill(muonVector[0].Pt(),  muonVector[0].Eta(),_PileupWeight);
    _hPtPhiu  -> Fill(muonVector[0].Pt(),  muonVector[0].Phi(),_PileupWeight);
    _hEtaPhiu -> Fill(muonVector[0].Eta(), muonVector[0].Phi(),_PileupWeight);

    _hChargeu->Fill(muonVector[0].GetCharge(),_PileupWeight);
   }    



// Fill MZemu
  double missingPhi = xy_to_phi(evc->missingEx, evc->missingEy); 
  if(electronVector.size()>0 && muonVector.size() > 0) {
      for(int it1=0; it1<electronVector.size(); it1++){
          for(int it2=0; it2<muonVector.size(); it2++){
	    TLorentzVector tlv_ele = electronVector[it1];
	    TLorentzVector tlv_muo = muonVector[it2];
              double Memu = (tlv_ele + tlv_muo).M();
              TLorentzVector nv1, nv2;
              double METx = evc->missingEx;
              double METy = evc->missingEy;
              soulveZtautau(tlv_ele,tlv_muo,METx,METy,nv1,nv2);
              double M_METcoLinearemu = (tlv_ele + tlv_muo + nv1 + nv2).M();
              _hMZemu -> Fill( Memu ,_PileupWeight);
              _hMemuMET -> Fill( Memu, evc->missingEt ,_PileupWeight);
              if(electronVector[it1].GetCharge() * muonVector[it2].GetCharge() < 0 ) {
                  _hM_METcoLinearemu_OS -> Fill( M_METcoLinearemu ,_PileupWeight);
                  _hM_METcoLinearemu_MET_OS -> Fill( M_METcoLinearemu, evc->missingEt ,_PileupWeight);
                  _hMZemu_OS -> Fill( Memu ,_PileupWeight);
                  _hMemu_OSMET -> Fill( Memu, evc->missingEt ,_PileupWeight);
	          if(evc->jets.size()==0) {
                      _hMZemu_0jet_OS -> Fill( Memu , _PileupWeight);
                      _hMemu_0jet_OSMET -> Fill( Memu, evc->missingEt , _PileupWeight);
	          } else if(evc->jets.size()==1) {
                      _hMZemu_1jet_OS -> Fill( Memu , _PileupWeight);
                      _hMemu_1jet_OSMET -> Fill( Memu, evc->missingEt , _PileupWeight);
	          } else if(evc->jets.size()==2) {
                      _hMZemu_2jet_OS -> Fill( Memu , _PileupWeight);
                      _hMemu_2jet_OSMET -> Fill( Memu, evc->missingEt , _PileupWeight);
	          } else if(evc->jets.size()==3) {
                      _hMZemu_3jet_OS -> Fill( Memu , _PileupWeight);
                      _hMemu_3jet_OSMET -> Fill( Memu, evc->missingEt , _PileupWeight);
	          } else if(evc->jets.size()==4) {
                      _hMZemu_4jet_OS -> Fill( Memu , _PileupWeight);
                      _hMemu_4jet_OSMET -> Fill( Memu, evc->missingEt , _PileupWeight);
                  }
	          if(evc->jets.size()>=1) {
                      _hMZemu_1plusjet_OS -> Fill( Memu , _PileupWeight);
                      _hMemu_1plusjet_OSMET -> Fill( Memu, evc->missingEt , _PileupWeight);
                  }
	          if(evc->jets.size()>=2) {
                      _hMZemu_2plusjet_OS -> Fill( Memu , _PileupWeight);
                      _hMemu_2plusjet_OSMET -> Fill( Memu, evc->missingEt , _PileupWeight);
                  }
	          if(evc->jets.size()>=3) {
                      _hMZemu_3plusjet_OS -> Fill( Memu , _PileupWeight);
                      _hMemu_3plusjet_OSMET -> Fill( Memu, evc->missingEt , _PileupWeight);
                  }
// Fill Memu vs Toplogical variable
	          double delta = -100;
                  _hmissingPhi ->Fill ( missingPhi ,_PileupWeight);
	          if(evc->jets.size()>0) {
	            delta = (dPhi(missingPhi,evc->jets[0].Phi())); 
	            _hMemu_OSdPhiMETJ1->Fill(Memu, delta,_PileupWeight); 
	            _hMemu_OSMET_projJ1->Fill(Memu, evc->missingEt * fabs(cos(delta)) ,_PileupWeight); 
	            _hMETdPhiMETJ1->Fill(evc->missingEt, delta,_PileupWeight); 
	            delta = (dPhi((TLorentzVector(electronVector[it1]) + TLorentzVector(muonVector[it2])).Phi(),evc->jets[0].Phi())); 
	            _hMemu_OSdPhiZtJet->Fill(Memu, delta,_PileupWeight); 
	            _hMETdPhiZtJet->Fill(evc->missingEt, delta,_PileupWeight); 
	          }
	          delta = (dPhi(missingPhi,electronVector[0].Phi()));
	          _hMemu_OSdPhiMETE1->Fill(Memu, delta,_PileupWeight);
	          _hMemu_OSMET_projE1->Fill(Memu, evc->missingEt * fabs(cos(delta)) ,_PileupWeight); 
	          _hMETdPhiMETEl->Fill(evc->missingEt, delta,_PileupWeight);
	          delta = (dPhi(missingPhi,(TLorentzVector(electronVector[it1]) + TLorentzVector(muonVector[it2])).Phi())); 
	          _hMemu_OSdPhiMETZt->Fill(Memu, delta,_PileupWeight); 
	          _hMemu_OSMET_projZt->Fill(Memu, evc->missingEt * fabs(cos(delta)) ,_PileupWeight); 
	          _hMETdPhiMETZt->Fill(evc->missingEt, delta,_PileupWeight); 
              } else {
                  _hMZemu_SS -> Fill( Memu ,_PileupWeight);
                  _hMemu_SSMET -> Fill( Memu, evc->missingEt ,_PileupWeight);
	      }
          }
      }
  } else {
      _hMZemu    -> Fill ( -8. ,_PileupWeight);
      _hMZemu_OS -> Fill ( -8. ,_PileupWeight);
      _hMZemu_0jet_OS -> Fill ( -8. ,_PileupWeight);
      _hMZemu_1jet_OS -> Fill ( -8. ,_PileupWeight);
      _hMZemu_2jet_OS -> Fill ( -8. ,_PileupWeight);
      _hMZemu_3jet_OS -> Fill ( -8. ,_PileupWeight);
      _hMZemu_4jet_OS -> Fill ( -8. ,_PileupWeight);
      _hMZemu_1plusjet_OS -> Fill ( -8. ,_PileupWeight);
      _hMZemu_2plusjet_OS -> Fill ( -8. ,_PileupWeight);
      _hMZemu_3plusjet_OS -> Fill ( -8. ,_PileupWeight);
      _hMZemu_SS -> Fill ( -8. ,_PileupWeight);
      _hMemuMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_OSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_0jet_OSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_1jet_OSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_2jet_OSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_3jet_OSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_4jet_OSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_1plusjet_OSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_2plusjet_OSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_3plusjet_OSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
      _hMemu_SSMET  -> Fill( -8, evc->missingEt ,_PileupWeight);
  }

  _hMET -> Fill( evc->missingEt ,_PileupWeight);
  //cout<<"End of HistogrammingZemu::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












