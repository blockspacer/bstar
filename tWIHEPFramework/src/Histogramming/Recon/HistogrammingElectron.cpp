/******************************************************************************
 * HistogrammingElectron.hpp                                                  *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Must pass a paramter describing the type of electron you want to plot      *
 *   Choices are All, Veto, Tight                                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingElectron()           -- Parameterized Constructor          *
 *    ~HistogrammingElectron()          -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingElectron"    *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of electron properties                               *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      18 Jan 2007 - P. Ryan - modified for all/tight/veto/                  *
 *      06 Feb 2008 - P. Ryan - add 2d histograms                             *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingElectron.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingElectron::HistogrammingElectron(EventContainer *obj)          *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingElectron::HistogrammingElectron(EventContainer *obj, TString electronTypePassed)
{
  // Check electronType parameter
  if( electronTypePassed.CompareTo("All") && electronTypePassed.CompareTo("Isolated") && electronTypePassed.CompareTo("UnIsolated") &&
      electronTypePassed.CompareTo("Tight") && electronTypePassed.CompareTo("Veto") && electronTypePassed.CompareTo("PtEtaCut") ){
    std::cout << "ERROR " << "<HistogrammingElectron::HistogrammingElectron()> " 
              << "Must pass Tight, Veto, or UnIsolated, PtEtaCut to constructor" << std::endl;
    exit(8);
  } //if
  electronType = electronTypePassed;

  SetEventContainer(obj);
} //HistogrammingElectron()

/******************************************************************************
 * HistogrammingElectron::~HistogrammingElectron()                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingElectron::~HistogrammingElectron()
{
  
} //HistogrammingElectron

/******************************************************************************
 * void HistogrammingElectron::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingElectron::BookHistogram(){
  
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
  // 1 Pt
  _hPtObj1 = DeclareTH1F("Ele1Pt", "Electron 1 P_{T}", ptBin, ptMin, ptMax);
  _hPtObj1->SetXAxisTitle("P_{T} (Ele 1) [GeV]");
  _hPtObj1->SetYAxisTitle("Events");
  // 1 Pt - Blown up
  _hPtBlowupObj1 = DeclareTH1F("Ele1PtBlowup","Electron 1 P_{T}",200,0.,100.);
  _hPtBlowupObj1->SetXAxisTitle("p_{T} (Ele 1) [GeV]");
  _hPtBlowupObj1->SetYAxisTitle("Events");
  // 1 Eta
  _hEtaObj1 = DeclareTH1F("Ele1Eta", "Electron 1 Pseudorapidity", etaBin, etaMin, etaMax);
  _hEtaObj1->SetXAxisTitle("#eta (Ele 1) ");
  _hEtaObj1->SetYAxisTitle("Events");
  // 1 Phi
  _hPhiObj1 = DeclareTH1F("Ele1Phi", "Electron 1 phi", phiBin, phiMin, phiMax);
  _hPhiObj1->SetXAxisTitle("#phi (Ele 1)");
  _hPhiObj1->SetYAxisTitle("Events");
  // 1 Pt-Eta
  _hPtEtaObj1 = DeclareTH2F("Ele1PtEta", "Electron 1 pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEtaObj1->SetXAxisTitle("P_{T} [GeV]");
  _hPtEtaObj1->SetYAxisTitle("#eta");
  // 1 Pt-Phi
  _hPtPhiObj1 = DeclareTH2F("Ele1PtPhi", "Electron 1 pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhiObj1->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhiObj1->SetYAxisTitle("#phi");
  // 1 Eta-Phi
  _hEtaPhiObj1 = DeclareTH2F("Ele1EtaPhi", "Electron 1 eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhiObj1->SetXAxisTitle("#eta [GeV]");
  _hEtaPhiObj1->SetYAxisTitle("#phi");

  // 2 Pt
  _hPtObj2 = DeclareTH1F("Ele2Pt", "Electron 2 P_{T}", ptBin, ptMin, ptMax);
  _hPtObj2->SetXAxisTitle("P_{T} (Ele 2) [GeV]");
  _hPtObj2->SetYAxisTitle("Events");
  // 2 Eta
  _hEtaObj2 = DeclareTH1F("Ele2Eta", "Electron 2 Pseudorapidity", etaBin, etaMin, etaMax);
  _hEtaObj2->SetXAxisTitle("#eta (Ele 2) ");
  _hEtaObj2->SetYAxisTitle("Events");
  // 2 Phi
  _hPhiObj2 = DeclareTH1F("Ele2Phi", "Electron 2 phi", phiBin, phiMin, phiMax);
  _hPhiObj2->SetXAxisTitle("#phi (Ele 2)");
  _hPhiObj2->SetYAxisTitle("Events");
  // 2 Pt-Eta
  _hPtEtaObj2 = DeclareTH2F("Ele2PtEta", "Electron 2 pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEtaObj2->SetXAxisTitle("P_{T} [GeV]");
  _hPtEtaObj2->SetYAxisTitle("#eta");
  // 2 Pt-Phi
  _hPtPhiObj2 = DeclareTH2F("Ele2PtPhi", "Electron 2 pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhiObj2->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhiObj2->SetYAxisTitle("#phi");
  // 2 Eta-Phi
  _hEtaPhiObj2 = DeclareTH2F("Ele2EtaPhi", "Electron 2 eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhiObj2->SetXAxisTitle("#eta [GeV]");
  _hEtaPhiObj2->SetYAxisTitle("#phi");


// All Electron jet Phi
  _hJetDeltaPhi = DeclareTH1F("EleJetPhi", "Leading Electron Jet Delta phi", 160, phiMin, phiMax);
  _hJetDeltaPhi->SetXAxisTitle("#phi (Ele)");
  _hJetDeltaPhi->SetYAxisTitle("Events");

// All Electron jet deltar
  _hJetDeltaR = DeclareTH1F("EleJetR", "Leading Electron Jet Delta R", 160, 0, 6);
  _hJetDeltaR->SetXAxisTitle("#deltaR (Ele)");
  _hJetDeltaR->SetYAxisTitle("Events");


  ////////////////////////////////////////////////////////////////
  // Kinematic Qunatities - All
  // All Electron Pt
  ////////////////////////////////////////////////////////////////
  _hPt = DeclareTH1F("ElePt", "Any Electron P_{T}", ptBin, ptMin, ptMax);
  _hPt->SetXAxisTitle("P_{T} (Ele) [GeV]");
  _hPt->SetYAxisTitle("Events");
  // All Electron Eta
  _hEta = DeclareTH1F("EleEta", "Any Electron Pseudorapidity", etaBin, etaMin, etaMax);
  _hEta->SetXAxisTitle("#eta (Ele) ");
  _hEta->SetYAxisTitle("Events");
  // All Electron Phi
  _hPhi = DeclareTH1F("ElePhi", "Any Electron phi", phiBin, phiMin, phiMax);
  _hPhi->SetXAxisTitle("#phi (Ele)");
  _hPhi->SetYAxisTitle("Events");
  // All Pt-Eta
  _hPtEta = DeclareTH2F("ElePtEta", "Any Electron pt-eta", ptBin2D, ptMin, ptMax, etaBin2D, etaMin, etaMax);
  _hPtEta->SetXAxisTitle("P_{T} [GeV]");
  _hPtEta->SetYAxisTitle("#eta");
  // All Pt-Phi
  _hPtPhi = DeclareTH2F("ElePtPhi", "Any Electron pt-phi", ptBin2D, ptMin, ptMax, phiBin2D, phiMin, phiMax);
  _hPtPhi->SetXAxisTitle("P_{T} [GeV]");
  _hPtPhi->SetYAxisTitle("#phi");
  // All Eta-Phi
  _hEtaPhi = DeclareTH2F("EleEtaPhi", "Any Electron eta-phi", etaBin2D, etaMin, etaMax, phiBin2D, phiMin, phiMax);
  _hEtaPhi->SetXAxisTitle("#eta");
  _hEtaPhi->SetYAxisTitle("#phi");

  // any Electron Charge
  _hCharge = DeclareTH1F("EleCharge","Any Electron charge",3,-1.5,1.5);
  _hCharge->SetXAxisTitle("Charge^{Ele}");
  _hCharge->SetYAxisTitle("Events");
  // any Electron E over P
  _hEoverP = DeclareTH1F("EleEoverP","Any Electron Energy over Momentum",50,0.,5.);
  _hEoverP->SetXAxisTitle("E/P^{Ele}");
  _hEoverP->SetYAxisTitle("Events");
  
  // selected Electron Etcone30/Pt
  _hEtCone30OverPt = DeclareTH1F("ElEtCone30OverPt", "Selected electron cone30 E_{T} over P_{T}", 100, 0., 1.);
  _hEtCone30OverPt->SetXAxisTitle("E_{T}^{cone30}/P_{T} (electron)");
  _hEtCone30OverPt->SetYAxisTitle("Events");
  // selected Electron Ptcone30/Pt
  _hPtCone30OverPt = DeclareTH1F("ElPtCone30OverPt", "Selected electron cone30 P_{T} over P_{T}", 100, 0., 1.);
  _hPtCone30OverPt->SetXAxisTitle("P_{T}^{cone30}/P_{T} (electron)");
  _hPtCone30OverPt->SetYAxisTitle("Events");
  //cout<<"end of HistogrammingElectron::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingElectron::Apply()                                         *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingElectron::Apply()
{
  //cout<<"Begin of HistogrammingElectron::Apply()"<<endl;
  
  // Get Event Tree
  EventContainer *evc = GetEventContainer();
  
  // Electrons class depends on 
  std::vector<Electron> electronVector;
  if(     "All"        == electronType) electronVector.assign(evc -> electrons.begin(),           evc -> electrons.end());
  else if("Veto"       == electronType) electronVector.assign(evc -> vetoElectrons.begin(),       evc -> vetoElectrons.end());
  else if("PtEtaCut"   == electronType) electronVector.assign(evc -> ptetaElectrons.begin(),      evc -> ptetaElectrons.end());
  else if("Tight"      == electronType) electronVector.assign(evc -> tightElectrons.begin(),      evc -> tightElectrons.end());
  else if("Isolated"   == electronType) electronVector.assign(evc -> isolatedElectrons.begin(),   evc -> isolatedElectrons.end());
  else if("UnIsolated" == electronType) electronVector.assign(evc -> unIsolatedElectrons.begin(), evc -> unIsolatedElectrons.end());
  else{
    std::cout << "ERROR " << "<HistogramminElectron::Apply()> " 
	      << "electronType must be All, Tight, Isolated, UnIsolated, or Veto, PtEtaCut" << std::endl;
    exit(8);
  } //else

  // Fill Histograms
  _hNObj -> Fill(electronVector.size());


  if(electronVector.size()>0) {

    // TMP
    //if( evc -> L1Electrons[0].GetDeltaRElectron() < 0.2 )   // DeltaR Recon
    //if( evc -> L1Electrons[0].GetDeltaRMCElectron() < 0.2 )  // DeltaR Truth
    //if( evc -> L1Electrons[0].GetDeltaRL2Electron() < 0.2 )  // DeltaR L2
    //if( evc -> L1Electrons[0].GetDeltaREvtFilterElectron() < 0.2 )  // DeltaR EF
    
    //if( evc -> L2Electrons[0].GetDeltaRElectron() < 0.2 )   // DeltaR Recon
    //if( evc -> L2Electrons[0].GetDeltaRMCElectron() < 0.2 )  // DeltaR Truth
    //if( evc -> L2Electrons[0].GetDeltaRL1Electron() < 0.2 )  // DeltaR L1
    //if( evc -> L2Electrons[0].GetDeltaREvtFilterElectron() < 0.2 )  // DeltaR EF
    
    //if( evc -> EFElectrons[0].GetDeltaRElectron() < 0.2 )   // DeltaR Recon
    //if( evc -> EFElectrons[0].GetDeltaRMCElectron() < 0.2 )  // DeltaR Truth
    //if( evc -> EFElectrons[0].GetDeltaRL2Electron() < 0.2 )  // DeltaR L2
    //if( evc -> EFElectrons[0].GetDeltaRL1Electron() < 0.2 )  // DeltaR L1
    //{


    _hPtObj1  -> Fill(electronVector[0].Pt());
    _hEtaObj1 -> Fill(electronVector[0].Eta());
    _hPhiObj1 -> Fill(electronVector[0].Phi());
    _hPtBlowupObj1  -> Fill(electronVector[0].Pt());
    
    _hPtEtaObj1  -> Fill(electronVector[0].Pt(),  electronVector[0].Eta());
    _hPtPhiObj1  -> Fill(electronVector[0].Pt(),  electronVector[0].Phi());
    _hEtaPhiObj1 -> Fill(electronVector[0].Eta(), electronVector[0].Phi());
    
    if(evc->jets.size()>0){
      _hJetDeltaR-> Fill(fabs(evc->jets[0].DeltaR(electronVector[0])));
      _hJetDeltaPhi-> Fill(fabs(evc->jets[0].DeltaPhi(electronVector[0])));
    }

    //      } //if DeltaR
  } // if size < 0
  if(electronVector.size()>1) {
    _hPtObj2  -> Fill(electronVector[1].Pt());
    _hEtaObj2 -> Fill(electronVector[1].Eta());
    _hPhiObj2 -> Fill(electronVector[1].Phi());

    _hPtEtaObj2  -> Fill(electronVector[1].Pt(),  electronVector[1].Eta());
    _hPtPhiObj2  -> Fill(electronVector[1].Pt(),  electronVector[1].Phi());
    _hEtaPhiObj2 -> Fill(electronVector[1].Eta(), electronVector[1].Phi());
  } // if size < 1

  // now loop over all electrons and histogram their properties
  for(int io=0;io<electronVector.size();io++) {
    _hPt  -> Fill(electronVector[io].Pt());
    _hEta -> Fill(electronVector[io].Eta());
    _hPhi -> Fill(electronVector[io].Phi());
    
    _hPtEta  -> Fill(electronVector[io].Pt(),  electronVector[io].Eta());
    _hPtPhi  -> Fill(electronVector[io].Pt(),  electronVector[io].Phi());
    _hEtaPhi -> Fill(electronVector[io].Eta(), electronVector[io].Phi());

    _hCharge->Fill(electronVector[io].GetCharge());
    //_hEoverP->Fill(electronVector[io].GetEoverP());

    //_hEtCone30OverPt->Fill(electronVector[io].GetEtCone30()/electronVector[io].Pt());
    //_hPtCone30OverPt->Fill(electronVector[io].GetPtCone30()/electronVector[io].Pt());
  }    

  //cout<<"End of HistogrammingElectron::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












