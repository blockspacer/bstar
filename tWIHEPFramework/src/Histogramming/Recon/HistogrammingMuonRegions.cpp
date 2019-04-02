/******************************************************************************
 * HistogramMuonRegions.cpp                                                   *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Books histograms in different pT, eta, and phi regions                     *
 * Only makes hists for Number of Events - used for Efficiency Calculations   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogramMuonRegions()    -- Parameterized Constructor                  *
 *    ~HistogramMuonRegions()   -- Destructor                                 *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                 -- Returns "HistogramMuonRegions"          *
 *                                                                            *
 * Private Data Members of HistogramMuonRegions class                         *
 * - lots of histograms of muon propreties                                    *
 *                                                                            *
 * History                                                                    *
 *      10 Dec 2007 - Created by P. Ryan                                      *
 *      22 Jan 2007 - Updated for new object definitions                      *
 *      2  Dec 2008 - P. Ryan - added 2D rejection histograms                 *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingMuonRegions.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingMuonRegions::HistogrammingMuonRegions(EventContainer *obj)    *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingMuonRegions::HistogrammingMuonRegions(EventContainer *obj, TString muonTypePassed)
{
  // Check muonType parameter
  if( muonTypePassed.CompareTo("All") && muonTypePassed.CompareTo("UnIsolated") && muonTypePassed.CompareTo("Isolated") &&
      muonTypePassed.CompareTo("Tight") && muonTypePassed.CompareTo("Veto") ){
    std::cout << "ERROR " << "<HistogrammingMuonRegions::HistogrammingMuonRegions()> " 
              << "Must pass Tight, Veto, Isolated, or UnIsolated to constructor" << std::endl;
    exit(8);
  } //if
  muonType = muonTypePassed;

  SetEventContainer(obj);

} //HistogrammingMuonRegions()

/******************************************************************************
 * HistogrammingMuonRegions::~HistogrammingMuonRegions()                      *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingMuonRegions::~HistogrammingMuonRegions()
{
  
} //HistogrammingMuonRegions

/******************************************************************************
 * void HistogrammingMuonRegions::BookHistogram()                             *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingMuonRegions::BookHistogram(){

  // Histogram bin sizes and ranges
  Int_t nBin    = 10;
  Double_t nMin = 0.0;
  Double_t nMax = 10.0;

  // For histogram names
  TString histName;
  TString histTitle;

  ////////////////////////////////////////////////////////////////
  // All Regions
  ////////////////////////////////////////////////////////////////
  _hNObj =  DeclareTH1F("MuN","Number of Muons", nBin, nMin, nMax);
  _hNObj -> SetXAxisTitle("N");
  _hNObj -> SetYAxisTitle("Events");


  ////////////////////////////////////////////////////////////////
  // pT Divisions
  ////////////////////////////////////////////////////////////////
  Double_t ptBase = 10.0;
  Double_t ptIncrement = 5.0;
  Int_t ptCutNumber = 15;

  ptCutsVector.clear();
  for(Int_t i = 0; i < ptCutNumber; i++) ptCutsVector.push_back(ptBase + i*ptIncrement); 

  // Separate pT Regions
  ptHistVector.clear();

  for(Int_t i = 0; i < ptCutNumber - 1; i++){
    // Histogram Names
    std::ostringstream histNameStream;
    histNameStream << "MuN_pt" << i+1;
    histName = histNameStream.str().c_str();
    // Histogram Titles
    std::ostringstream histTitleStream;
    histTitleStream << "Number of Muons - pt" << i+1;
    histTitle = histNameStream.str().c_str();
    // Book Histogram
    ptHistVector.push_back( (myTH1F*)(DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax)) );
    ptHistVector[i] -> SetXAxisTitle("N");
    ptHistVector[i] -> SetYAxisTitle("Events");
  } //for

  // pT Turn-on
  ptTurnOnHistVector.clear();

  for(Int_t i = 0; i < ptCutNumber ; i++){
    // Histogram Names
    std::ostringstream histNameStream;
    histNameStream << "MuN_ptTurnOn" << i+1;
    histName = histNameStream.str().c_str();
    // Histogram Titles
    std::ostringstream histTitleStream;
    histTitleStream << "Number of Muons - pt Turn On" << i+1;
    histTitle = histNameStream.str().c_str();
    // Book Histogram
    ptTurnOnHistVector.push_back( (myTH1F*)(DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax)) );
    ptTurnOnHistVector[i] -> SetXAxisTitle("N");
    ptTurnOnHistVector[i] -> SetYAxisTitle("Events");
  } //for

  ////////////////////////////////////////////////////////////////
  // Eta Divisions
  ////////////////////////////////////////////////////////////////
  Double_t etaBase = -2.5;
  Double_t etaIncrement = 0.5;
  Int_t etaCutNumber = 11;

  etaHistVector.clear();
  etaCutsVector.clear();

  for(Int_t i = 0; i < etaCutNumber; i++) etaCutsVector.push_back(etaBase + i*etaIncrement); 

  for(Int_t i = 0; i < etaCutNumber - 1; i++){
    // Histogram Names
    std::ostringstream histNameStream;
    histNameStream << "MuN_eta" << i+1;
    histName = histNameStream.str().c_str();
    // Histogram Titles
    std::ostringstream histTitleStream;
    histTitleStream << "Number of Muons - eta" << i+1;
    histTitle = histNameStream.str().c_str();
    // Book Histogram
    etaHistVector.push_back( (myTH1F*)(DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax)) );
    etaHistVector[i] -> SetXAxisTitle("N");
    etaHistVector[i] -> SetYAxisTitle("Events");
  } //for

  ////////////////////////////////////////////////////////////////
  // Phi Divisions
  ////////////////////////////////////////////////////////////////
  Double_t phiBase = 0.0;
  Double_t phiIncrement = TMath::Pi() / 6.0;
  Int_t phiCutNumber = 13;

  phiHistVector.clear();
  phiCutsVector.clear();

  for(Int_t i = 0; i < phiCutNumber; i++)  phiCutsVector.push_back(phiBase + i*phiIncrement - TMath::Pi()); 

  for(Int_t i = 0; i < phiCutNumber - 1; i++){
    // Histogram Names
    std::ostringstream histNameStream;
    histNameStream << "MuN_phi" << i+1;
    histName = histNameStream.str().c_str();
    // Histogram Titles
    std::ostringstream histTitleStream;
    histTitleStream << "Number of Muons - phi" << i+1;
    histTitle = histNameStream.str().c_str();
    // Book Histogram
    phiHistVector.push_back( (myTH1F*)(DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax)) );
    phiHistVector[i] -> SetXAxisTitle("N");
    phiHistVector[i] -> SetYAxisTitle("Events");
  } //for
  
  ////////////////////////////////////////////////////////////////
  // pT-Eta Divisions
  ////////////////////////////////////////////////////////////////
  ptEtaHistVectorVector.clear();

  for(Int_t i = 0; i < ptCutNumber - 1; i++){
    ptEtaHistVector.clear();
    for(Int_t j = 0; j < etaCutNumber - 1; j++){

      // Histogram Names
      std::ostringstream histNameStream;
      histNameStream << "MuN_pt" << i+1 << "_eta" << j+1;
      histName = histNameStream.str().c_str();

      // Histogram Titles
      std::ostringstream histTitleStream;
      histTitleStream << "Number of Muons - pt " << i+1 << " eta " << j+1;
      histTitle = histNameStream.str().c_str();

      // Book Histogram
      ptEtaHistVector.push_back( DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax) );

    } //for j
    // Push Vectors into Vector of Vectors
    ptEtaHistVectorVector.push_back( ptEtaHistVector );
  } //for i

  for(Int_t i = 0; i < ptCutNumber - 1; i++){
    for(Int_t j = 0; j < etaCutNumber - 1; j++){
      ptEtaHistVectorVector[i][j] -> SetXAxisTitle("N");
      ptEtaHistVectorVector[i][j] -> SetYAxisTitle("Events");
    } //for j
  } //for i


  ////////////////////////////////////////////////////////////////
  // pT-Phi Divisions
  ////////////////////////////////////////////////////////////////
  ptPhiHistVectorVector.clear();

  for(Int_t i = 0; i < ptCutNumber - 1; i++){
    ptPhiHistVector.clear();
    for(Int_t j = 0; j < phiCutNumber - 1; j++){

      // Histogram Names
      std::ostringstream histNameStream;
      histNameStream << "MuN_pt" << i+1 << "_phi" << j+1;
      histName = histNameStream.str().c_str();

      // Histogram Titles
      std::ostringstream histTitleStream;
      histTitleStream << "Number of Muons - pt " << i+1 << " phi " << j+1;
      histTitle = histNameStream.str().c_str();

      // Book Histogram
      ptPhiHistVector.push_back( DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax) );

    } //for j
    // Push Vectors into Vector of Vectors
    ptPhiHistVectorVector.push_back( ptPhiHistVector );
  } //for i

  for(Int_t i = 0; i < ptCutNumber - 1; i++){
    for(Int_t j = 0; j < phiCutNumber - 1; j++){
      ptPhiHistVectorVector[i][j] -> SetXAxisTitle("N");
      ptPhiHistVectorVector[i][j] -> SetYAxisTitle("Events");
    } //for j
  } //for i


  ////////////////////////////////////////////////////////////////
  // eta-Phi Divisions
  ////////////////////////////////////////////////////////////////
  etaPhiHistVectorVector.clear();

  for(Int_t i = 0; i < etaCutNumber - 1; i++){
    etaPhiHistVector.clear();
    for(Int_t j = 0; j < phiCutNumber - 1; j++){

      // Histogram Names
      std::ostringstream histNameStream;
      histNameStream << "MuN_eta" << i+1 << "_phi" << j+1;
      histName = histNameStream.str().c_str();

      // Histogram Titles
      std::ostringstream histTitleStream;
      histTitleStream << "Number of Muons - eta " << i+1 << " phi " << j+1;
      histTitle = histNameStream.str().c_str();

      // Book Histogram
      etaPhiHistVector.push_back( DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax) );

    } //for j
    // Push Vectors into Vector of Vectors
    etaPhiHistVectorVector.push_back( etaPhiHistVector );
  } //for i

  for(Int_t i = 0; i < etaCutNumber - 1; i++){
    for(Int_t j = 0; j < phiCutNumber - 1; j++){
      etaPhiHistVectorVector[i][j] -> SetXAxisTitle("N");
      etaPhiHistVectorVector[i][j] -> SetYAxisTitle("Events");
    } //for j
  } //for i


} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingMuonRegions::Apply()                                   *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingMuonRegions::Apply()
{
  //cout<<"Begin of HistogrammingMuonRegions::Apply()"<<endl;

  // Get Event Tree
  EventContainer *evc = GetEventContainer();

  ////////////////////////////////////////////////////////////////
  // Get Multiplicity, pT, eta, phi
  ////////////////////////////////////////////////////////////////
  Int_t multiplicity = evc -> tightMuons.size();
  Double_t pt  = -999.0;
  Double_t eta = -999.0;
  Double_t phi = -999.0;

  if( multiplicity > 0 ){
    //pt  = evc->muons[0].Pt();
    //eta = evc->muons[0].Eta();
    //phi = evc->muons[0].Phi();

    pt  = evc->tightMuons[0].Pt();
    eta = evc->tightMuons[0].Eta();
    phi = evc->tightMuons[0].Phi();
  } //if

  ////////////////////////////////////////////////////////////////
  // Fill histograms according to regions
  ////////////////////////////////////////////////////////////////
  for(Int_t i = 0; i < (ptCutsVector.size()  - 1); i++) if( (pt  > ptCutsVector[i]) && (pt  < ptCutsVector[i+1])   )  ptHistVector[i]        -> Fill(multiplicity);
  for(Int_t i = 0; i < (etaCutsVector.size() - 1); i++) if( (eta > etaCutsVector[i]) && (eta < etaCutsVector[i+1]) )  etaHistVector[i]       -> Fill(multiplicity);
  for(Int_t i = 0; i < (phiCutsVector.size() - 1); i++) if( (phi > phiCutsVector[i]) && (phi < phiCutsVector[i+1]) )  phiHistVector[i]       -> Fill(multiplicity);

  // pt-eta
  for(Int_t i = 0; i < (ptCutsVector.size() -1); i++){
    for(Int_t j = 0; j < (etaCutsVector.size() -1); j++){
      if( (pt  > ptCutsVector[i]) && (pt  < ptCutsVector[i+1]) && (eta > etaCutsVector[j]) && (eta < etaCutsVector[j+1])  ) ptEtaHistVectorVector[i][j] -> Fill(multiplicity);
    } //for j
  } //for i

  // pt-phi
  for(Int_t i = 0; i < (ptCutsVector.size() -1); i++){
    for(Int_t j = 0; j < (phiCutsVector.size() -1); j++){
      if( (pt  > ptCutsVector[i]) && (pt  < ptCutsVector[i+1]) && (phi > phiCutsVector[j]) && (phi < phiCutsVector[j+1])  ) ptPhiHistVectorVector[i][j] -> Fill(multiplicity);
    } //for j
  } //for i

  // eta-phi
  for(Int_t i = 0; i < (etaCutsVector.size() -1); i++){
    for(Int_t j = 0; j < (phiCutsVector.size() -1); j++){
      if( (eta  > etaCutsVector[i]) && (eta  < etaCutsVector[i+1]) && (phi > phiCutsVector[j]) && (phi < phiCutsVector[j+1])  ) etaPhiHistVectorVector[i][j] -> Fill(multiplicity);
    } //for j
  } //for i

  ////////////////////////////////////////////////////////////////
  // Fill histograms for turn-on curves
  ////////////////////////////////////////////////////////////////
  //  Int_t multiplicityTurnOn = evc -> isolatedMuons.size();
  Int_t multiplicityTurnOn = evc -> tightMuons.size();
  Double_t ptTurnOn  = -999.0;
  Double_t etaTurnOn = -999.0;
  Double_t phiTurnOn = -999.0;

  if( multiplicityTurnOn > 0 ){
    ptTurnOn  = evc -> tightMuons[0].Pt();
    etaTurnOn = evc -> tightMuons[0].Eta();
    phiTurnOn = evc -> tightMuons[0].Phi();
  } //if

  // TMP
  //if(evc -> L1Muons.size() > 0){
  // if( evc -> L1Muons[0].GetDeltaRMuon() < 0.2 )   // DeltaR Recon
      //if( evc -> L1Muons[0].GetDeltaRMCMuon() < 0.2 )  // DeltaR Truth
      //if( evc -> L1Muons[0].GetDeltaRL2Muon() < 0.2 )  // DeltaR L2
      //if( evc -> L1Muons[0].GetDeltaREvtFilterMuon() < 0.2 )  // DeltaR EF
  //  {for(Int_t i = 0; i < (ptCutsVector.size()); i++) if( (ptTurnOn  > ptCutsVector[i]) )  ptTurnOnHistVector[i]  -> Fill(multiplicityTurnOn);} //if DeltaR
  // } //if size
      

  //if(evc -> L2Muons.size() > 0){
    //if( evc -> L2Muons[0].GetDeltaRMuon() < 0.2 )   // DeltaR Recon
    //if( evc -> L2Muons[0].GetDeltaRMCMuon() < 0.2 )  // DeltaR Truth
    //if( evc -> L2Muons[0].GetDeltaRL1Muon() < 0.2 )  // DeltaR L1
    //if( evc -> L2Muons[0].GetDeltaREvtFilterMuon() < 0.2 )  // DeltaR EF
    //{for(Int_t i = 0; i < (ptCutsVector.size()); i++) if( (ptTurnOn  > ptCutsVector[i]) )  ptTurnOnHistVector[i]  -> Fill(multiplicityTurnOn);} //if DeltaR
  //} //if size  

  //if(evc -> EFMuons.size() > 0){
    //    if(evc -> L1Muons[0].GetDeltaRMuon() > 0.2) std::cout << "@@@@ " << evc -> L1Muons[0].GetDeltaRMuon() << std::endl;
    //if( evc -> EFMuons[0].GetDeltaRMuon() < 0.2 )   // DeltaR Recon
    //if( evc -> EFMuons[0].GetDeltaRMCMuon() < 0.2 )  // DeltaR Truth
    //if( evc -> EFMuons[0].GetDeltaRL2Muon() < 0.2 )  // DeltaR L2
    //if( evc -> EFMuons[0].GetDeltaRL1Muon() < 0.2 )  // DeltaR L1
    //{for(Int_t i = 0; i < (ptCutsVector.size()); i++) if( (ptTurnOn  > ptCutsVector[i]) )  ptTurnOnHistVector[i]  -> Fill(multiplicityTurnOn);} //if DeltaR
  //} //if size

  //  else{//    std::cout << "#### " << evc -> L1Muons[0].GetDeltaRMuon() << std::endl; // tmp}

  // Old - histograms events with 0 multiplicity
  //for(Int_t i = 0; i < (ptCutsVector.size()); i++){
  // if( (ptTurnOn  > ptCutsVector[i]) || -999.0 == ptTurnOn ){
  //  ptTurnOnHistVector[i]  -> Fill(multiplicityTurnOn);
      //      if( -999.0 == ptTurnOn) std::cout << "WARNING: <HistogrammingMuonRegions:Apply()> Muon Multiplicity is 0" << std::endl;
  // } //if
  //} //for

  //  if( -999.0 == ptTurnOn) std::cout << "WARNING: <HistogrammingMuonRegions:Apply()> Muon Multiplicity is 0.  Event not histogrammed." << std::endl;
  for(Int_t i = 0; i < (ptCutsVector.size()); i++){
    if( (ptTurnOn  > ptCutsVector[i]) && -999.0 != ptTurnOn ){
      ptTurnOnHistVector[i]  -> Fill(multiplicityTurnOn);
    } //if
  } //for

  // Fill
  _hNObj -> Fill(multiplicity);
  
  ////////////////tmp
  //for(Int_t i = 0; i < (ptCutsVector.size() ); i++) std::cout << "Histogram " << i+1 << "    pT > " << ptCutsVector[i] << std::endl;
  ////////////// end tmp
  
  //cout<<"End of HistogrammingMuonRegions::Apply()"<<endl;
  return kTRUE;  
  
} //Apply











