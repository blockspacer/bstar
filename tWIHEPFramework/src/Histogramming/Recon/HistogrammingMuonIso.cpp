/**********************************************************************************
 * HistogrammingMuonIso.pp                                                        *
 *                                                                                *
 * Books and fills histograms                                                     *
 * Books histograms related to isolation                                          *
 *                                                                                *
 * Derived from HistoCut which is in turn derived from BaseCut                    *
 *                                                                                *
 * Public Member Functions of AnalysisMain class                                  *
 *    HistogrammingMuonRegions()    -- Parameterized Constructor                  *
 *    ~HistogrammingMuonRegions()   -- Destructor                                 *
 *    BookHistogramming()                   -- Book histograms                    *
 *    Apply()                               -- Fill histograms only (No Cuts)     *
 *    GetCutName()                 -- Returns "HistogrammingMuonRegions"          *
 *                                                                                *
 * Private Data Members of HistogrammingMuonRegions class                         *
 * - lots of histograms of muon propreties                                        *
 *                                                                                *
 * History                                                                        *
 *      17 June 2015 - Created by Huaqiao ZHANG                                   *
 *                                                                                *
 **********************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingMuonIso.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingMuonIso::HistogrammingMuonIso(EventContainer *obj)            *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingMuonIso::HistogrammingMuonIso(EventContainer *obj, TString muonTypePassed)
{
  // Check muonType parameter
  if( muonTypePassed.CompareTo("All") && muonTypePassed.CompareTo("UnIsolated") && muonTypePassed.CompareTo("Isolated") &&
      muonTypePassed.CompareTo("Tight") && muonTypePassed.CompareTo("Veto") ){
    std::cout << "ERROR " << "<HistogrammingMuonIso::HistogrammingMuonIso()> " 
              << "Must pass Tight, Veto, Isolated, or UnIsolated to constructor" << std::endl;
    exit(8);
  } //if
  muonType = muonTypePassed;

  SetEventContainer(obj);

} //HistogrammingMuonIso()

/******************************************************************************
 * HistogrammingMuonIso::~HistogrammingMuonIso()                              *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingMuonIso::~HistogrammingMuonIso()
{
  
} //HistogrammingMuonIso

/******************************************************************************
 * void HistogrammingMuonIso::BookHistogram()                                 *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingMuonIso::BookHistogram(){

  // Histogram bin sizes and ranges
  Int_t nBin    = 10;
  Double_t nMin = 0.0;
  Double_t nMax = 10.0;

  Int_t EtConeBin    = 25;
  Double_t EtConeMin = -5.0;
  Double_t EtConeMax = 20.0;

  Int_t DeltaRBin    = 16;
  Double_t DeltaRMin = 0.0;
  Double_t DeltaRMax = 8.0;

  Int_t DeltaRFineBin    = 50;
  Double_t DeltaRFineMin = 0.0;
  Double_t DeltaRFineMax = 0.5;

  Int_t EtaBin    = 50;
  Double_t EtaMin = -2.5;
  Double_t EtaMax = 2.5;

  // For histogram names
  TString histName;
  TString histTitle;

  ////////////////////////////////////////////////////////////////
  // All Regions
  ////////////////////////////////////////////////////////////////
  _hNObj =  DeclareTH1F("MuN","Number of Muons", nBin, nMin, nMax);
  _hNObj -> SetXAxisTitle("N");
  _hNObj -> SetYAxisTitle("Events");

  _hEtCone =  DeclareTH1F("MuEtCone","Muon EtCone", EtConeBin, EtConeMin, EtConeMax);
  _hEtCone -> SetXAxisTitle("EtCone20");
  _hEtCone -> SetYAxisTitle("Events");

  _hDeltaR =  DeclareTH1F("MuDeltaR","Muon DeltaR", DeltaRBin, DeltaRMin, DeltaRMax);
  _hDeltaR -> SetXAxisTitle("#DeltaR Muon - Closest Jet");
  _hDeltaR -> SetYAxisTitle("Events");

  _hDeltaRFine =  DeclareTH1F("MuDeltaRFine","Muon DeltaR Fine", DeltaRFineBin, DeltaRFineMin, DeltaRFineMax);
  _hDeltaRFine -> SetXAxisTitle("#DeltaR Muon - Closest Jet");
  _hDeltaRFine -> SetYAxisTitle("Events");

  _hEta =  DeclareTH1F("MuEta","Muon Eta", EtaBin, EtaMin, EtaMax);
  _hEta -> SetXAxisTitle("Muon #eta");
  _hEta -> SetYAxisTitle("Events");

  _hEtaAbs =  DeclareTH1F("MuEtaAbs","Abs Muon Eta", EtaBin, EtaMin, EtaMax);
  _hEtaAbs -> SetXAxisTitle("Muon |#eta|");
  _hEtaAbs -> SetYAxisTitle("Events");


  ////////////////////////////////////////////////////////////////
  // EtCone Divisions
  ////////////////////////////////////////////////////////////////
  Double_t EtConeBase      = -5.0;
  Double_t EtConeIncrement = 1.0;
  Int_t EtConeCutNumber    = 26;

  EtConeCutsVector.clear();
  for(Int_t i = 0; i < EtConeCutNumber; i++) EtConeCutsVector.push_back(EtConeBase + i * EtConeIncrement); 

  // Separate pT Regions
  EtConeHistVector.clear();

  for(Int_t i = 0; i < EtConeCutNumber; i++){
    // Histogram Names
    std::ostringstream histNameStream;
    histNameStream << "MuN_EtCone" << i+1;
    histName = histNameStream.str().c_str();
    // Histogram Titles
    std::ostringstream histTitleStream;
    histTitleStream << "Number of Muons - EtCone" << i+1;
    histTitle = histNameStream.str().c_str();
    // Book Histogram
    EtConeHistVector.push_back( (myTH1F*)(DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax)) );
    EtConeHistVector[i] -> SetXAxisTitle("N");
    EtConeHistVector[i] -> SetYAxisTitle("Events");
  } //for

  ////////////////////////////////////////////////////////////////
  // DeltaR Divisions
  ////////////////////////////////////////////////////////////////
  Double_t DeltaRBase      = 0.0;
  Double_t DeltaRIncrement = 0.5;
  Int_t DeltaRCutNumber    = 17;

  DeltaRHistVector.clear();
  DeltaRCutsVector.clear();

  for(Int_t i = 0; i < DeltaRCutNumber; i++) DeltaRCutsVector.push_back(DeltaRBase + i*DeltaRIncrement); 

  for(Int_t i = 0; i < DeltaRCutNumber; i++){
    // Histogram Names
    std::ostringstream histNameStream;
    histNameStream << "MuN_DeltaR" << i+1;
    histName = histNameStream.str().c_str();
    // Histogram Titles
    std::ostringstream histTitleStream;
    histTitleStream << "Number of Muons - DeltaR" << i+1;
    histTitle = histNameStream.str().c_str();
    // Book Histogram
    DeltaRHistVector.push_back( (myTH1F*)(DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax)) );
    DeltaRHistVector[i] -> SetXAxisTitle("N");
    DeltaRHistVector[i] -> SetYAxisTitle("Events");
  } //for

  ////////////////////////////////////////////////////////////////
  // DeltaRFine Divisions
  ////////////////////////////////////////////////////////////////
  Double_t DeltaRFineBase      = 0.0;
  Double_t DeltaRFineIncrement = 0.02;
  Int_t DeltaRFineCutNumber    = 21;

  DeltaRFineHistVector.clear();
  DeltaRFineCutsVector.clear();

  for(Int_t i = 0; i < DeltaRFineCutNumber; i++) DeltaRFineCutsVector.push_back(DeltaRFineBase + i*DeltaRFineIncrement); 

  for(Int_t i = 0; i < DeltaRFineCutNumber; i++){
    // Histogram Names
    std::ostringstream histNameStream;
    histNameStream << "MuN_DeltaRFine" << i+1;
    histName = histNameStream.str().c_str();
    // Histogram Titles
    std::ostringstream histTitleStream;
    histTitleStream << "Number of Muons - DeltaRFine" << i+1;
    histTitle = histNameStream.str().c_str();
    // Book Histogram
    DeltaRFineHistVector.push_back( (myTH1F*)(DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax)) );
    DeltaRFineHistVector[i] -> SetXAxisTitle("N");
    DeltaRFineHistVector[i] -> SetYAxisTitle("Events");
  } //for


  ////////////////////////////////////////////////////////////////
  // EtaAbs Divisions
  ////////////////////////////////////////////////////////////////
  Double_t EtaAbsBase      = 0.0;
  Double_t EtaAbsIncrement = 0.1;
  Int_t EtaAbsCutNumber    = 26;

  EtaAbsHistVector.clear();
  EtaAbsCutsVector.clear();

  for(Int_t i = 0; i < EtaAbsCutNumber; i++) EtaAbsCutsVector.push_back(EtaAbsBase + i*EtaAbsIncrement); 

  for(Int_t i = 0; i < EtaAbsCutNumber; i++){
    // Histogram Names
    std::ostringstream histNameStream;
    histNameStream << "MuN_EtaAbs" << i+1;
    histName = histNameStream.str().c_str();
    // Histogram Titles
    std::ostringstream histTitleStream;
    histTitleStream << "Number of Muons - EtaAbs" << i+1;
    histTitle = histNameStream.str().c_str();
    // Book Histogram
    EtaAbsHistVector.push_back( (myTH1F*)(DeclareTH1F(histName.Data(), histTitle.Data(), nBin, nMin, nMax)) );
    EtaAbsHistVector[i] -> SetXAxisTitle("N");
    EtaAbsHistVector[i] -> SetYAxisTitle("Events");
  } //for

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingMuonIso::Apply()                                       *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingMuonIso::Apply()
{
  //cout<<"Begin of HistogrammingMuonIso::Apply()"<<endl;

  // Get Event Tree
  EventContainer *evc = GetEventContainer();

  ////////////////////////////////////////////////////////////////
  // Get Multiplicity, pT, eta, phi
  ////////////////////////////////////////////////////////////////
  Int_t multiplicity = evc -> muons.size();
  Double_t EtCone = -999.0;
  Double_t DeltaR = -999.0;
  Double_t Eta    = -999.0;
  Double_t EtaAbs = -999.0;

  if( multiplicity > 0 ){
    // Only dong EtCone20
    Eta    = evc -> muons[0].Eta();
    EtaAbs = fabs(Eta);
  } //if

  // Fill Histograms which cover entire range
  _hNObj -> Fill(multiplicity);
  if( multiplicity > 0 ){
    _hEtCone     -> Fill(EtCone);
    _hDeltaR     -> Fill(DeltaR);
    _hDeltaRFine -> Fill(DeltaR);
    _hEta        -> Fill(Eta);
    _hEtaAbs     -> Fill(EtaAbs);
  } //if

  ////////////////////////////////////////////////////////////////
  // Fill histograms according to regions
  ////////////////////////////////////////////////////////////////
  if( multiplicity > 0 ){
    for(Int_t i = 0; i < (EtConeCutsVector.size() - 1); i++) if( (EtCone  < EtConeCutsVector[i]) ) EtConeHistVector[i] -> Fill(multiplicity);
    for(Int_t i = 0; i < (DeltaRCutsVector.size() - 1); i++) if( (DeltaR  < DeltaRCutsVector[i]) ) DeltaRHistVector[i] -> Fill(multiplicity);
    for(Int_t i = 0; i < (EtaAbsCutsVector.size() - 1); i++) if( (EtaAbs  < EtaAbsCutsVector[i]) ) EtaAbsHistVector[i] -> Fill(multiplicity);

    for(Int_t i = 0; i < (DeltaRFineCutsVector.size() - 1); i++) if( (DeltaR < DeltaRFineCutsVector[i]) ) DeltaRFineHistVector[i] -> Fill(multiplicity);
  } //if
  
  return kTRUE;  
  
} //Apply











