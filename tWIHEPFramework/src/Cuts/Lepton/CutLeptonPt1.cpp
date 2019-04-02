/******************************************************************************
 * CutLeptonPt1.hpp                                                               *
 *                                                                            *
 * Cuts on muon Number (can cut on All, UnIsolated, Tight, or Veto Muons)     *
 * Must pass All, UnIsolated, Tight, or Veto to constructor                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutLeptonPt1 class                                  *
 *    CutLeptonPt1()                     -- Parameterized Constructor             *
 *    ~CutLeptonPt1()                    -- Destructor                            *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutLeptonPt1"                    *
 *                                                                            *
 * Private Data Members of CutLeptonPt1 class                                     *
 *    myTH1F* _hMuonNumberBefore;    -- Hist mu Number before cut             *
 *    myTH1F* _hMuonNumberAfter;     -- Hist mu Number of jets after cut      *
 *                                                                            *
 *    Int_t _MuonNumberMin;          -- Minimum Muon Number                   *
 *    Int_t _MuonNumberMax;          -- Maximum Muon Number                   *
 *                                                                            *
 * History                                                                    *
 *      15 Jan 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Cuts/Lepton/CutLeptonPt1.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * CutLeptonPt1::CutLeptonPt1(EventContainer *EventContainerObj, TString muonType)    *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
CutLeptonPt1::CutLeptonPt1(EventContainer *EventContainerObj, TString leptonType)
{
  // Check muonType parameter
  if( leptonType.CompareTo("All") && leptonType.CompareTo("UnIsolated") && leptonType.CompareTo("Isolated") && 
      leptonType.CompareTo("Tight") && leptonType.CompareTo("PtEtaCut") && leptonType.CompareTo("Veto") ){
    std::cout << "ERROR " << "<CutLeptonPt1::CutLeptonPt1()> " 
	      << "Must pass All, Tight, PtEtaCut, Veto, Isolated, or UnIsolated to constructor" << std::endl;
    exit(8);
  } //if
  _leptonType = leptonType;

  // Set Event Container
  SetEventContainer(EventContainerObj);
} // CutLeptonPt1


/******************************************************************************
 * CutLeptonPt1::~CutLeptonPt1()                                                *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
CutLeptonPt1::~CutLeptonPt1()
{
  
}//~CutLeptonPt1

/******************************************************************************
 * void CutLeptonPt1::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutLeptonPt1::BookHistogram(){
  
  // ***********************************************
  // Make Strings for histogram titles and labels
  // ***********************************************  

  // Histogram Before Cut
  std::ostringstream histNameBeforeStream;
  histNameBeforeStream << _leptonType << "LeadingLeptonPtBefore";
  TString histNameBefore = histNameBeforeStream.str().c_str();

  std::ostringstream histTitleBeforeStream;
  histTitleBeforeStream << _leptonType << "Leading Lepton Pt Before Cut";
  TString histTitleBefore = histTitleBeforeStream.str().c_str();

  // Histogram After Cut
  std::ostringstream histNameAfterStream;
  histNameAfterStream << _leptonType << "LeadingLeptonPtAfter";
  TString histNameAfter = histNameAfterStream.str().c_str();

  std::ostringstream histTitleAfterStream;
  histTitleAfterStream << _leptonType << "Leading Lepton Pt After Cut";
  TString histTitleAfter = histTitleAfterStream.str().c_str();

  // ***********************************************
  // Book Histograms
  // ***********************************************  

  // Histogram before cut
  _hLeadingLeptonPtBefore =  DeclareTH1F(histNameBefore.Data(), histTitleBefore.Data(), 100, 0.0, 200.);
  _hLeadingLeptonPtBefore -> SetXAxisTitle("p_{T,leading}");
  _hLeadingLeptonPtBefore -> SetYAxisTitle("Events");

  // Histogram after cut
  _hLeadingLeptonPtAfter =  DeclareTH1F(histNameAfter.Data(), histTitleAfter.Data(), 100, 0.0, 200.);
  _hLeadingLeptonPtAfter -> SetXAxisTitle("p_{T,leading}");
  _hLeadingLeptonPtAfter -> SetYAxisTitle("Events");

  // ***********************************************
  // Get cuts from configuration file
  // ***********************************************  

  // Get configuration file
  EventContainer *EventContainerObj = GetEventContainer();
  TEnv *config = EventContainerObj -> GetConfig();

  // Make name in configuration file depending upon muon type (all, veto, tight, unIsolated)
  std::ostringstream configMinStream;
  configMinStream << "Cut.Muon." << _leptonType.Data() << ".LeadingPt";
  TString configMin = configMinStream.str().c_str();

  std::ostringstream configMaxStream;
  configMaxStream << "Cut.Electron." << _leptonType.Data() << ".LeadingPt";
  TString configMax = configMaxStream.str().c_str();

  // Use configuration to set min and max number of jets to cut on
  _LeadingMuonPtCut  = config -> GetValue(configMin.Data(), 999);
  _LeadingElectronPtCut  = config -> GetValue(configMax.Data(), 999);
  
  // ***********************************************
  // Add these cuts to the cut flow table
  // ***********************************************
  ostringstream cutFlowTitleStream;
  ostringstream cutFlowNameStream;
  TString cutFlowTitle;
  TString cutFlowName;

  cutFlowTitleStream << _leptonType.Data() << " Leading lepton : P_{T,#mu} > " << _LeadingMuonPtCut << " p_{T,e} > " << _LeadingElectronPtCut;
  cutFlowTitle = cutFlowTitleStream.str().c_str();

  cutFlowNameStream << _leptonType.Data() << "Lepton.LeadingPt";
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName.Data(), cutFlowTitle.Data());

}//BookHistograms()

/******************************************************************************
 * Bool_t CutMuonN::Apply()                                                *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Cut on the number of electrons and the p_T of these electrons              *
 * Apply separate p_T cuts to the first and second (don't worry about others) *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t CutLeptonPt1::Apply()
{

  // ********************************
  // Get the leading Pt of the leptons
  // ********************************

  // Leading lepton pt
  Float_t LeadingLeptonPt = 0.;
  
  //Flag
  Bool_t PassesLeadingLetonPt = kFALSE;
  
  //Get Event Container
  EventContainer *EventContainerObj = GetEventContainer();

  std::vector<Muon> muonVector;
  std::vector<Electron> electronVector;

  //Assign collections to the above defined vectors
  if(     "All"        == _leptonType) {
    muonVector.assign(EventContainerObj -> muons.begin(),           EventContainerObj -> muons.end());
    electronVector.assign(EventContainerObj -> electrons.begin(),           EventContainerObj -> electrons.end());
  }
  else if("Veto"       == _leptonType) {
    muonVector.assign(EventContainerObj -> vetoMuons.begin(),       EventContainerObj -> vetoMuons.end());
    electronVector.assign(EventContainerObj -> vetoElectrons.begin(),       EventContainerObj -> vetoElectrons.end());
  }
  else if("PtEtaCut"   == _leptonType) {
    muonVector.assign(EventContainerObj -> ptetaMuons.begin(),      EventContainerObj -> ptetaMuons.end());
    electronVector.assign(EventContainerObj -> ptetaElectrons.begin(),      EventContainerObj -> ptetaElectrons.end());
  }
  else if("Tight"      == _leptonType) {
    muonVector.assign(EventContainerObj -> tightMuons.begin(),      EventContainerObj -> tightMuons.end());
    electronVector.assign(EventContainerObj -> tightElectrons.begin(),      EventContainerObj -> tightElectrons.end());
  }
  else if("Isolated"   == _leptonType) {
    muonVector.assign(EventContainerObj -> isolatedMuons.begin(),   EventContainerObj -> isolatedMuons.end());
    electronVector.assign(EventContainerObj -> isolatedElectrons.begin(),   EventContainerObj -> isolatedElectrons.end());
  }
  else if("UnIsolated" == _leptonType) {
    muonVector.assign(EventContainerObj -> unIsolatedMuons.begin(), EventContainerObj -> unIsolatedMuons.end());
    electronVector.assign(EventContainerObj -> unIsolatedElectrons.begin(), EventContainerObj -> unIsolatedElectrons.end());
  }
  else{
    std::cout << "ERROR " << "<HistogramminMuon::Apply()> "
	      << "muonType must be All, Tight, Veto, Isolated, or UnIsolated, PtEtaCut" << std::endl;
    exit(8);
  } //else                     

  // Now go through and see if there's a lepton passing this cut
  for (auto const muon : muonVector){
    if (muon.Pt() > LeadingLeptonPt) LeadingLeptonPt = muon.Pt();  
    if (muon.Pt() > _LeadingMuonPtCut){
      PassesLeadingLetonPt = kTRUE;
    }
  }

  for (auto const electron : electronVector){
    if (electron.Pt() > LeadingLeptonPt) LeadingLeptonPt = electron.Pt();  
    if (electron.Pt() > _LeadingElectronPtCut){
      PassesLeadingLetonPt = kTRUE;
    }
  }
  
  _hLeadingLeptonPtBefore->Fill(LeadingLeptonPt);

  ostringstream cutFlowNameStream;
  TString cutFlowName;

  cutFlowNameStream << _leptonType.Data() << "Lepton.LeadingPt";
  cutFlowName = cutFlowNameStream.str().c_str();
  
  if (PassesLeadingLetonPt) {
    _hLeadingLeptonPtAfter->Fill(LeadingLeptonPt);
    GetCutFlowTable()->PassCut(cutFlowName.Data());
  }
  else {
    GetCutFlowTable()->FailCut(cutFlowName.Data());
  }

  return PassesLeadingLetonPt;

 
} //Apply












