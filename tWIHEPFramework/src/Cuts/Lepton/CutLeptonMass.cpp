/******************************************************************************
 * CutLeptonMass.hpp                                                             *
 *                                                                            *
 * Cuts on lepton Number (can cut on All, UnIsolated, Tight, or Veto Leptons) *
 * Must pass All, UnIsolated, Tight, or Veto to constructor                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutLeptonMass class                                *
 *    CutLeptonMass()                     -- Parameterized Constructor           *
 *    ~CutLeptonMass()                    -- Destructor                          *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutLeptonMass"                  *
 *                                                                            *
 * Private Data Members of CutLeptonMass class                                   *
 *    myTH1F* _hLeptonNumberBefore;    -- Hist lepton Number before cut       *
 *    myTH1F* _hLeptonNumberAfter;     -- Hist lepton Number of jets after cut*
 *                                                                            *
 *    Int_t _LeptonNumberMin;          -- Minimum Lepton Number               *
 *    Int_t _LeptonNumberMax;          -- Maximum Lepton Number               *
 *                                                                            *
 * History                                                                    *
 *      24 Mar 2009 - Created by P. Ryan                                      *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Cuts/Lepton/CutLeptonMass.hpp"
#include<iostream>

using namespace std;

/*********************************************************************************
 * CutLeptonMass::CutLeptonMass(EventContainer *EventContainerObj, TString leptonType) *
 *                                                                               *
 * Parameterized Constructor                                                     *
 *                                                                               *
 * Input:  Event Object class                                                    *
 * Output: None                                                                  *
 ******************************************************************************/
CutLeptonMass::CutLeptonMass(EventContainer *EventContainerObj, TString leptonTypePassed, TString massCutPost)
{
  // Check leptonType parameter
  if( leptonTypePassed.CompareTo("All") && leptonTypePassed.CompareTo("UnIsolated") && leptonTypePassed.CompareTo("Isolated") && 
      leptonTypePassed.CompareTo("Tight") && leptonTypePassed.CompareTo("Veto") ){
    std::cout << "ERROR " << "<CutLeptonMass::CutLeptonMass()> " 
	      << "Must pass All, Tight, Veto, Isolated, or UnIsolated to constructor" << std::endl;
    exit(8);
  } //if
  leptonType = leptonTypePassed;
  massCutPostfix = massCutPost;

  // Set Event Container
  SetEventContainer(EventContainerObj);
} // CutLeptonMass


/******************************************************************************
 * CutLeptonMass::~CutLeptonMass()                                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
CutLeptonMass::~CutLeptonMass()
{
  
}//~CutLeptonMass

/******************************************************************************
 * void CutLeptonMass::BookHistogram()                                           *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutLeptonMass::BookHistogram(){
  
  // ***********************************************
  // Make Strings for histogram titles and labels
  // ***********************************************  

  // Histogram Before Cut
  std::ostringstream histNameBeforeStream;
  histNameBeforeStream << leptonType << "DileptonMassBefore" << massCutPostfix;
  TString histNameBefore = histNameBeforeStream.str().c_str();

  std::ostringstream histTitleBeforeStream;
  histTitleBeforeStream << leptonType << "Dilepton Mass Before Cut " << massCutPostfix;
  TString histTitleBefore = histTitleBeforeStream.str().c_str();

  // Histogram After Cut
  std::ostringstream histNameAfterStream;
  histNameAfterStream << leptonType << "DileptonMassAfter" << massCutPostfix;
  TString histNameAfter = histNameAfterStream.str().c_str();

  std::ostringstream histTitleAfterStream;
  histTitleAfterStream << leptonType << "Dilepton Mass After Cut " << massCutPostfix;
  TString histTitleAfter = histTitleAfterStream.str().c_str();

  // ***********************************************
  // Book Histograms
  // ***********************************************  

  // Histogram before cut
  _hLeptonMassBefore =  DeclareTH1F(histNameBefore.Data(), histTitleBefore.Data(), 10, 0.0, 10.0);
  _hLeptonMassBefore -> SetXAxisTitle("m_{ll}");
  _hLeptonMassBefore -> SetYAxisTitle("Events");

  // Histogram after cut
  _hLeptonMassAfter =  DeclareTH1F(histNameAfter.Data(), histTitleAfter.Data(), 10, 0.0, 10);
  _hLeptonMassAfter -> SetXAxisTitle("m_{ll}");
  _hLeptonMassAfter -> SetYAxisTitle("Events");

  // ***********************************************
  // Get cuts from configuration file
  // ***********************************************  

  // Get configuration file
  EventContainer *EventContainerObj = GetEventContainer();
  TEnv *config = EventContainerObj -> GetConfig();

  // Make name in configuration file depending upon lepton type (all, veto, tight, unIsolated)
  std::ostringstream configMinStream;
  configMinStream << "Cut.Dilepton." << leptonType.Data() << ".Mass.Min" << massCutPostfix;
  TString configMin = configMinStream.str().c_str();

  std::ostringstream configMaxStream;
  configMaxStream << "Cut.Dilepton." << leptonType.Data() << ".Mass.Max" << massCutPostfix;
  TString configMax = configMaxStream.str().c_str();

  std::ostringstream configRejectStream;
  configRejectStream << "Cut.Dilepton." << leptonType.Data() << ".Mass.Reject" << massCutPostfix;
  TString configReject = configRejectStream.str().c_str();

  // Use configuration to set min and max number of jets to cut on
  _LeptonMassMin  = config -> GetValue(configMin.Data(), 999.);
  _LeptonMassMax  = config -> GetValue(configMax.Data(), 999.);
  
  // Reject mass region if true, accept if false
  _LeptonMassReject = config -> GetValue(configReject.Data(), true);
  
  // ***********************************************
  // Add these cuts to the cut flow table
  // ***********************************************
  ostringstream cutFlowTitleStream;
  ostringstream cutFlowNameStream;
  TString cutFlowTitle;
  TString cutFlowName;

  TString minReject = _LeptonMassReject?" <= ":" >= ";
  TString maxReject = _LeptonMassReject?" >= ":" <= ";

  // Min cut
  cutFlowTitleStream << leptonType.Data() << " Dilepton" << massCutPostfix << " : Mass" << minReject << _LeptonMassMin;
  cutFlowTitle = cutFlowTitleStream.str().c_str();

  cutFlowNameStream << leptonType.Data() << "Dilepton.Mass.Min" << massCutPostfix;
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName.Data(), cutFlowTitle.Data());

  // Max cut
  cutFlowTitleStream.str("");
  cutFlowTitleStream << leptonType.Data() << " Dilepton" << massCutPostfix << " : Mass" << maxReject << _LeptonMassMax;
  cutFlowTitle = cutFlowTitleStream.str().c_str();
  
  cutFlowNameStream.str("");
  cutFlowNameStream << leptonType.Data() << "Dilepton.Mass.Max" << massCutPostfix;
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName.Data(), cutFlowTitle.Data());

  // Min + Max cut
  cutFlowTitleStream.str("");
  cutFlowTitleStream << leptonType.Data() << " Dilepton" << massCutPostfix << " : Mass" << minReject << _LeptonMassMin << " and " << "Mass" << maxReject << _LeptonMassMax;
  cutFlowTitle = cutFlowTitleStream.str().c_str();

  cutFlowNameStream.str("");
  cutFlowNameStream << leptonType.Data() << "Dilepton.Mass.All" << massCutPostfix;
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName, cutFlowTitle);

  
}//BookHistograms()

/******************************************************************************
 * Bool_t CutLeptonMass::Apply()                                                 *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Cut on the number of electrons and the p_T of these electrons              *
 * Apply separate p_T cuts to the first and second (don't worry about others) *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t CutLeptonMass::Apply()
{

  // ***********************************************
  // Get Number of Leptons and fill histograms
  // ***********************************************
  
  // Initialize number of leptons
  Int_t MuonNumber      = 0;       
  Int_t ElectronNumber  = 0;       
  Float_t LeptonMass    = 0;       

  // Flags 
  Bool_t LeptonMassMinPass = kTRUE;   // Event passes min pT lepton cuts
  Bool_t LeptonMassMaxPass = kTRUE;   // Event passes max pT lepton cuts

  // Get Event Container
  EventContainer *EventContainerObj = GetEventContainer();

  std::vector<Muon> muonVector;
  std::vector<Electron> electronVector;

  //Assign collections to the above defined vectors
  if(     "All"        == leptonType) {
    muonVector.assign(EventContainerObj -> muons.begin(),           EventContainerObj -> muons.end());
    electronVector.assign(EventContainerObj -> electrons.begin(),           EventContainerObj -> electrons.end());
  }
  else if("Veto"       == leptonType) {
    muonVector.assign(EventContainerObj -> vetoMuons.begin(),       EventContainerObj -> vetoMuons.end());
    electronVector.assign(EventContainerObj -> vetoElectrons.begin(),       EventContainerObj -> vetoElectrons.end());
  }
  else if("PtEtaCut"   == leptonType) {
    muonVector.assign(EventContainerObj -> ptetaMuons.begin(),      EventContainerObj -> ptetaMuons.end());
    electronVector.assign(EventContainerObj -> ptetaElectrons.begin(),      EventContainerObj -> ptetaElectrons.end());
  }
  else if("Tight"      == leptonType) {
    muonVector.assign(EventContainerObj -> tightMuons.begin(),      EventContainerObj -> tightMuons.end());
    electronVector.assign(EventContainerObj -> tightElectrons.begin(),      EventContainerObj -> tightElectrons.end());
  }
  else if("Isolated"   == leptonType) {
    muonVector.assign(EventContainerObj -> isolatedMuons.begin(),   EventContainerObj -> isolatedMuons.end());
    electronVector.assign(EventContainerObj -> isolatedElectrons.begin(),   EventContainerObj -> isolatedElectrons.end());
  }
  else if("UnIsolated" == leptonType) {
    muonVector.assign(EventContainerObj -> unIsolatedMuons.begin(), EventContainerObj -> unIsolatedMuons.end());
    electronVector.assign(EventContainerObj -> unIsolatedElectrons.begin(), EventContainerObj -> unIsolatedElectrons.end());
  }
  else{
    std::cout << "ERROR " << "<HistogramminMuon::Apply()> "
	      << "muonType must be All, Tight, Veto, Isolated, or UnIsolated, PtEtaCut" << std::endl;
    exit(8);
  } //else                                                                                                          
  Particle leadingLep;
  Particle subLeadingLep;
  leadingLep.SetPtEtaPhiE(0.,0.,0.,0.);
  subLeadingLep.SetPtEtaPhiE(0.,0.,0.,0.);
  Bool_t leadingIsMuon;
  Bool_t subLeadingIsMuon;

  // Grab the leading and sub leading leptons
  for (auto mu : muonVector){
    if (mu.Pt() > leadingLep.Pt()){
      if (leadingLep.Pt() > subLeadingLep.Pt()){
	subLeadingLep = leadingLep;
	subLeadingIsMuon = leadingIsMuon;
      }
      leadingLep = mu;
      leadingIsMuon = kTRUE;
    }
    else if (mu.Pt() > subLeadingLep.Pt()){
      subLeadingLep = mu;
      subLeadingIsMuon = kTRUE;
    }
  }
  for (auto ele : electronVector){
    if (ele.Pt() > leadingLep.Pt()){
      if (leadingLep.Pt() > subLeadingLep.Pt()){
	subLeadingLep = leadingLep;
	subLeadingIsMuon = leadingIsMuon;
      }
      leadingLep = ele;
      leadingIsMuon = kFALSE;
    }
    else if (ele.Pt() > subLeadingLep.Pt()){
      subLeadingLep = ele;
      subLeadingIsMuon = kFALSE;
    }
  }
  //If the configuration doesn't match that required by the channel, return false
  if (EventContainerObj->GetChannelName() == "mumu" && !(leadingIsMuon && subLeadingIsMuon)) return false;
  if (EventContainerObj->GetChannelName() == "ee" && (leadingIsMuon || subLeadingIsMuon)) return false;
  if (EventContainerObj->GetChannelName() == "emu" && (leadingIsMuon == subLeadingIsMuon)) {
    std::cout << leadingLep.Pt() << " " << subLeadingLep.Pt() << " " << leadingIsMuon << " " << subLeadingIsMuon << std::endl;
    return false;
  }
  LeptonMass = (leadingLep+subLeadingLep).M();

  // Fill the histograms before the cuts
  _hLeptonMassBefore    -> Fill(LeptonMass);
  
  // ***********************************************
  // Fill cut flow table
  // ***********************************************
  
  // Names for Cut Flow Table
  ostringstream cutFlowNameMinStream;
  ostringstream cutFlowNameMaxStream;
  ostringstream cutFlowNameAllStream;
  
  TString cutFlowNameMin;
  TString cutFlowNameMax;
  TString cutFlowNameAll;
  
  cutFlowNameMinStream << leptonType.Data() << "Dilepton.Mass.Min" << massCutPostfix;
  cutFlowNameMin = cutFlowNameMinStream.str().c_str();
  
  cutFlowNameMaxStream << leptonType.Data() << "Dilepton.Mass.Max" << massCutPostfix;
  cutFlowNameMax = cutFlowNameMaxStream.str().c_str();
  
  cutFlowNameAllStream << leptonType.Data() << "Dilepton.Mass.All" << massCutPostfix;
  cutFlowNameAll = cutFlowNameAllStream.str().c_str();
  
  // Cut on Min Number of Leptons
  // 999 value for Min means there is no Min cut
  if( (_LeptonMassMin != 999.) && ((LeptonMass < _LeptonMassMin) != _LeptonMassReject )){
    LeptonMassMinPass = kFALSE;
    GetCutFlowTable()->FailCut(cutFlowNameMin.Data());
  }//if
  else{
    GetCutFlowTable() -> PassCut(cutFlowNameMin.Data());
  } //else
  
  // Cut on Max Number Leptons
  // 999 value for Min means there is no Min cut
  if( (_LeptonMassMax != 999.) && ((LeptonMass > _LeptonMassMax) != _LeptonMassReject)){
    LeptonMassMaxPass = kFALSE;
    GetCutFlowTable()->FailCut(cutFlowNameMax.Data());
  } //if
  else{
    GetCutFlowTable() -> PassCut(cutFlowNameMax.Data());
  } //else
  
  // Cut on Min and Max Number of Leptons
  if((LeptonMassMinPass != LeptonMassMaxPass) == _LeptonMassReject){
    _hLeptonMassAfter -> Fill(LeptonMass);
    GetCutFlowTable() -> PassCut(cutFlowNameAll.Data());
  } //if
  else { 
    GetCutFlowTable() -> FailCut(cutFlowNameAll.Data());
    //    std::cout << leadingLep.Pt() << " " << subLeadingLep.Pt() << " " << LeptonMass << std::endl;
  }

  // ***********************************************
  // Return result of Min and Max Cut
  // ***********************************************
  
  return((LeptonMassMinPass != LeptonMassMaxPass) == _LeptonMassReject);
 
} //Apply












