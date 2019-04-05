/******************************************************************************
 * CutMuonN.hpp                                                               *
 *                                                                            *
 * Cuts on muon Number (can cut on All, UnIsolated, Tight, or Veto Muons)     *
 * Must pass All, UnIsolated, Tight, or Veto to constructor                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutMuonN class                                  *
 *    CutMuonN()                     -- Parameterized Constructor             *
 *    ~CutMuonN()                    -- Destructor                            *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutMuonN"                    *
 *                                                                            *
 * Private Data Members of CutMuonN class                                     *
 *    myTH1F* _hMuonNumberBefore;    -- Hist mu Number before cut             *
 *    myTH1F* _hMuonNumberAfter;     -- Hist mu Number of jets after cut      *
 *                                                                            *
 *    Int_t _MuonNumberMin;          -- Minimum Muon Number                   *
 *    Int_t _MuonNumberMax;          -- Maximum Muon Number                   *
 *                                                                            *
 * History                                                                    *
 *      15 Jan 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Cuts/Muon/CutMuonN.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * CutMuonN::CutMuonN(EventContainer *EventContainerObj, TString muonType)    *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
CutMuonN::CutMuonN(EventContainer *EventContainerObj, TString muonTypePassed)
{
  // Check muonType parameter
  if( muonTypePassed.CompareTo("All") && muonTypePassed.CompareTo("UnIsolated") && muonTypePassed.CompareTo("Isolated") && 
      muonTypePassed.CompareTo("Tight") && muonTypePassed.CompareTo("PtEtaCut") && muonTypePassed.CompareTo("Veto") ){
    std::cout << "ERROR " << "<CutMuonN::CutMuonN()> " 
	      << "Must pass All, Tight, PtEtaCut, Veto, Isolated, or UnIsolated to constructor" << std::endl;
    exit(8);
  } //if
  muonType = muonTypePassed;

  // Set Event Container
  SetEventContainer(EventContainerObj);
} // CutMuonN


/******************************************************************************
 * CutMuonN::~CutMuonN()                                                *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
CutMuonN::~CutMuonN()
{
  
}//~CutMuonN

/******************************************************************************
 * void CutMuonN::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutMuonN::BookHistogram(){
  
  // ***********************************************
  // Make Strings for histogram titles and labels
  // ***********************************************  

  // Histogram Before Cut
  std::ostringstream histNameBeforeStream;
  histNameBeforeStream << muonType << "MuonNumberBefore";
  TString histNameBefore = histNameBeforeStream.str().c_str();

  std::ostringstream histTitleBeforeStream;
  histTitleBeforeStream << muonType << "Muon Number Before Cut";
  TString histTitleBefore = histTitleBeforeStream.str().c_str();

  // Histogram After Cut
  std::ostringstream histNameAfterStream;
  histNameAfterStream << muonType << "MuonNumberAfter";
  TString histNameAfter = histNameAfterStream.str().c_str();

  std::ostringstream histTitleAfterStream;
  histTitleAfterStream << muonType << "Muon Number After Cut";
  TString histTitleAfter = histTitleAfterStream.str().c_str();

  // ***********************************************
  // Book Histograms
  // ***********************************************  

  // Histogram before cut
  _hMuonNumberBefore =  DeclareTH1F(histNameBefore.Data(), histTitleBefore.Data(), 10, 0.0, 10.0);
  _hMuonNumberBefore -> SetXAxisTitle("N_{#mu}");
  _hMuonNumberBefore -> SetYAxisTitle("Events");

  // Histogram after cut
  _hMuonNumberAfter =  DeclareTH1F(histNameAfter.Data(), histTitleAfter.Data(), 10, 0.0, 10);
  _hMuonNumberAfter -> SetXAxisTitle("N_{#mu}");
  _hMuonNumberAfter -> SetYAxisTitle("Events");

  // ***********************************************
  // Get cuts from configuration file
  // ***********************************************  

  // Get configuration file
  EventContainer *EventContainerObj = GetEventContainer();
  TEnv *config = EventContainerObj -> GetConfig();

  // Make name in configuration file depending upon muon type (all, veto, tight, unIsolated)
  std::ostringstream configMinStream;
  configMinStream << "Cut.Muon." << muonType.Data() << ".Number.Min";
  TString configMin = configMinStream.str().c_str();

  std::ostringstream configMaxStream;
  configMaxStream << "Cut.Muon." << muonType.Data() << ".Number.Max";
  TString configMax = configMaxStream.str().c_str();

  // Use configuration to set min and max number of jets to cut on
  _MuonNumberMin  = config -> GetValue(configMin.Data(), 999);
  _MuonNumberMax  = config -> GetValue(configMax.Data(), 999);
  
  // ***********************************************
  // Add these cuts to the cut flow table
  // ***********************************************
  ostringstream cutFlowTitleStream;
  ostringstream cutFlowNameStream;
  TString cutFlowTitle;
  TString cutFlowName;

  // Min cut
  cutFlowTitleStream << muonType.Data() << " Muon : " << "N >= " << _MuonNumberMin;
  cutFlowTitle = cutFlowTitleStream.str().c_str();

  cutFlowNameStream << muonType.Data() << "Muon.Number.Min";
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName.Data(), cutFlowTitle.Data());

  // Max cut
  cutFlowTitleStream.str("");
  cutFlowTitleStream << muonType.Data() << " Muon : " << "N <= " << _MuonNumberMax;
  cutFlowTitle = cutFlowTitleStream.str().c_str();
  
  cutFlowNameStream.str("");
  cutFlowNameStream << muonType.Data() << "Muon.Number.Max";
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName.Data(), cutFlowTitle.Data());

  // Min + Max cut
  cutFlowTitleStream.str("");
  cutFlowTitleStream << muonType.Data() << " Muon : " << _MuonNumberMin << " <= N <= " << _MuonNumberMax;
  cutFlowTitle = cutFlowTitleStream.str().c_str();

  cutFlowNameStream.str("");
  cutFlowNameStream << muonType.Data() << "Muon.Number.All";
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName, cutFlowTitle);

  
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
Bool_t CutMuonN::Apply()
{

  // ***********************************************
  // Get Number of Muons and fill histograms
  // ***********************************************
  
  // Initialize number of muons
  Int_t MuonNumber  = 0;       

  // Flags 
  Bool_t MuonNumberMinPass = kTRUE;   // Event passes min pT muon cuts
  Bool_t MuonNumberMaxPass = kTRUE;   // Event passes max pT muon cuts

  // Get Event Container
  EventContainer *EventContainerObj = GetEventContainer();

  // Get variables from EventContainerObj - Depens on muon Type
  if(      ! muonType.CompareTo("Veto") )       MuonNumber = EventContainerObj -> vetoMuons.size();
  else if( ! muonType.CompareTo("Tight") )      MuonNumber = EventContainerObj -> tightMuons.size();
  else if( ! muonType.CompareTo("PtEtaCut") )   MuonNumber = EventContainerObj -> ptetaMuons.size();
  else if( ! muonType.CompareTo("UnIsolated") ) MuonNumber = EventContainerObj -> unIsolatedMuons.size();
  else if( ! muonType.CompareTo("Isolated") )   MuonNumber = EventContainerObj -> isolatedMuons.size();
  else if( ! muonType.CompareTo("All") )        MuonNumber = EventContainerObj -> muons.size();
  else{
    std::cout << "ERROR " << "<CutMuonN::Apply()> " << "muonType has an incorrect value of: " << muonType.Data() << std::endl;
    exit(8);
  } //else

  // Fill the histograms before the cuts
  _hMuonNumberBefore    -> Fill(MuonNumber);
  
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
  
  cutFlowNameMinStream << muonType.Data() << "Muon.Number.Min";
  cutFlowNameMin = cutFlowNameMinStream.str().c_str();
  
  cutFlowNameMaxStream << muonType.Data() << "Muon.Number.Max";
  cutFlowNameMax = cutFlowNameMaxStream.str().c_str();
  
  cutFlowNameAllStream << muonType.Data() << "Muon.Number.All";
  cutFlowNameAll = cutFlowNameAllStream.str().c_str();
  
  // Cut on Min Number of Muons
  // 999 value for Min means there is no Min cut
  if( (_MuonNumberMin != 999) && (MuonNumber < _MuonNumberMin) ){
    MuonNumberMinPass = kFALSE;
    GetCutFlowTable()->FailCut(cutFlowNameMin.Data());
  }//if
  else{
    GetCutFlowTable() -> PassCut(cutFlowNameMin.Data());
  } //else
  
  // Cut on Max Number Muons
  // 999 value for Min means there is no Min cut
  if( (_MuonNumberMax != 999) && (MuonNumber > _MuonNumberMax) ){
    MuonNumberMaxPass = kFALSE;
    GetCutFlowTable()->FailCut(cutFlowNameMax.Data());
  } //if
  else{
    GetCutFlowTable() -> PassCut(cutFlowNameMax.Data());
  } //else
  
  // Cut on Min and Max Number of Muons
  if(MuonNumberMinPass && MuonNumberMaxPass){
    _hMuonNumberAfter -> Fill(MuonNumber);
    GetCutFlowTable() -> PassCut(cutFlowNameAll.Data());
  } //if
  else GetCutFlowTable() -> FailCut(cutFlowNameAll.Data());

  // ***********************************************
  // Return result of Min and Max Cut
  // ***********************************************
  //  cout<<"MuonNCut: Event: "<<EventContainerObj ->eventNumber<<endl;
  // if(MuonNumberMinPass && MuonNumberMaxPass) cout<<EventContainerObj ->eventNumber<<endl;
  return(MuonNumberMinPass && MuonNumberMaxPass);
 
} //Apply












