/******************************************************************************
 * CutElectronN.cpp                                                           *
 *                                                                            *
 * Cuts on electron Number (can cut on Tight or Veto Electrons)               *
 * Must pass Tight or Veto to constructor                                     *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutElectronN class                              *
 *    CutElectronN()                     -- Parameterized Constructor         *
 *    ~CutElectronN()                    -- Destructor                        *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutElectronN"                *
 *                                                                            *
 * Private Data Members of CutElectronN class                                 *
 *    myTH1F* _hElectronNumberBefore;    -- Hist el Number before cut         *
 *    myTH1F* _hElectronNumberAfter;     -- Hist el Number of jets after cut  *
 *                                                                            *
 *    Int_t _ElectronNumberMin;          -- Minimum Electron Number           *
 *    Int_t _ElectronNumberMax;          -- Maximum Electron Number           *
 *                                                                            *
 * History                                                                    *
 *      15 Jan 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Cuts/Electron/CutElectronN.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * CutElectronN::CutElectronN(EventContainer *EventContainerObj, TString electronType)    *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
CutElectronN::CutElectronN(EventContainer *EventContainerObj, TString electronTypePassed)
{
  // Check electronType parameter
  if( electronTypePassed.CompareTo("All") && electronTypePassed.CompareTo("UnIsolated") && electronTypePassed.CompareTo("Isolated") && 
      electronTypePassed.CompareTo("Tight") && electronTypePassed.CompareTo("PtEtaCut") && electronTypePassed.CompareTo("Veto") ){
    std::cout << "ERROR " << "<CutElectronN::CutElectronN()> " 
	      << "Must pass Tight, PtEtaCut, Veto, All, UnIsolated, or Isolated to constructor" << std::endl;
    exit(8);
  } //if
  electronType = electronTypePassed;

  // Set Event Container
  SetEventContainer(EventContainerObj);
} // CutElectronN


/******************************************************************************
 * CutElectronN::~CutElectronN()                                              *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
CutElectronN::~CutElectronN()
{
  
}//~CutElectronN

/******************************************************************************
 * void CutElectronN::BookHistogram()                                         *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutElectronN::BookHistogram(){
  
  // ***********************************************
  // Make Strings for histogram titles and labels
  // ***********************************************  

  // Histogram Before Cut
  std::ostringstream histNameBeforeStream;
  histNameBeforeStream << electronType << "ElectronNumberBefore";
  TString histNameBefore = histNameBeforeStream.str().c_str();

  std::ostringstream histTitleBeforeStream;
  histTitleBeforeStream << electronType << "Electron Number Before Cut";
  TString histTitleBefore = histTitleBeforeStream.str().c_str();

  // Histogram After Cut
  std::ostringstream histNameAfterStream;
  histNameAfterStream << electronType << "ElectronNumberAfter";
  TString histNameAfter = histNameAfterStream.str().c_str();

  std::ostringstream histTitleAfterStream;
  histTitleAfterStream << electronType << "Electron Number After Cut";
  TString histTitleAfter = histTitleAfterStream.str().c_str();

  // ***********************************************
  // Book Histograms
  // ***********************************************  

  // Histogram before cut
  _hElectronNumberBefore =  DeclareTH1F(histNameBefore.Data(), histTitleBefore.Data(), 10, 0.0, 10.0);
  _hElectronNumberBefore -> SetXAxisTitle("N_{el}");
  _hElectronNumberBefore -> SetYAxisTitle("Events");

  // Histogram after cut
  _hElectronNumberAfter =  DeclareTH1F(histNameAfter.Data(), histTitleAfter.Data(), 10, 0.0, 10);
  _hElectronNumberAfter -> SetXAxisTitle("N_{el}");
  _hElectronNumberAfter -> SetYAxisTitle("Events");

  // ***********************************************
  // Get cuts from configuration file
  // ***********************************************  

  // Get configuration file
  EventContainer *EventContainerObj = GetEventContainer();
  TEnv *config = EventContainerObj -> GetConfig();

  // Make name in configuration file depending upon electron type (veto, tight)
  std::ostringstream configMinStream;
  configMinStream << "Cut.Electron." << electronType.Data() << ".Number.Min";
  TString configMin = configMinStream.str().c_str();

  std::ostringstream configMaxStream;
  configMaxStream << "Cut.Electron." << electronType.Data() << ".Number.Max";
  TString configMax = configMaxStream.str().c_str();

  // Use configuration to set min and max number of jets to cut on
  _ElectronNumberMin  = config -> GetValue(configMin.Data(), 999);
  _ElectronNumberMax  = config -> GetValue(configMax.Data(), 999);
  
  // ***********************************************
  // Add these cuts to the cut flow table
  // ***********************************************
  ostringstream cutFlowTitleStream;
  ostringstream cutFlowNameStream;
  TString cutFlowTitle;
  TString cutFlowName;

  // Min cut
  cutFlowTitleStream << electronType.Data() << " Electron : " << "N >= " << _ElectronNumberMin;
  cutFlowTitle = cutFlowTitleStream.str().c_str();

  cutFlowNameStream << electronType.Data() << "Electron.Number.Min";
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName.Data(), cutFlowTitle.Data());

  // Max cut
  cutFlowTitleStream.str("");
  cutFlowTitleStream << electronType.Data() << " Electron : " << "N <= " << _ElectronNumberMax;
  cutFlowTitle = cutFlowTitleStream.str().c_str();
  
  cutFlowNameStream.str("");
  cutFlowNameStream << electronType.Data() << "Electron.Number.Max";
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName.Data(), cutFlowTitle.Data());

  // Min + Max cut
  cutFlowTitleStream.str("");
  cutFlowTitleStream << electronType.Data() << " Electron : " << _ElectronNumberMin << " <= N <= " << _ElectronNumberMax;
  cutFlowTitle = cutFlowTitleStream.str().c_str();

  cutFlowNameStream.str("");
  cutFlowNameStream << electronType.Data() << "Electron.Number.All";
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName, cutFlowTitle);

  
}//BookHistograms()

/******************************************************************************
 * Bool_t CutElectronN::Apply()                                               *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Cut on the number of electrons and the p_T of these electrons              *
 * Apply separate p_T cuts to the first and second (don't worry about others) *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t CutElectronN::Apply()
{

  // ***********************************************
  // Get Number of Electrons and fill histograms
  // ***********************************************
  
  // Initialize number of electrons
  Int_t ElectronNumber  = 0;       

  // Flags 
  Bool_t ElectronNumberMinPass = kTRUE;   // Event passes min pT electron cuts
  Bool_t ElectronNumberMaxPass = kTRUE;   // Event passes max pT electron cuts

  // Get Event Container
  EventContainer *EventContainerObj = GetEventContainer();

  // Get variables from EventContainerObj - Depens on electron Type
  if(      ! electronType.CompareTo("Veto") )       ElectronNumber = EventContainerObj -> vetoElectrons.size();
  else if( ! electronType.CompareTo("Tight") )      ElectronNumber = EventContainerObj -> tightElectrons.size();
  else if( ! electronType.CompareTo("PtEtaCut") )   ElectronNumber = EventContainerObj -> ptetaElectrons.size();
  else if( ! electronType.CompareTo("UnIsolated") ) ElectronNumber = EventContainerObj -> unIsolatedElectrons.size();
  else if( ! electronType.CompareTo("Isolated") )   ElectronNumber = EventContainerObj -> isolatedElectrons.size();
  else if( ! electronType.CompareTo("All") )        ElectronNumber = EventContainerObj -> electrons.size();
  else{
    std::cout << "ERROR " << "<CutElectronN::Apply()> " << "electronType has an incorrect value of: " << electronType.Data() << std::endl;
    exit(8);
  } //else

  // Fill the histograms before the cuts
  _hElectronNumberBefore    -> Fill(ElectronNumber);
  //  std::cout << electronType << " " << ElectronNumber << std::endl;
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
  
  cutFlowNameMinStream << electronType.Data() << "Electron.Number.Min";
  cutFlowNameMin = cutFlowNameMinStream.str().c_str();
  
  cutFlowNameMaxStream << electronType.Data() << "Electron.Number.Max";
  cutFlowNameMax = cutFlowNameMaxStream.str().c_str();
  
  cutFlowNameAllStream << electronType.Data() << "Electron.Number.All";
  cutFlowNameAll = cutFlowNameAllStream.str().c_str();
  
  // Cut on Min Number of Electrons
  // 999 value for Min means there is no Min cut
  if( (_ElectronNumberMin != 999) && (ElectronNumber < _ElectronNumberMin) ){
    //    std::cout << "FalseMin " << ElectronNumber << " < " << _ElectronNumberMin << std::endl;
    ElectronNumberMinPass = kFALSE;
    GetCutFlowTable()->FailCut(cutFlowNameMin.Data());
  }//if
  else{
    GetCutFlowTable() -> PassCut(cutFlowNameMin.Data());
  } //else
  
  // Cut on Max Number Electrons
  // 999 value for Max means there is no Max cut
  if( (_ElectronNumberMax != 999) && (ElectronNumber > _ElectronNumberMax) ){
    //    std::cout << "FalseMax " << ElectronNumber << " > " << _ElectronNumberMax << std::endl;
    ElectronNumberMaxPass = kFALSE;
    GetCutFlowTable()->FailCut(cutFlowNameMax.Data());
  } //if
  else{
    GetCutFlowTable() -> PassCut(cutFlowNameMax.Data());
  } //else
  
  // Cut on Min and Max Number of Electrons
  if(ElectronNumberMinPass && ElectronNumberMaxPass){
    _hElectronNumberAfter -> Fill(ElectronNumber);
    GetCutFlowTable() -> PassCut(cutFlowNameAll.Data());
  } //if
  else GetCutFlowTable() -> FailCut(cutFlowNameAll.Data());

  // ***********************************************
  // Return result of Min and Max Cut
  // ***********************************************

  return(ElectronNumberMinPass && ElectronNumberMaxPass);
 
} //Apply












