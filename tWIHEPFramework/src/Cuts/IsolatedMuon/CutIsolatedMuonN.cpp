/******************************************************************************
 * CutIsolatedMuonN.hpp                                                       *
 *                                                                            *
 * Cuts on number of isolated muons in an event                               *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutIsolatedMuonN class                          *
 *    CutIsolatedMuonN()                   -- Parameterized Constructor       *
 *    ~CutIsolatedMuonN()                  -- Destructor                      *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutIsolatedMuonN"         *
 *                                                                            *
 * Private Data Members of CutIsolatedMuonN                                   *
 *    myTH1F* _hMuonNumberBefore         -- Hist of Num of muons before cuts  *
 *    myTH1F* _hMuonNumberAfter          -- Hist of Num of muons before cuts  *
 *    Int_t _MuonNumberMin;              -- Minimum number of muons to require*
 *    Int_t _MuonNumberMax;              -- Maximum number of muons allowed   *
 *                                                                            *
 * History                                                                    *
 *      17 Jul 2007 - created by J. Holzbauer based on CutMuonN               *
 *****************************************************************************/


#include "SingleTopRootAnalysis/Cuts/IsolatedMuon/CutIsolatedMuonN.hpp"
#include <iostream>
#include <sstream>

using namespace std;


/******************************************************************************
 * CutIsolatedMuonN::CutIsolatedMuonN(EventContainer *EventContainerObj)      *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
CutIsolatedMuonN::CutIsolatedMuonN(EventContainer *EventContainerObj)
{
  SetEventContainer(EventContainerObj);
} // CutIsolatedMuonN


/******************************************************************************
 * CutIsolatedMuonN::~CutIsolatedMuonN()                                      *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
CutIsolatedMuonN::~CutIsolatedMuonN()
{
  
}//~CutIsolatedMuonN

/******************************************************************************
 * void CutIsolatedMuonN::BookHistogram()                                     *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutIsolatedMuonN::BookHistogram(){
  
  // Histogram of events before cut on number of muons
  _hMuonNumberBefore =  DeclareTH1F("IsolatedMuonNumberBefore","Number of Isolated Muons before the cut",10, 0.0, 10.0);
  _hMuonNumberBefore -> SetXAxisTitle("N_{Muon}^{Isolated}");
  _hMuonNumberBefore -> SetYAxisTitle("Events");

  // Histogram of events after cut on number of muons
  _hMuonNumberAfter =  DeclareTH1F("IsolatedMuonNumberAfter","Number of Isolated Muons after the cut",10, 0.0, 10.0);
  _hMuonNumberAfter -> SetXAxisTitle("N_{Muon}^{Isolated}");
  _hMuonNumberAfter -> SetYAxisTitle("Events");
  
  // Get configuration
  EventContainer *EventContainerObj = GetEventContainer();
  TEnv *config = EventContainerObj -> GetConfig();

  // Use configuration to set min and max number of muons to cut on
  _MuonNumberMin = config -> GetValue("Cut.Isolated.Muon.Number.Min", 999);
  _MuonNumberMax = config -> GetValue("Cut.Isolated.Muon.Number.Max", 999);

  //
  // also add these two cuts to the cut flow table
  ostringstream titleStr;
  titleStr <<"Isolated Muon : " << "N >=" << _MuonNumberMin;
  TString tTitle(titleStr.str().c_str());
  GetCutFlowTable()->AddCutToFlow("Isolated.Muon.Number.Min",tTitle);
  titleStr.str("");
  titleStr <<"Isolated Muon : N <= "<<_MuonNumberMax;
  tTitle=titleStr.str().c_str();
  GetCutFlowTable()->AddCutToFlow("Isolated.Muon.Number.Max",tTitle);
  // also add the overal pass rate
  titleStr.str("");
  titleStr <<"Isolated Muon : "<< _MuonNumberMin<<" <= N <= "<<_MuonNumberMax;
  tTitle=titleStr.str().c_str();
  GetCutFlowTable()->AddCutToFlow("Isolated.Muon.Number.All",tTitle);
  
}//BookHistograms()

/******************************************************************************
 * Bool_t CutIsolatedMuonN::Apply()                                           *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t CutIsolatedMuonN::Apply()
{
  //cout<<"Begin of CutIsolatedMuonN::Apply()"<<endl;

  // Get Event Tree
  EventContainer *evObj = GetEventContainer();

  Int_t MuonNumber = evObj->tightMuons.size();

  // Fill the histogram before the cut
  _hMuonNumberBefore -> Fill(MuonNumber);

  // Flags 
  Bool_t MuonNumberMinPass    = kTRUE;   // Event passes min cut
  Bool_t MuonNumberMaxPass    = kTRUE;   // Event passes max cut

  // Cut on Min Pt of Muon
  // Negative cut value for Min means there is no Min cut
  if( (_MuonNumberMin != 999) && (MuonNumber < _MuonNumberMin) ){
     MuonNumberMinPass = kFALSE;
     GetCutFlowTable()->FailCut("Isolated.Muon.Number.Min");
  }//if
  else{
    GetCutFlowTable() -> PassCut("Isolated.Muon.Number.Min");
  } //else

  // Cut on Min Pt of Muon
  // Negative cut value for Min means there is no Min cut
  if( (_MuonNumberMax != 999) && (MuonNumber > _MuonNumberMax) ){
     MuonNumberMaxPass = kFALSE;
     GetCutFlowTable()->FailCut("Isolated.Muon.Number.Max");
  } //if
  else{
    GetCutFlowTable() -> PassCut("Isolated.Muon.Number.Max");
  } //else

  // Fill Histogram after cuts
  if(MuonNumberMinPass && MuonNumberMaxPass){
    _hMuonNumberAfter -> Fill(MuonNumber);
    GetCutFlowTable() -> PassCut("Isolated.Muon.Number.All");
  } //if
  else GetCutFlowTable() -> FailCut("Isolated.Muon.Number.All");

  ///////////////  TMP ///////////////////
  EventTree *evTree = evObj -> GetEventTree();
  if(MuonNumberMinPass && MuonNumberMaxPass){
    //  std::cout << "AAA" << evTree -> eventNumber << " " << evTree -> eventWeightMCatNLO << std::endl;
  } //if
  ///////////////  TMP ///////////////////

  return(MuonNumberMinPass && MuonNumberMaxPass);

} //Apply











