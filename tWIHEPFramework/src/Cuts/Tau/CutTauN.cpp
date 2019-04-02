/******************************************************************************
 * CutTauN.hpp                                                                *
 *                                                                            *
 * Cuts on tau Number                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTauN class                                   *
 *    CutTauN()                     -- Parameterized Constructor              *
 *    ~CutTauN()                    -- Destructor                             *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutTauN"                     *
 *                                                                            *
 * Private Data Members of CutTauN class                                      *
 *    myTH1F* _hTauNumberBefore;    -- Hist mu Number before cut              *
 *    myTH1F* _hTauNumberAfter;     -- Hist mu Number of jets after cut       *
 *                                                                            *
 *    Int_t _TauNumberMin;          -- Minimum Tau Number                     *
 *    Int_t _TauNumberMax;          -- Maximum Tau Number                     *
 *                                                                            *
 * History                                                                    *
 *      22 Aug 2007 - Created by P. Ryan                                      *
 *****************************************************************************/


#include "SingleTopRootAnalysis/Cuts/Tau/CutTauN.hpp"
#include<iostream>

using namespace std;


/******************************************************************************
 * CutTauN::CutTauN(EventContainer *EventContainerObj)                        *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
CutTauN::CutTauN(EventContainer *EventContainerObj)
{
  SetEventContainer(EventContainerObj);
} // CutTauN


/******************************************************************************
 * CutTauN::~CutTauN()                                                        *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
CutTauN::~CutTauN()
{
  
}//~CutTauN

/******************************************************************************
 * void CutTauN::BookHistogram()                                              *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutTauN::BookHistogram(){
  
  // Histogram of events before cut on number of jets
  // Number of Electrons
  _hTauNumberBefore =  DeclareTH1F("TauNumberBefore","N_{#mu} before cut",10, 0, 10.0);
  _hTauNumberBefore -> SetXAxisTitle("N_{#mu}");
  _hTauNumberBefore -> SetYAxisTitle("Events");

  // Histogram of events after cut on number of electrons
  // Number of Electrons
  _hTauNumberAfter =  DeclareTH1F("TauNumberAfter","N_{#mu} after cut",10, 0, 10);
  _hTauNumberAfter -> SetXAxisTitle("N_{#mu}");
  _hTauNumberAfter -> SetYAxisTitle("Events");

  // Get configuration
  EventContainer *EventContainerObj = GetEventContainer();
  TEnv *config = EventContainerObj -> GetConfig();

  // Use configuration to set min and max number of jets to cut on
  _TauNumberMin  = config -> GetValue("Cut.Tau.Number.Min", 999);
  _TauNumberMax  = config -> GetValue("Cut.Tau.Number.Max", 999);
  
  // Add these cuts to the cut flow table
  ostringstream titleStr;
  titleStr << "Tau : " << "N >= " << _TauNumberMin;
  TString tTitle(titleStr.str().c_str());
  GetCutFlowTable()->AddCutToFlow("Tau.Number.Min", tTitle);

  titleStr.str("");
  titleStr << "Tau : " << "N <= " << _TauNumberMax;
  tTitle = titleStr.str().c_str();
  GetCutFlowTable()->AddCutToFlow("Tau.Number.Max", tTitle);

  titleStr.str("");
  titleStr << "Tau : " << _TauNumberMin << " <= N <= " << _TauNumberMax;
  tTitle = titleStr.str().c_str();
  GetCutFlowTable()->AddCutToFlow("Tau.Number.All", tTitle);

  
}//BookHistograms()

/******************************************************************************
 * Bool_t CutTauN::Apply()                                                    *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Cut on the number of electrons and the p_T of these electrons              *
 * Apply separate p_T cuts to the first and second (don't worry about others) *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t CutTauN::Apply()
{

  // Electron Variables on which to cut
  Int_t TauNumber  = 0;         // Number of electrons

  // Flags 
  Bool_t TauNumberMinPass    = kTRUE;   // Event passes min pT electron cuts
  Bool_t TauNumberMaxPass    = kTRUE;   // Event passes max pT electron cuts

  // Get Event Container
  EventContainer *EventContainerObj = GetEventContainer();

  // Get Variables from electron vector in Particles Object
  TauNumber = EventContainerObj -> taus.size();

  // Fill the histograms before the cuts
  _hTauNumberBefore    -> Fill(TauNumber);

  // Cut on Min Pt of Electron
  // Negative cut value for Min means there is no Min cut
  if( (_TauNumberMin != 999) && (TauNumber < _TauNumberMin) ){
     TauNumberMinPass = kFALSE;
     GetCutFlowTable()->FailCut("Tau.Number.Min");
  }//if
  else{
    GetCutFlowTable() -> PassCut("Tau.Number.Min");
  } //else

  // Cut on Min Pt of Electron
  // Negative cut value for Min means there is no Min cut
  if( (_TauNumberMax != 999) && (TauNumber > _TauNumberMax) ){
     TauNumberMaxPass = kFALSE;
     GetCutFlowTable()->FailCut("Tau.Number.Max");
  } //if
  else{
    GetCutFlowTable() -> PassCut("Tau.Number.Max");
  } //else

  // Fill Histogram after cuts
  if(TauNumberMinPass && TauNumberMaxPass){
    _hTauNumberAfter -> Fill(TauNumber);
    GetCutFlowTable() -> PassCut("Tau.Number.All");
  } //if
  else GetCutFlowTable() -> FailCut("Tau.Number.All");

  return(TauNumberMinPass && TauNumberMaxPass);


} //Apply












