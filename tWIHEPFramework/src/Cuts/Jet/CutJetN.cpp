/******************************************************************************
 * CutJetN.hpp                                                                *
 *                                                                            *
 * Cuts on number of jets in an event                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutNJet class                                   *
 *    CutJetN()                         -- Parameterized Constructor          *
 *    ~CutJetN()                        -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetN"                  *
 *                                                                            *
 * Private Data Members of CutNJet class                                      *
 *    myTH1F* _hJetNumberBefore         -- Hist of Num of jets before cuts    *
 *    myTH1F* _hJetNumberAfter          -- Hist of Num of jets before cuts    *
 *    Int_t _JetNumberMin;              -- Minimum number of jets to require  *
 *    Int_t _JetNumberMax;              -- Maximum number of jets allowed     *
 *                                                                            *
 * History                                                                    *
 *       7 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      10 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      17 Nov 2006 - RS: add interface to cutFlowTable.                      *
 *      11 Jan 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Cuts/Jet/CutJetN.hpp"
#include <iostream>
#include <sstream>

using namespace std;


/******************************************************************************
 * CutJetN::CutJetN(EventContainer *EventContainerObj)                      *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
CutJetN::CutJetN(EventContainer *EventContainerObj, Int_t nJetsDefault)
{
  SetEventContainer(EventContainerObj);
  _JetNumberMin = nJetsDefault; //Allows for an over-riding of the default nJets option so we can not have hundreds of config files
  _JetNumberMax = nJetsDefault;
} // CutJetN


/******************************************************************************
 * CutJetN::~CutJetN()                                                      *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
CutJetN::~CutJetN()
{
  
}//~CutJetN

/******************************************************************************
 * void CutJetN::BookHistogram()                                             *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutJetN::BookHistogram(){
  
  // Histogram of events before cut on number of jets
  _hJetNumberBefore =  DeclareTH1F("JetNumberBefore","Number of Jets before the cut",20, 0.0, 20.0);
  _hJetNumberBefore -> SetXAxisTitle("N_{Jet}");
  _hJetNumberBefore -> SetYAxisTitle("Events");

  // Histogram of events after cut on number of jets
  _hJetNumberAfter =  DeclareTH1F("JetNumberAfter","Number of Jets after the cut",20, 0.0, 20.0);
  _hJetNumberAfter -> SetXAxisTitle("N_{Jet}");
  _hJetNumberAfter -> SetYAxisTitle("Events");
  
  // Get configuration
  EventContainer *EventContainerObj = GetEventContainer();
  TEnv *config = EventContainerObj -> GetConfig();

  // Use configuration to set min and max number of jets to cut on, if we haven't over-ridden this from the command line
  if (_JetNumberMax < 0){
    _JetNumberMin = config -> GetValue("Cut.Jet.Number.Min", 999);
    _JetNumberMax = config -> GetValue("Cut.Jet.Number.Max", 999);
  }
  //
  // also add these two cuts to the cut flow table
  ostringstream titleStr;
  titleStr <<"Jet : " << "N >=" << _JetNumberMin;
  TString tTitle(titleStr.str().c_str());
  GetCutFlowTable()->AddCutToFlow("Jet.Number.Min",tTitle);
  titleStr.str("");
  titleStr <<"Jet : N <= "<<_JetNumberMax;
  tTitle=titleStr.str().c_str();
  GetCutFlowTable()->AddCutToFlow("Jet.Number.Max",tTitle);
  // also add the overal pass rate
  titleStr.str("");
  titleStr <<"Jet : "<< _JetNumberMin<<" <= N <= "<<_JetNumberMax;
  tTitle=titleStr.str().c_str();
  GetCutFlowTable()->AddCutToFlow("Jet.Number.All",tTitle);

  
}//BookHistograms()

/******************************************************************************
 * Bool_t CutJetN::Apply()                                                   *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t CutJetN::Apply()
{
  //cout<<"Begin of CutJetN::Apply()"<<endl;

  // Get Event Tree
  EventContainer *evObj = GetEventContainer();

  Int_t JetNumber = evObj->jets.size();
  // Fill the histogram before the cut
  _hJetNumberBefore -> Fill(JetNumber);

  // Flags 
  Bool_t JetNumberMinPass    = kTRUE;   // Event passes min cut
  Bool_t JetNumberMaxPass    = kTRUE;   // Event passes max cut

  // Cut on Min Pt of Jet
  // Negative cut value for Min means there is no Min cut
  if( (_JetNumberMin != 999) && (JetNumber < _JetNumberMin) ){
     JetNumberMinPass = kFALSE;
     GetCutFlowTable()->FailCut("Jet.Number.Min");
  }//if
  else{
    GetCutFlowTable() -> PassCut("Jet.Number.Min");
  } //else

  // Cut on Min Pt of Jet
  // Negative cut value for Min means there is no Min cut
  if( (_JetNumberMax != 999) && (JetNumber > _JetNumberMax) ){
     JetNumberMaxPass = kFALSE;
     GetCutFlowTable()->FailCut("Jet.Number.Max");
  } //if
  else{
    GetCutFlowTable() -> PassCut("Jet.Number.Max");
  } //else

  // Fill Histogram after cuts
  if(JetNumberMinPass && JetNumberMaxPass){
    _hJetNumberAfter -> Fill(JetNumber);
    GetCutFlowTable() -> PassCut("Jet.Number.All");
  } //if
  else GetCutFlowTable() -> FailCut("Jet.Number.All");

  return(JetNumberMinPass && JetNumberMaxPass);

} //Apply












