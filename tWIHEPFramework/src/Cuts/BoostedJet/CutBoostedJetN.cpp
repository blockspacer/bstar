/******************************************************************************
 * CutBoostedJetN.hpp                                                                *
 *                                                                            *
 * Cuts on number of jets in an event                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutNBoostedJet class                                   *
 *    CutBoostedJetN()                         -- Parameterized Constructor          *
 *    ~CutBoostedJetN()                        -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutBoostedJetN"                  *
 *                                                                            *
 * Private Data Members of CutNBoostedJet class                                      *
 *    myTH1F* _hBoostedJetNumberBefore         -- Hist of Num of jets before cuts    *
 *    myTH1F* _hBoostedJetNumberAfter          -- Hist of Num of jets before cuts    *
 *    Int_t _BoostedJetNumberMin;              -- Minimum number of jets to require  *
 *    Int_t _BoostedJetNumberMax;              -- Maximum number of jets allowed     *
 *                                                                            *
 * History                                                                    *
 *       7 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      10 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      17 Nov 2006 - RS: add interface to cutFlowTable.                      *
 *      11 Jan 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Cuts/BoostedJet/CutBoostedJetN.hpp"
#include <iostream>
#include <sstream>

using namespace std;


/******************************************************************************
 * CutBoostedJetN::CutBoostedJetN(EventContainer *EventContainerObj)                      *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
CutBoostedJetN::CutBoostedJetN(EventContainer *EventContainerObj, Int_t nBoostedJetsDefault)
{
  SetEventContainer(EventContainerObj);
  _BoostedJetNumberMin = nBoostedJetsDefault; //Allows for an over-riding of the default nBoostedJets option so we can not have hundreds of config files
  _BoostedJetNumberMax = nBoostedJetsDefault;
} // CutBoostedJetN


/******************************************************************************
 * CutBoostedJetN::~CutBoostedJetN()                                                      *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
CutBoostedJetN::~CutBoostedJetN()
{
  
}//~CutBoostedJetN

/******************************************************************************
 * void CutBoostedJetN::BookHistogram()                                             *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutBoostedJetN::BookHistogram(){
  
  // Histogram of events before cut on number of jets
  _hBoostedJetNumberBefore =  DeclareTH1F("BoostedJetNumberBefore","Number of BoostedJets before the cut",20, 0.0, 20.0);
  _hBoostedJetNumberBefore -> SetXAxisTitle("N_{BoostedJet}");
  _hBoostedJetNumberBefore -> SetYAxisTitle("Events");

  // Histogram of events after cut on number of jets
  _hBoostedJetNumberAfter =  DeclareTH1F("BoostedJetNumberAfter","Number of BoostedJets after the cut",20, 0.0, 20.0);
  _hBoostedJetNumberAfter -> SetXAxisTitle("N_{BoostedJet}");
  _hBoostedJetNumberAfter -> SetYAxisTitle("Events");
  
  // Get configuration
  EventContainer *EventContainerObj = GetEventContainer();
  TEnv *config = EventContainerObj -> GetConfig();

  // Use configuration to set min and max number of jets to cut on, if we haven't over-ridden this from the command line
  if (_BoostedJetNumberMax < 0){
    _BoostedJetNumberMin = config -> GetValue("Cut.BoostedJet.Number.Min", 999);
    _BoostedJetNumberMax = config -> GetValue("Cut.BoostedJet.Number.Max", 999);
  }
  //
  // also add these two cuts to the cut flow table
  ostringstream titleStr;
  titleStr <<"BoostedJet : " << "N >=" << _BoostedJetNumberMin;
  TString tTitle(titleStr.str().c_str());
  GetCutFlowTable()->AddCutToFlow("BoostedJet.Number.Min",tTitle);
  titleStr.str("");
  titleStr <<"BoostedJet : N <= "<<_BoostedJetNumberMax;
  tTitle=titleStr.str().c_str();
  GetCutFlowTable()->AddCutToFlow("BoostedJet.Number.Max",tTitle);
  // also add the overal pass rate
  titleStr.str("");
  titleStr <<"BoostedJet : "<< _BoostedJetNumberMin<<" <= N <= "<<_BoostedJetNumberMax;
  tTitle=titleStr.str().c_str();
  GetCutFlowTable()->AddCutToFlow("BoostedJet.Number.All",tTitle);

  
}//BookHistograms()

/******************************************************************************
 * Bool_t CutBoostedJetN::Apply()                                                   *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t CutBoostedJetN::Apply()
{
  //cout<<"Begin of CutBoostedJetN::Apply()"<<endl;

  // Get Event Tree
  EventContainer *evObj = GetEventContainer();

  Int_t BoostedJetNumber = evObj->boostedjets.size();
  //Int_t BoostedJetNumber = evObj->jets.size();
  // Fill the histogram before the cut
  _hBoostedJetNumberBefore -> Fill(BoostedJetNumber);

  // Flags 
  Bool_t BoostedJetNumberMinPass    = kTRUE;   // Event passes min cut
  Bool_t BoostedJetNumberMaxPass    = kTRUE;   // Event passes max cut

  // Cut on Min Pt of BoostedJet
  // Negative cut value for Min means there is no Min cut
  if( (_BoostedJetNumberMin != 999) && (BoostedJetNumber < _BoostedJetNumberMin) ){
     BoostedJetNumberMinPass = kFALSE;
     GetCutFlowTable()->FailCut("BoostedJet.Number.Min");
  }//if
  else{
    GetCutFlowTable() -> PassCut("BoostedJet.Number.Min");
  } //else

  // Cut on Min Pt of BoostedJet
  // Negative cut value for Min means there is no Min cut
  if( (_BoostedJetNumberMax != 999) && (BoostedJetNumber > _BoostedJetNumberMax) ){
     BoostedJetNumberMaxPass = kFALSE;
     GetCutFlowTable()->FailCut("BoostedJet.Number.Max");
  } //if
  else{
    GetCutFlowTable() -> PassCut("BoostedJet.Number.Max");
  } //else

  // Fill Histogram after cuts
  if(BoostedJetNumberMinPass && BoostedJetNumberMaxPass){
    _hBoostedJetNumberAfter -> Fill(BoostedJetNumber);
    GetCutFlowTable() -> PassCut("BoostedJet.Number.All");
  } //if
  else GetCutFlowTable() -> FailCut("BoostedJet.Number.All");

  return(BoostedJetNumberMinPass && BoostedJetNumberMaxPass);

} //Apply












