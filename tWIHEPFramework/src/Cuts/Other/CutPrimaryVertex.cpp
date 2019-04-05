/******************************************************************************
 * CutPrimaryVertex.cpp                                                       *
 *                                                                            *
 * Cuts on the requirement of a good primary vertex                           *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of CutPrimaryVertex class                              *
 *    CutPrimaryVertex()                     -- Parameterized Constructor         *
 *    ~CutPrimaryVertex()                    -- Destructor                        *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutPrimaryVertex"                *
 *                                                                            *
 * Private Data Members of CutPrimaryVertex class                                 *
 *    myTH1F* _hPrimaryVertexBefore;    -- Hist of PV before cut        *
 *    myTH1F* _hPrimaryVertexAfter;     -- Hist of PV after cut         *
 *                                                                            *
 * History                                                                    *
 *      5th July 2016 - Created by Duncan Leggat                              *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Cuts/Other/CutPrimaryVertex.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * CutPrimaryVertex::CutPrimaryVertex(EventContainer *EventContainerObj, TString electronType)    *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
CutPrimaryVertex::CutPrimaryVertex(EventContainer *EventContainerObj)
{
  // Set Event Container
  SetEventContainer(EventContainerObj);
} // CutPrimaryVertex


/******************************************************************************
 * CutPrimaryVertex::~CutPrimaryVertex()                                              *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
CutPrimaryVertex::~CutPrimaryVertex()
{
  
}//~CutPrimaryVertex

/******************************************************************************
 * void CutPrimaryVertex::BookHistogram()                                         *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutPrimaryVertex::BookHistogram(){
  
  // ***********************************************
  // Make Strings for histogram titles and labels
  // ***********************************************  

  // Histogram Before Cut
  std::ostringstream histNameBeforeStream;
  histNameBeforeStream << "PrimaryVertexBefore";
  TString histNameBefore = histNameBeforeStream.str().c_str();

  std::ostringstream histTitleBeforeStream;
  histTitleBeforeStream << "PrimaryVertexBefore";
  TString histTitleBefore = histTitleBeforeStream.str().c_str();

  // Histogram After Cut
  std::ostringstream histNameAfterStream;
  histNameAfterStream << "PrimaryVertexAfter";
  TString histNameAfter = histNameAfterStream.str().c_str();

  std::ostringstream histTitleAfterStream;
  histTitleAfterStream << "PrimaryVertexAfter";
  TString histTitleAfter = histTitleAfterStream.str().c_str();

  // ***********************************************
  // Book Histograms
  // ***********************************************  

  // Histogram before cut
  _hPrimaryVertexBefore =  DeclareTH1F(histNameBefore.Data(), histTitleBefore.Data(), 100, 0.0, 100.);
  _hPrimaryVertexBefore -> SetXAxisTitle("Primary Vertex");
  _hPrimaryVertexBefore -> SetYAxisTitle("Events");

  // Histogram after cut
  _hPrimaryVertexAfter =  DeclareTH1F(histNameAfter.Data(), histTitleAfter.Data(), 100, 0.0, 100.);
  _hPrimaryVertexAfter -> SetXAxisTitle("Primary Vertex");
  _hPrimaryVertexAfter -> SetYAxisTitle("Events");

  //Additional histogram containing the number of true interactions
  _hNTrueInteractions = DeclareTH1F("nTrueInteractions","Number of true interactions",100,0.0,100.);
  _hNTrueInteractions -> SetXAxisTitle("N True Interactions");
  _hNTrueInteractions -> SetYAxisTitle("Events");

  _hNTrueIntsUnweighted = DeclareTH1F("nTrueIntsUnweighted","Number of true interactions unweighted",100,0.0,100.);
  _hNTrueIntsUnweighted -> SetXAxisTitle("N True Interactions");
  _hNTrueIntsUnweighted -> SetYAxisTitle("Events");

  // ***********************************************
  // Add these cuts to the cut flow table
  // ***********************************************
  ostringstream cutFlowTitleStream;
  ostringstream cutFlowNameStream;
  TString cutFlowTitle;
  TString cutFlowName;

  // Min cut
  cutFlowTitleStream << "Primary Vertex Cut";
  cutFlowTitle = cutFlowTitleStream.str().c_str();

  cutFlowNameStream << "PV";
  cutFlowName = cutFlowNameStream.str().c_str();

  GetCutFlowTable()->AddCutToFlow(cutFlowName.Data(), cutFlowTitle.Data());


  // ***********************************************
  // Get the PV IDs from the config file. Perhaps these should be included as an object
  // ***********************************************

  //First, get the config file
  EventContainer *EventContainerObj = GetEventContainer();
  TEnv * config = EventContainerObj->GetConfig();

  //Then set the cuts here.
  _ndofCut = config -> GetValue("ObjectID.Vertex.Ndof",0.0);
  _zCut = config -> GetValue("ObjectID.Vertex.Z",999.);
  _dxyCut = config -> GetValue("ObjectID.Vertex.Dxy",999.);
  


}//BookHistograms()

/******************************************************************************
 * Bool_t CutPrimaryVertex::Apply()                                               *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Cut on the selected trigger                                                *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t CutPrimaryVertex::Apply()
{

  EventTree *EventContainerObj = GetEventContainer()->GetEventTree();

  Bool_t hasPV = kFALSE;

  Bool_t passesMETFilterFlags = kTRUE;

  //  if ( !EventContainerObj->Flag_goodVertices || !EventContainerObj->Flag_HBHENoiseFilter || !EventContainerObj->Flag_HBHENoiseIsoFilter || !EventContainerObj->Flag_EcalDeadCellTriggerPrimitiveFilter || !EventContainerObj->Flag_globalTightHalo2016Filter) passesMETFilterFlags = kFALSE; //2016 version. 2016 halo filter not currently in full production, will be added next crab round.
  if ( !EventContainerObj->Flag_goodVertices || !EventContainerObj->Flag_HBHENoiseFilter || !EventContainerObj->Flag_HBHENoiseIsoFilter || !EventContainerObj->Flag_EcalDeadCellTriggerPrimitiveFilter || !EventContainerObj->Flag_CSCTightHaloFilter) passesMETFilterFlags = kFALSE; //old halo filter. Will be changed after re-crab.

  //  if (!EventContainerObj->Flag_METFilters) passesMETFilterFlags = kFALSE;

  tmpCounter = 0;

  for (unsigned int  i = 0; i < EventContainerObj->pvertex_ndof->size() ; i++){
    //    cout << EventContainerObj->pvertex_ndof->at(i) << " " << TMath::Abs(EventContainerObj->pvertex_z->at(i)) << " " << EventContainerObj->pvertex_dxy->at(i);
    if (EventContainerObj->pvertex_ndof->at(i) < _ndofCut) continue;
    if (TMath::Abs(EventContainerObj->pvertex_z->at(i)) > _zCut) continue;
    if (EventContainerObj->pvertex_dxy->at(i) > _dxyCut) continue;
    tmpCounter ++;
    hasPV = kTRUE;
  }
  //cout << endl;
  GetEventContainer()->SetnPvt(tmpCounter);

  // Fill the histograms before the cuts
  _hPrimaryVertexBefore->Fill(tmpCounter);

  // ***********************************************
  // Fill cut flow table                            
  // ***********************************************
                                                    
  // Names for Cut Flow Table                       

  ostringstream cutFlowNameStream;
                                
  TString cutFlowName;            

  cutFlowNameStream << "PV";
  cutFlowName = cutFlowNameStream.str().c_str();

  _hNTrueInteractions->Fill(EventContainerObj->trueInteractions);
  _hNTrueIntsUnweighted->FillWithoutWeight(EventContainerObj->trueInteractions);

  if (hasPV && passesMETFilterFlags) {
    _hPrimaryVertexAfter->Fill(tmpCounter);
    GetCutFlowTable()->PassCut(cutFlowName);
    return kTRUE;
  }
  else {
    GetCutFlowTable()->FailCut(cutFlowName);
    return kFALSE;
  }

} //Apply












