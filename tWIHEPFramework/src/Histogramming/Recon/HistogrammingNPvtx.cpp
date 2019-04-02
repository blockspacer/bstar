/******************************************************************************
 * HistogrammingNPvtx.cpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingNPvtx()              -- Parameterized Constructor            *
 *    ~HistogrammingNPvtx()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingNPvtx"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *    myTH1F* _hNPvtx                   -- Hist of MET                          *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      20 Nov 2006 - Modified by Bernard Pope                                *
 *      21 Mar 2007 - RS: Fill from event container, add sumET, mex, mey      *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingNPvtx.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingNPvtx::HistogrammingNPvtx(EventContainer *obj)                     *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Particle class                                                     *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingNPvtx::HistogrammingNPvtx(EventContainer *obj)
{
  SetEventContainer(obj);
} //HistogrammingNPvtx()

/******************************************************************************
 * HistogrammingNPvtx::~HistogrammingNPvtx()                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingNPvtx::~HistogrammingNPvtx()
{
  
} //HistogrammingNPvtx

/******************************************************************************
 * void HistogrammingNPvtx::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingNPvtx::BookHistogram(){
  
  // Histogram of MET
  _hNPvtx = DeclareTH1F("nPvtx","Number of primary vertex",50,0.,50.);
  _hNPvtx->SetXAxisTitle("nPvtx");
  _hNPvtx->SetYAxisTitle("Events");

  // Histogram of MEX 
  _hNPvtx_nJet = DeclareTH2F("nPvtx_nJet","Number of primary vertex vs number of Jets",50,0.,50., 50,0,50);
  _hNPvtx_nJet->SetXAxisTitle("nPvtx");
  _hNPvtx_nJet->SetYAxisTitle("nJet");

  // Histogram of MET
  _hNPvtx_nElectron = DeclareTH2F("nPvtx_nElectron","Number of primary vertex vs number of Electrons",50,0.,50.,50,0,50);
  _hNPvtx_nElectron->SetXAxisTitle("vPvtx");
  _hNPvtx_nElectron->SetYAxisTitle("nElectron");

  // Histogram of SUM ET
  _hNPvtx_nMuon = DeclareTH2F("nPvtx_nMuon","Number of primary vertex vs number of Muon",50,0.,50.,50,0,50);
  _hNPvtx_nMuon->SetXAxisTitle("vPvtx");
  _hNPvtx_nMuon->SetYAxisTitle("nMuon");

  // Histogram of MET
  _hNTrueInteractions = DeclareTH1F("nTrueInteractions","Number of true vertices",50,0.,50.);
  _hNTrueInteractions->SetXAxisTitle("nTrueInteractions");
  _hNTrueInteractions->SetYAxisTitle("Events");

  // Histogram of MEX 
  _hNTrueInteractions_nJet = DeclareTH2F("nTrueInteractions_nJet","Number of true vertices vs number of Jets",50,0.,50., 50,0,50);
  _hNTrueInteractions_nJet->SetXAxisTitle("nTrueInteractions");
  _hNTrueInteractions_nJet->SetYAxisTitle("nJet");

  // Histogram of MET
  _hNTrueInteractions_nElectron = DeclareTH2F("nTrueInteractions_nElectron","Number of true vertices vs number of Electrons",50,0.,50.,50,0,50);
  _hNTrueInteractions_nElectron->SetXAxisTitle("vTrueInteractions");
  _hNTrueInteractions_nElectron->SetYAxisTitle("nElectron");

  // Histogram of SUM ET
  _hNTrueInteractions_nMuon = DeclareTH2F("nTrueInteractions_nMuon","Number of true vertices vs number of Muon",50,0.,50.,50,0,50);
  _hNTrueInteractions_nMuon->SetXAxisTitle("vTrueInteractions");
  _hNTrueInteractions_nMuon->SetYAxisTitle("nMuon");

  // Histogram of MET
  _hNPUVertices = DeclareTH1F("nPUVertices","Number of PU vertices",50,0.,50.);
  _hNPUVertices->SetXAxisTitle("nPUVertices");
  _hNPUVertices->SetYAxisTitle("Events");

  // Histogram of MEX 
  _hNPUVertices_nJet = DeclareTH2F("nPUVertices_nJet","Number of PU vertices vs number of Jets",50,0.,50., 50,0,50);
  _hNPUVertices_nJet->SetXAxisTitle("nPUVertices");
  _hNPUVertices_nJet->SetYAxisTitle("nJet");

  // Histogram of MET
  _hNPUVertices_nElectron = DeclareTH2F("nPUVertices_nElectron","Number of PU vertices vs number of Electrons",50,0.,50.,50,0,50);
  _hNPUVertices_nElectron->SetXAxisTitle("vPUVertices");
  _hNPUVertices_nElectron->SetYAxisTitle("nElectron");

  // Histogram of SUM ET
  _hNPUVertices_nMuon = DeclareTH2F("nPUVertices_nMuon","Number of PU vertices vs number of Muon",50,0.,50.,50,0,50);
  _hNPUVertices_nMuon->SetXAxisTitle("vPUVertices");
  _hNPUVertices_nMuon->SetYAxisTitle("nMuon");


  //cout<<"end of HistogrammingNPvtx::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingNPvtx::Apply()                                         *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingNPvtx::Apply()
{
  //cout<<"Begin of HistogrammingNPvtx::Apply()"<<endl;
  
  EventContainer *evc = GetEventContainer();
  
  // Fill the histogram before the cut
  //cout<<"missing et is "<<MissingEt<<endl;

  // Fill Histograms
  //  if((evc->tightElectrons.size() >0)||(evc->tightMuons.size() >0)){
  _hNPvtx -> Fill(evc->nPvtx);
  _hNPvtx_nJet -> Fill(evc->nPvtx,evc->jets.size()); // njets
  _hNPvtx_nElectron -> Fill(evc->nPvtx,evc->tightElectrons.size()); // tight electron only
  _hNPvtx_nMuon -> Fill(evc->nPvtx,evc->tightMuons.size()); // tight muon only

  _hNTrueInteractions -> Fill(evc->trueInteractions);
  _hNTrueInteractions_nJet -> Fill(evc->trueInteractions,evc->jets.size()); // njets
  _hNTrueInteractions_nElectron -> Fill(evc->trueInteractions,evc->tightElectrons.size()); // tight electron only
  _hNTrueInteractions_nMuon -> Fill(evc->trueInteractions,evc->tightMuons.size()); // tight muon only

  _hNPUVertices -> Fill(evc->npuVertices);
  _hNPUVertices_nJet -> Fill(evc->npuVertices,evc->jets.size()); // njets
  _hNPUVertices_nElectron -> Fill(evc->npuVertices,evc->tightElectrons.size()); // tight electron only
  _hNPUVertices_nMuon -> Fill(evc->npuVertices,evc->tightMuons.size()); // tight muon only
  

  //  }
  //cout<<"End of HistogrammingNPvtx::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












