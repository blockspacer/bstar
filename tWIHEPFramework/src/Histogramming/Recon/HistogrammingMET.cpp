/******************************************************************************
 * HistogrammingMET.cpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMET()              -- Parameterized Constructor            *
 *    ~HistogrammingMET()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMET"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *    myTH1F* _hMET                   -- Hist of MET                          *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      20 Nov 2006 - Modified by Bernard Pope                                *
 *      21 Mar 2007 - RS: Fill from event container, add sumET, mex, mey      *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingMET.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingMET::HistogrammingMET(EventContainer *obj)                     *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Particle class                                                     *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingMET::HistogrammingMET(EventContainer *obj)
{
  SetEventContainer(obj);
} //HistogrammingMET()

/******************************************************************************
 * HistogrammingMET::~HistogrammingMET()                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingMET::~HistogrammingMET()
{
  
} //HistogrammingMET

/******************************************************************************
 * void HistogrammingMET::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingMET::BookHistogram(){
  
  // Histogram of MET
  _hMET = DeclareTH1F("MET","Missing transverse energy",200,0.,200.);
  _hMET->SetXAxisTitle("MET [GeV]");
  _hMET->SetYAxisTitle("Events");

  // Histogram of MEX 
  _hMEX = DeclareTH1F("MEX","X component of missing transverse energy",400,-200.,200.);
  _hMEX->SetXAxisTitle("MEX [GeV]");
  _hMEX->SetYAxisTitle("Events");

  // Histogram of MEY
  _hMEY = DeclareTH1F("MEY","Y component of missing transverse energy",400,-200.,200.);
  _hMEY->SetXAxisTitle("MEY [GeV]");
  _hMEY->SetYAxisTitle("Events");

  // Histogram of MET phi
  _hMETPhi = DeclareTH1F("METPhi","Phi of missing transverse energy",50,-3.3,3.3);
  _hMETPhi->SetXAxisTitle("#Phi_{MET} [GeV]");
  _hMETPhi->SetYAxisTitle("Events");

  // Histogram of MET corrected for xy
  _hMET_xy = DeclareTH1F("MET_xy","Missing transverse energy xy corrected",200,0.,200.);
  _hMET_xy->SetXAxisTitle("MET [GeV]");
  _hMET_xy->SetYAxisTitle("Events");

  // Histogram of MEX  corrected for xy
  _hMEX_xy = DeclareTH1F("MEX_xy","X component of missing transverse energy xy corrected",400,-200.,200.);
  _hMEX_xy->SetXAxisTitle("MEX [GeV]");
  _hMEX_xy->SetYAxisTitle("Events");

  // Histogram of MEY corrected for xy
  _hMEY_xy = DeclareTH1F("MEY_xy","Y component of missing transverse energy xy corrected",400,-200.,200.);
  _hMEY_xy->SetXAxisTitle("MEY [GeV]");
  _hMEY_xy->SetYAxisTitle("Events");

  // Histogram of MET phi corrected for xy
  _hMETPhi_xy = DeclareTH1F("METPhi_xy","Phi of missing transverse energy xy corrected",50,-3.3,3.3);
  _hMETPhi_xy->SetXAxisTitle("#Phi_{MET} [GeV]");
  _hMETPhi_xy->SetYAxisTitle("Events");

  // Histogram of SUM ET
  _hSumEt = DeclareTH1F("SumEt","Summed transverse energy",200,0.,2000.);
  _hSumEt->SetXAxisTitle("SumEt [GeV]");
  _hSumEt->SetYAxisTitle("Events");


  //cout<<"end of HistogrammingMET::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingMET::Apply()                                         *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingMET::Apply()
{
  //cout<<"Begin of HistogrammingMET::Apply()"<<endl;
  
  EventContainer *evc = GetEventContainer();
  
  // Fill the histogram before the cut
  //cout<<"missing et is "<<MissingEt<<endl;

  // Fill Histograms
  //  if((evc->tightElectrons.size() >0)||(evc->tightMuons.size() >0)){
  _hMET    -> Fill(evc->missingEt);
  _hMEX    -> Fill(evc->missingEx);
  _hMEY    -> Fill(evc->missingEy);
  _hMETPhi -> Fill(evc->missingPhi);

  _hMET_xy    -> Fill(evc->missingEt_xy);
  _hMEX_xy    -> Fill(evc->missingEx_xy);
  _hMEY_xy    -> Fill(evc->missingEy_xy);
  _hMETPhi_xy -> Fill(evc->missingPhi_xy);

  _hSumEt  -> Fill(evc -> sumEt);
  //  }
  //cout<<"End of HistogrammingMET::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












