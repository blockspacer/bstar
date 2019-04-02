/******************************************************************************
 * HistogrammingEventList.cpp                                                 *
 *                                                                            *
 * Format write out selected Events: runNumber	EventNumber                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingEventList()              -- Parameterized Constructor      *
 *    ~HistogrammingEventList()             -- Destructor                     *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingEventList"   *
 *                                                                            *
 * Private Data Members of this class                                         *
 * - lots of histograms of muon properties                                    *
 *                                                                            *
 * History                                                                    *
 *      16 Mar 2011 - Created by Huaqiao ZHANG                                *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Other/HistogrammingEventList.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingEventList::HistogrammingEventList(EventContainer *obj)        *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingEventList::HistogrammingEventList(EventContainer *obj, string evtFile)
{
  SetEventContainer(obj);
  evtListFile.open(evtFile.c_str());
} //HistogrammingEventList()

/******************************************************************************
 * HistogrammingEventList::~HistogrammingEventList()                                    *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingEventList::~HistogrammingEventList()
{
  
} //HistogrammingEventList

/******************************************************************************
 * void HistogrammingEventList::BookHistogram()                                    *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingEventList::BookHistogram(){
  
  ////////////////////////////////////////////////////////////////
  // Histogram
  ////////////////////////////////////////////////////////////////
  _hNone = DeclareTH1F("Nothing","Nothing", 2, 0.0, 10.0);
  _hNone->SetXAxisTitle("N");
  _hNone->SetYAxisTitle("Events");

  //cout<<"end of HistogrammingEventList::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingEventList::Apply()                                     *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingEventList::Apply()
{
  //cout<<"Begin of HistogrammingEventList::Apply()"<<endl;

  // Get Event Tree
  EventContainer *evc = GetEventContainer();
  // Fill Histograms
  _hNone -> Fill(1);
  evtListFile << (int)(evc->runNumber) << "\t" << (int)(evc->eventNumber) << endl;
  //cout<<"End of HistogrammingEventList::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












