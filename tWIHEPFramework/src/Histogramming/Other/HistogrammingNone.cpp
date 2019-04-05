/******************************************************************************
 * HistogrammingNone.cpp                                                      *
 *                                                                            *
 * Used as a placeholder.  Does not histogram anything                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingNone()              -- Parameterized Constructor           *
 *    ~HistogrammingNone()             -- Destructor                          *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingNone"        *
 *                                                                            *
 * Private Data Members of this class                                         *
 * - lots of histograms of muon properties                                    *
 *                                                                            *
 * History                                                                    *
 *      30 Mar 2009 - Created by P. Ryan                                      *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Other/HistogrammingNone.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingNone::HistogrammingNone(EventContainer *obj)                  *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event container                                                    *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingNone::HistogrammingNone(EventContainer *obj)
{
  SetEventContainer(obj);
} //HistogrammingNone()

/******************************************************************************
 * HistogrammingNone::~HistogrammingNone()                                    *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingNone::~HistogrammingNone()
{
  
} //HistogrammingNone

/******************************************************************************
 * void HistogrammingNone::BookHistogram()                                    *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingNone::BookHistogram(){
  
  ////////////////////////////////////////////////////////////////
  // Histogram
  ////////////////////////////////////////////////////////////////
  _hNone = DeclareTH1F("Nothing","Nothing", 2, 0.0, 10.0);
  _hNone->SetXAxisTitle("N");
  _hNone->SetYAxisTitle("Events");

  //cout<<"end of HistogrammingNone::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingNone::Apply()                                          *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingNone::Apply()
{
  //cout<<"Begin of HistogrammingNone::Apply()"<<endl;

  // Fill Histograms
  _hNone -> Fill(1);

  //cout<<"End of HistogrammingNone::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












