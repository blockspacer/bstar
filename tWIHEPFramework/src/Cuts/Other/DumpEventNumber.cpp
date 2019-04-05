/******************************************************************************
 * DumpEventNumber.cpp                                                       *
 *                                                                            *
 * Cuts on the requirement of a good primary vertex                           *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of DumpEventNumber class                              *
 *    DumpEventNumber()                     -- Parameterized Constructor         *
 *    ~DumpEventNumber()                    -- Destructor                        *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "DumpEventNumber"                *
 *                                                                            *
 * Private Data Members of DumpEventNumber class                                 *
 *    myTH1F* _hPrimaryVertexBefore;    -- Hist of PV before cut        *
 *    myTH1F* _hPrimaryVertexAfter;     -- Hist of PV after cut         *
 *                                                                            *
 * History                                                                    *
 *      5th July 2016 - Created by Duncan Leggat                              *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Cuts/Other/DumpEventNumber.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * DumpEventNumber::DumpEventNumber(EventContainer *EventContainerObj, TString electronType)    *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
DumpEventNumber::DumpEventNumber(EventContainer *EventContainerObj, TString fileNameToPrint)
{
  // Set Event Container
  SetEventContainer(EventContainerObj);

  _fileToPrint.open(fileNameToPrint);
  
} // DumpEventNumber


/******************************************************************************
 * DumpEventNumber::~DumpEventNumber()                                              *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
DumpEventNumber::~DumpEventNumber()
{
  _fileToPrint.close();
}//~DumpEventNumber

/******************************************************************************
 * void DumpEventNumber::BookHistogram()                                         *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void DumpEventNumber::BookHistogram(){
  


}//BookHistograms()

/******************************************************************************
 * Bool_t DumpEventNumber::Apply()                                               *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Cut on the selected trigger                                                *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t DumpEventNumber::Apply()
{

  //Get event container
  EventContainer *EventContainerObj = GetEventContainer();

  _fileToPrint << EventContainerObj->eventNumber << std::endl;
 
  return true;

} //Apply












