/******************************************************************************
 * HistoCut.cpp                                                               *
 *                                                                            *
 * Abstract class which defines wrappers for TH1F, TH2F, and TGraph.          *
 * Wrappers are of type myTH1F, myTH2F, and TGraph.  The booking and filling  *
 * of the histograms is done in the derived classes.                          *
 *                                                                            *
 * Derived from BaseCut                                                       *
 * Parent to classes such as CutJetN, CutElectron, etc                       *
 *                                                                            *
 * Public Member Functions of CutNJet class                                   *
 *    HistoCut()                        -- Constructor                        *
 *    ~HistoCut()                       -- Destructor                         *
 *    DeclareTH1F()                     -- Wrapper for TH1F                   *
 *    DeclareTH2F()                     -- Wrapper for TH2F                   *
 *    DeclareTGraph()                   -- Wrapper for TGraph                 *
 *                                                                            *
 * Private Data Members of CutNJet class                                      *
 *    TObjArray _histograms              -- Collection of histograms in array *
 *    TObjArray _graphs                  -- Collection of graphs in array     *
 *                                                                            *
 * History                                                                    *
 *      10 Oct 2006 - RS: Ported D0/theory original to Atlas.                 *
 *      14 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"

using namespace std;

// Integrate classes into the Root system
ClassImp (HistoCut);

/******************************************************************************
 * HistoCut::HistoCut(EventContainer *EventContainerObj)                      *
 *                                                                            *
 * Constructor                                                                *
 *                                                                            *
 * Input:  Event Container object or nothing                                  *
 * Output: None                                                               *
 ******************************************************************************/
HistoCut::HistoCut(EventContainer *EventContainerObj)
{
  SetEventContainer(EventContainerObj);
} // HistoCut
  
/******************************************************************************
 * virtual HistoCut::~HistoCut()                                              *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Deletes all histograms and graphs in their respective TObjArrays           *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistoCut::~HistoCut(){
  for (int i=0;i<_histograms.GetEntries();i++) delete _histograms[i];
  for (int i=0;i<_graphs.GetEntries();i++) delete _graphs[i];
}  //HistoCut

/******************************************************************************
 * virtual void HistoCut::BookHistogram()                                     *
 *                                                                            *
 * Book Histograms - Do nothing for this class                                *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistoCut::BookHistogram()
{ 
    
} //BookHistogram
  
    
/******************************************************************************
 * virtual Bool_t CutJetN::Apply()                                           *
 *                                                                            *
 * Apply cuts and fill histograms - Do nothing for this class                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistoCut::Apply(){
    
  return kTRUE;

} //Apply
  

