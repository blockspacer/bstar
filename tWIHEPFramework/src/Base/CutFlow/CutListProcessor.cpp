/******************************************************************************
 * CutListProcessor.cpp                                                       *
 *                                                                            *
 * Process the list of cuts and histogramming classes                         *
 *                                                                            *
 * Derived from BaseCut                                                       *
 * Parent to AnalysisMain                                                     *
 *                                                                            *
 * Public Member Functions of CutNJet class                                   *
 *    CutListProcessor()                -- Parameterized Constructor          *
 *    ~CutListProcessor()               -- Destructor                         *
 *    AddCut()                          -- Add a cut to the list to process   *
 *    BookHistogram()                   -- Book Histograms                    *
 *    Apply()                           -- Make cuts and fill histograms      *
 *    GetCutName()                      -- Get name of cut                    *
 *    PrintCut()                        -- Print name of cut                  *
 *    SetDebug()                        -- Set Debug on or off                *
 *    DeclareTH2F()                     -- Wrapper for TH2F                   *
 *    AddDefaultBaseCut                 -- Add template to every CutListProc  *
 *    template<class T>                 -- Template to be added to every list *
 *    GetCurrentPackage()               -- Returns current package            *
 *    GetCurrentEventContainer()        -- Returns event read                 *
 *                                                                            *
 * Private Data Members of CutNJet class                                      *
 *    string name                         -- Name of CutListProcessor         *
 *    vector<BaseCut *> _vectorCutAndHist -- Vector of cuts and histograms    *
 *    TH1F *_hCutAndHist                  -- Number of times each cut runs    *
 *    Bool_t _alwaysPass                  -- Always returns kFALSE            *
 *    CleanUpRun()                        -- Clean up after loop              *
 *    Bool_t _gDegub                      -- If debug is on                   *
 *                                                                            *
 *    template <class T> class AutoInsert_t                                   *
 *    template <class T> static AutoInsert_t<T> *AutoInsert_t_maker           *
 *    static std::vector<AutoInsert*> &defaultCuts(void)                      *
 *    static BaseCut *_currentBaseCut;                                        *
 *    static EventContainer *_currentEventContainer;                          *
 *                                                                            *
 * History                                                                    *
 *      10 Oct 2006 - RS: Ported D0/theory original to Atlas.                 *
 *      14 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      18 Jul 2007 - J. Holzbauer - added Cut Efficiency Histogram           *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Base/CutFlow/CutListProcessor.hpp"

#include "SingleTopRootAnalysis/Base/Histograms/myTH1F.hpp"
#include "TDirectory.h"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <functional>

using namespace std;


// Set Debug Flag
Bool_t CutListProcessor::_gDebug = kFALSE;

// Set _currentEventContainer and _currentBaseCut to NULL
EventContainer *CutListProcessor::_currentEventContainer = NULL;
BaseCut *CutListProcessor::_currentBaseCut = NULL;

// Set default cuts (applied to all CutListProcessors)
vector<CutListProcessor::AutoInsert*> &CutListProcessor::defaultCuts()
{
  static vector<CutListProcessor::AutoInsert*> gCuts;
  return gCuts;
} //defaultCuts


/**********************************************************************************
 * CutListProcessor::CutListProcessor (const string &name, Bool_t alwaysPass) *
 *                                                                                 *
 * Parameterized Constructor                                                       *
 *                                                                                 *
 * Input:  Name of class and a flag                                                *
 * Output: None                                                                    *
 ***********************************************************************************/
CutListProcessor::CutListProcessor (const string &name, const TString& condition,EventContainer *obj)
  : _name(name), _alwaysPass(false), _doOr(false), _cutFlowDetails("GlobalCutFlow")
{
  for_each (defaultCuts().begin(), defaultCuts().end(),
		 bind2nd (mem_fun (&AutoInsert::AddCuts), this));
  if(condition=="alwaysPass" || condition=="Pass" || condition=="PASS") _alwaysPass=true;
  if(condition=="doOr" || condition=="Or" || condition=="OR") _doOr=true;
    SetEventContainer(obj);
} //CutListProcessor()

/******************************************************************************
 * CutListProcessor::GetCutName()                                             *
 *                                                                            *
 * Return name of CutListProcessor                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: Returns name                                                       *
 ******************************************************************************/
string CutListProcessor::GetCutName()
{
  return _name;
} //GetCutName

//
// BookHistogram
//
//  Pass the booking on down! We don't have any of them ourselves for now!
//

/******************************************************************************
 * void CutListBrocessor::BookHistogram()                                     *
 *                                                                            *
 * Book Histograms                                                            *
 * Book the histogram giving the number of times the cuts and filling ran     *
 * Also book histograms                                                       *
 * The booking of the other histograms is done in the derived classes         *
 *                                                                            *
 * Input:  None                                                               *
 * Output: Prints each cuts/histogram class processes                         *
 ******************************************************************************/
void CutListProcessor::BookHistogram ()
{
  //
  // Let em know what happened. :-)
  //
  //  Watch for the case that we have no sub-cuts running. In this case
  // don't generate the who ran histo.
  //
  
  // Directory name of cut/hist to be constructed from steram
  ostringstream directoryNameSteam;
  // Directory name of cut/hist as a string
  string directoryName;

  // Only book histograms if cuts/histogram classes exist in vector
  if (_vectorCutAndHist.size() > 0) {

    // Book histograms giving the number of time the cuts and filling classes run
    _hCutAndHist = new TH1F ("CutAndHist", "Number of times each cut and histogram class ran",
			     _vectorCutAndHist.size(), 0, _vectorCutAndHist.size());
    _hCutAndHist -> SetXTitle("Cut/Hist");
    _hCutAndHist -> SetYTitle("Events");

    //
    // Now, book all the guys below that....
    //
    
    // Initialize number of Cut/Hist classes
    Int_t numCutAndHist = 0;
    
    // Write to screen
    cout << "<CutListProcessor::BookHistogram>      --------------------------------  " << endl;
    cout << "<CutListProcessor::BookHistogram>     |   Cut & Histo Classes Called   | " << endl;
    cout << "<CutListProcessor::BookHistogram>      --------------------------------  " << endl;
    cout.setf(ios::left);

    // Loop over cuts/hists in Vector
    for(vector<BaseCut *>::iterator itask = _vectorCutAndHist.begin(); itask != _vectorCutAndHist.end(); ++itask){
      // Set axis label
      _hCutAndHist -> GetXaxis() -> SetBinLabel( numCutAndHist+1, (*itask) -> GetCutName().c_str() );

      // Save directory name as string
      directoryNameSteam.str("");
      directoryNameSteam << setw(2) << setfill('0') << numCutAndHist << setw(0) << "_" << (*itask)->GetCutName();
      directoryName = directoryNameSteam.str();
      
      // Write directory name to screen
      cout << "<CutListProcessor::BookHistogram>     |  " <<  setw(27) << directoryName << "   |" << endl;

      // Get directory of cut/hist and switch to that directory
      TDirectory *directory = GetDirectory() -> mkdir(directoryName.c_str()," ");
      directory -> cd();

      // Set directory of current cut/hist and book that histogram
      (*itask) -> SetDirectory(directory);
      (*itask) -> BookHistogram();

      // Increment counter
      numCutAndHist++;
    } //for

    cout.setf(ios::right);
    cout << "<CutListProcessor::BookHistogram>      -------------------------------- " << endl;

  } //if

 // if this is an OR processor, also add this to the cutflow table:
  if(_doOr) {
    GetCutFlowTable()->AddCutToFlow(_name+"_OR",_name+", CutListprocessor with OR");
  }

} //BookHistogram()

/******************************************************************************
 * void CutListBrocessor::SaveHistogram()                                     *
 *                                                                            *
 * Call SaveHistogram once for each histogram class                           *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutListProcessor::SaveHistogram ()
{
  // Loop over cuts/hists in Vector
  for(vector<BaseCut *>::iterator itask = _vectorCutAndHist.begin(); itask != _vectorCutAndHist.end(); ++itask){
    // Debug info
    if (_gDebug) {
      cout << "<CutListProcessor::SaveHistogram> " << "DEBUG " << (*itask) -> GetCutName() << endl;
    } //if
    
      // Clean up and return kFALSE if no cut/hist exists in vector
    (*itask) -> SaveHistogram();
    
  }//for
  
} //SaveHistogram()

/******************************************************************************
 * void CutListBrocessor::PrintCut()                                          *
 *                                                                            *
 * Print cut info                                                             *
 * Calls the derived classes for each cut/hist                                *
 *                                                                            *
 * Input:  String used to indent                                              *
 * Output: Prints info for cuts                                               *
 ******************************************************************************/
void CutListProcessor::PrintCut (const string &indent)
{
  // Call printing method at the end of the run. 
  // The PrintCut() method is usefully implemented in Selection class so far.
  //  cout << indent << "<CutListProcessor::PrintCut> Results from Cut List '" << _name << "':" << endl;
  for(vector<BaseCut *>::iterator itask=_vectorCutAndHist.begin(); itask!=_vectorCutAndHist.end(); ++itask){
    (*itask)->PrintCut(indent + "  ");
  } //for
} //PrintCut()


/******************************************************************************
 * void CutListBrocessor::Apply()                                             *
 *                                                                            *
 * Fills histogram for cut/hist information                                   *
 * Cleans up if necessary                                                     *
 *                                                                            *
 * Input:  EventContainer Object                                              *
 * Output: Flag kTRUE if successful, kFALSE if not                            *
 ******************************************************************************/
Bool_t CutListProcessor::Apply (EventContainer *EventContainerObj)
{

  if(EventContainerObj)   SetEventContainer(EventContainerObj);
  else EventContainerObj=GetEventContainer();

  // Copy current values to OLD
  BaseCut *oldBaseCut = _currentBaseCut;
  EventContainer *oldEventContainer = _currentEventContainer;


  try {
    
    // Set new to current
    _currentEventContainer = EventContainerObj;
    _cutFlowDetails.SetEventContainer(EventContainerObj);

    // Loop on all cuts/histograms added
    Int_t count = 0;
    // if this is an OR processor, we need to keep track of whether any cut passed:
    bool passIfOr=false;

    for(vector<BaseCut *>::iterator itask = _vectorCutAndHist.begin(); itask != _vectorCutAndHist.end(); ++itask){
      // Debug info
      if (_gDebug) {
	cout << "<CutListProcessor::Apply> " << "DEBUG: Apply " << (*itask) -> GetCutName() << endl;
      } //if
      // Change directory and get current task
      (*itask) -> GetDirectory() -> cd();
      _currentBaseCut = *itask;
      
      // Fill histogram
      _hCutAndHist -> Fill(count);
      //
      // if this is an OR processor:
      if(_doOr) {
	passIfOr|= (*itask) -> Apply();
      }
      else {
	// this is the usual AND of the conditions, 
	// and we quit if this cut fails.
	// Clean up and return kFALSE if no cut/hist exists in vector
	if( !(*itask) -> Apply() ) {
	  CleanUpRun (_vectorCutAndHist.begin(), itask+1);
	  return _alwaysPass;
	} //if
      } // end of else: not an OR
      // Increment counter
      count++;
      
    }//for
    

    CleanUpRun(_vectorCutAndHist.begin(), _vectorCutAndHist.end());
    if(_doOr) {
      if(passIfOr) GetCutFlowTable()->PassCut(_name+"_OR");
      else GetCutFlowTable()->FailCut(_name+"_OR");
      return passIfOr;
    }
    else {
      // this is a regular AND processor
      return kTRUE;
    }
    
  } //try
  
  // Catch exception
  catch (...) {
    // Reset current to old
    _currentEventContainer = oldEventContainer;
    _currentBaseCut        = oldBaseCut;
    throw;
  } //catch


} //Apply()

/******************************************************************************
 * void CutListBrocessor::CleanUpRun()                                        *
 *                                                                            *
 * Calls DoneProcessing() for each element in Vector                          *
 *                                                                            *
 * Input:  Start and End of Vector                                            *
 * Output: None                                                               *
 ******************************************************************************/
void CutListProcessor::CleanUpRun(vector<BaseCut *>::iterator start,
				     vector<BaseCut *>::iterator end)
{
  while (start != end) {
    end--;
    (*end) -> DoneProcessingLoop();
  } //while
} //CleanUpRun

/******************************************************************************
 * void CutListBrocessor::AddCut()                                            *
 *                                                                            *
 * Add a cut/hist to Vector                                                   *
 *                                                                            *
 * Input:  Cut/hist class                                                     *
 * Output: None                                                               *
 ******************************************************************************/
void CutListProcessor::AddCut (BaseCut *CutAndHistBaseCut)
{
  _vectorCutAndHist.push_back(CutAndHistBaseCut);
  // also tell this cutandhist what cut flow table to use
  CutAndHistBaseCut->SetCutFlowTable(&_cutFlowDetails);
} //AddCut()

/******************************************************************************
 * void CutListBrocessor::PrintCutDetails()                                   *
 *                                                                            *
 * Fill Cut Efficiency Histogram, print cut details                           *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void CutListProcessor::PrintCutDetails()  { 
  gFile->cd(); 
  _cutFlowDetails.PrintCutFlow(); 
  _cutLines = _cutFlowDetails.GetCutLines();
  NumCutLines =  _cutLines.size();
  int histoparam = NumCutLines + 10;
   
  if(NumCutLines != 0){
    _hCutEff = new TH1F ("CutEff", "CutEfficiency, if Cut Efficiency is Not 1.0",
			 histoparam, 0, (double)histoparam);
    _hCutEff -> SetXTitle("Cuts");
    _hCutEff -> SetYTitle("yieldpassed/yieldexposed");

    _hCutEffNoEffect = new TH1F ("CutEffNoEffect", "CutEfficiency, if Cut Efficiency is 1.0",
			 histoparam, 0, (double)histoparam);
    _hCutEffNoEffect -> SetXTitle("Cuts");
    _hCutEffNoEffect -> SetYTitle("yieldpassed/yieldexposed");
    
    double startExposedYield = 0;
    double passedYield;
    for(int iCut=0; iCut<NumCutLines; iCut++) {
      TString name = _cutLines[iCut].second.GetTitle();
      double exposedYield = _cutLines[iCut].second.GetExposedYield();
      if(iCut == 0){
	startExposedYield = exposedYield;
      }
      passedYield = _cutLines[iCut].second.GetPassYield();
      double eff = passedYield/exposedYield;
      // cout<<name<<"  "<<exposedYield<<"  "<<passedYield<<"  eff= "<<eff<<endl;
      if(eff != 1){
	_hCutEff->Fill(name, eff);
      }else{
	_hCutEffNoEffect->Fill(name, eff);
      }
    }
    double totaleff = passedYield/startExposedYield;
    cout<<endl<<"<CutListProcessor::PrintCutDetails> Total Efficiency: "<<passedYield<<"/"<<startExposedYield<<" = "<<totaleff<<endl << endl;
    _hCutEff->Fill("all cuts", totaleff);
    
    _hCutEff->LabelsOption("u", "X");
    _hCutEff->LabelsDeflate("X");
    _hCutEff->SetMaximum(1.5);

    _hCutEffNoEffect->LabelsOption("u", "X");
    _hCutEffNoEffect->LabelsDeflate("X");
    _hCutEffNoEffect->SetMaximum(1.5);
  }
}
