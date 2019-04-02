/******************************************************************************
 * AnalysisMain.hpp                                                           *
 *                                                                            *
 * This class controls the steering of the analysis.                          *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    AnalysisMain()                    -- Default Constructo                 *
 *    ~AnalysisMain()                   -- Destructor                         *
 *    OpenHistogramFile()               -- Open file for output histograms    *
 *    WriteHistogramFile()              -- Write output histograms            *
 *    CloseHistogramFile()              -- Close output histogram file        *
 *    GetHistogramFileName()            -- Get file name                      *
 *    Loop()                            -- Loop over events                   *
 *    BookHistogram()                   -- Not needed for this class          *
 *    Apply()                           -- Not needed for this class          *
 *                                                                            *
 * Private Data Members of AnalysisMain class                                 *
 *     TFile *_histogramFile            -- Output histogram file              *
 *     string _histogramFileName        -- Output histogram file name         *
 *     Int_t _eventLimit               -- Number of events to loop over       *
 *     Bool_t _doSkim                   -- Flag for skimming                  *
 *     TFile *_skimFile                 -- Output skim file                   *
 *     string _skimFileName             -- Name of the skim file              *
 *     TTree* _skimCollectionTree       -- Collection tree in the skim file   *
 *     TTree* _skimEventTree            -- Event tree in the skim file        *
 *     TTree* _skimTruthTree            -- Truth tree in the skim file        *
 *     TTree* _skimTriggerTree          -- Trigger tree in the skim file      *
 *                                                                            *
 * History                                                                    *
 *     17 June 2015  -  Huaqiao ZHANG Port from STAR                          *
 *****************************************************************************/

#ifndef AnalysisMain_h
#define AnalysisMain_h

#include <string>
#include <sstream>

#include "TFile.h"

#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include "SingleTopRootAnalysis/Base/CutFlow/CutListProcessor.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/AdditionalVarsProcessor.hpp"

class AnalysisMain: public EventContainer, public CutListProcessor, public AdditionalVarsProcessor
{
 public:
  // Default Constructor
  AnalysisMain();
  // Destructor
  ~AnalysisMain();

  // Open the output histogram file
  inline void OpenHistogramFile(const std::string& name) {
    _histogramFile = new TFile(name.c_str(), "RECREATE"); 
    _histogramFileName = name;

    // Check to see that file opened ok
    if( _histogramFile -> IsZombie() ){
      std::cout << "<AnalysisMain::OpenHistogramFile> " << "WARNING: ";
      std::cout << "Could not open " << _histogramFileName << std::endl;
    } //if
  } //OpenHistogramFile

  // Write histograms to open output file
  inline void WriteHistogramFile() { 
    if( _histogramFile ) _histogramFile -> Write(); 
  }

  // Close the output histogram file
  inline void CloseHistogramFile() { 
    //    if( _histogramFile -> IsOpen() ) _histogramFile -> Close(); 
    if(_histogramFile) {
      _histogramFile -> Close(); 
      delete _histogramFile; 
      _histogramFile = NULL; 
    } //if
  }//CloseHistogramFile

  // Finish up looping, write and close files
  inline void Finish() {
    SaveHistogram();
    WriteHistogramFile();
    CloseHistogramFile();

 if(chainConfig)    delete chainConfig;
  chainConfig       = NULL;
  
  }

  // Get the name of the histogram file
  inline std::string GetHistogramFileName() { return _histogramFileName; }
  

  // Same thing for the skim file
  // Open the output skim file
  inline void OpenSkimFile(const std::string& name) {
    _skimFile = new TFile(name.c_str(), "RECREATE"); 
    _skimFileName = name; 
  }

  // Write histograms to open output skim file
  inline void WriteSkimFile() { 
    if(_skimFile) { 
      //TDirectory *AODBookkeepingDir = _skimFile->mkdir("AODBookkeeping","AODBookkeepingD");
      //if (_hevents!=NULL){
      //  _hevents->Scale(1.0/_SkimMax);
      //  AODBookkeepingDir->Add(_hevents);
      //}
      //TDirectory *LooseTopInputsMetaDir = _skimFile->mkdir("LooseTopInputsTreeMeta","LooseTopInputsTreeMetaD");
      //if (_skimLooseTopInputsTreeMeta!=NULL){
      //  cout << _skimLooseTopInputsTreeMeta->GetName()<<" "<<GetTitle()<<endl;
      //  LooseTopInputsMetaDir->Add(_skimLooseTopInputsTreeMeta);
      //}
      //if(NULL != _skimDecisionTree){
      //	_skimFile = _skimEventTree -> GetCurrentFile();
      //}else if (NULL != _skimDecisionTree){
      //	_skimFile = _skimDecisionTree -> GetCurrentFile();
      //}else{
      //	cout<<"ERROR: Cannot write skim file.  skimEventTree and skimDecisionTree does not exist!!"<<endl;
      //}
      _skimFile -> Write();  
    } //if
    else{
      std::cout<<"<AnalysisMainhpp>: PROBLEM with skim file writing"<<std::endl;
    }
  }

  // Close the output skim histogram file
  inline void CloseSkimFile() { 
    if (_skimFile) {
      // _skimFile=_skimEventTree->GetCurrentFile(); 
      _skimFile->Close(); 
     delete _skimFile; 
     _skimFile=NULL; 
    } //if
  }

  // Get the name of the skim file
  std::string GetSkimFileName() { return _skimFileName; }

  // Get total event counts with MCatNLO weights
  Double_t GetTotalMCatNLOEvents() { return _totalMCatNLOEvents; }

  // Parse the command line to extract the names of the root files 
  // And put them into chains
  //Int_t ParseCmdLine(Int_t argc, char **argv, TChain *chainCollection, TChain *chainEVO, TChain *chainTruth);
  Int_t ParseCmdLine(Int_t argc, char **argv, TChain *chainEVO, TChain *chainTruth, TChain *chainTrigger, 
		     TChain *chainInfo,TChain *chainDecision, TChain *chainFastSim, TChain *chainBgkd);

  // Loop over events
  void Loop();

  // Declare these two because they are virtual in the class we are inheriting from
  // They are not needed for AnalysisMain class
  void BookHistogram() {};
  bool Apply() { return true; }
 
 private:  
  TFile *_histogramFile;            // Output histogram file
  std::string _histogramFileName;   // Output histogram file name
  Int_t _eventLimit;                 // Number of events to loop over from cmd line (set to 0 to loop over all events)
  Float_t  _EventNBeforePreselb;
  TBranch *_newBranchHFORb;
  TBranch *_newBranchNBPb;

  TChain *chainConfig;

  Int_t _HFORb;
  Bool_t _doSkim;                   // Should we write out selected events to a new skim?        
  TFile *_skimFile;                 // Output skim file
  std::string _skimFileName;        // Name of the skim file
  TTree* _skimEventTree;            // Event tree in the skim file
  TTree* _skimTruthTree;            // Truth tree in the skim file
  TTree* _skimTriggerTree;          // Trigger tree in the skim file
  TTree* _skimInfoTree;             // Info tree in the skim file
  TTree* _skimFastSimTree;          // FastSim tree in the skim file
  TTree* _skimDecisionTree;             // Info tree in the skim file
  TTree* _skimBkgdTree;             // bkgd tree in the skim file (for use with any trees hard coded into Analysis Main if statement)
  TTree* _skimConfigTreeMeta;  // ConfigMeta tree in skim file
  TTree* metaTree; 

  std::vector<int> _eventsToRunOn; //A vector of interesting event numbers. Mostly used for synch.

  Double_t _totalEvents;
  Double_t _totalMCatNLOEvents;
  TH1 *_hevents;
  TObjString *_skimLooseTopInputsTreeMeta;
  int _SkimMax;
 
   //Trig::D3PDBunchCrossingToolSA* m_tool;

  // Copy Constructor 
  // Will generate a compile-time error if called
  AnalysisMain(const AnalysisMain& oldAnalysisMain);

  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(AnalysisMain,0)
  ////////////////////// Do Not Put Anything Below this Line /////////////////////
};


#endif
