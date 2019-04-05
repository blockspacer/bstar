/******************************************************************************
 * AnalysisMain.cpp                                                           *
 *                                                                            *
 * This class controls the steering of the analysis.                          *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    AnalysisMain()                    -- Default Constructor                *
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
 *     TTree* _skimTruthAll             -- TruthAll tree in the skim file     *
 *     TTree* _skimTriggerTree          -- Trigger tree in the skim file      *
 *                                                                            *
 *                                                                            *
 * History                                                                    *
 *     4 Jul 2005 - Created by R. Schwienhorst at D0                          *
 *     7 Nov 2006 - Modified by R. Schwienhorst for ATLAS                     *
 *     9 Nov 2006 - Modified by P. Ryan - cleanup and reorganized             *
 *    21 Dec 2006 - Modified by R. Schwienhorst: Take out the Chain friends   *
 *                  so that it can run over more than a single file.          *
 *    19 Jan 2007 - Modified by R. Schwienhorst: Add trigger tree.            *
 *    4  May 2007 - P. Ryan - added TruthAll tree                             *
 *    27 Jun 2007 - P. Ryan - only fill Truth tree if it exists in root file  *
 *    28 Jul 2009 - J. Holzbauer - accept config file or command line inputs  *
 *   28  Oct 2009 - P. Ryan - update for v15051 - trees have changed          *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Base/Dictionary/AnalysisMain.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

#include <TStopwatch.h>
#include <TTreeIndex.h>
#include <TSystem.h>
#include <TChain.h>
#include <TGraphAsymmErrors.h>
#include "SingleTopRootAnalysis/Base/Histograms/utils.hpp"

using namespace std;

// Integrate classes into the Root system
ClassImp(AnalysisMain)
 

/******************************************************************************
 * AnalysisMain::AnalysisMain()                                               *
 *                                                                            *
 * Default Constructor                                                        *
 * Initialize private data members of AnalysisMain                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None
 ******************************************************************************/
AnalysisMain::AnalysisMain(): CutListProcessor("cutter"),            // New Instance of CutListProcessor class
                              AdditionalVarsProcessor(),
			      _histogramFile(NULL),                  // Output file points to NULL
			      _eventLimit(0),                        // Initialize max number of events processed
			      _histogramFileName("histograms.root"), // Initialize output file name_skimFileName("NONE"),                 // Initialize skim output file name
			      _doSkim(kFALSE),                       // Don't skim by default
			      //_skimCollectionTree(NULL),           // skimTree points to NULL
			      _skimEventTree(NULL),                  // skimTree points to NULL
			      _skimTruthTree(NULL),                  // skimTree points to NULL
			      _skimTriggerTree(NULL),                // skimTree points to NULL
			      _skimInfoTree(NULL),                   // skimTree points to NULL
			      _skimDecisionTree(NULL),                   // skimTree points to NULL
			      _skimFastSimTree(NULL),                 // skimTree points to NULL
			      _skimBkgdTree(NULL),                 // skimTree points to NULL
			      _hevents(NULL),                         // events histo points to NULL
			      _totalEvents(0),                        //events double is 0
                              _totalMCatNLOEvents(0),                  //events double is 0
			      _EventNBeforePreselb(-999),
			      _newBranchHFORb(NULL),
                              _newBranchNBPb(NULL),
			      _HFORb(-999),
			      
			      chainConfig(NULL),
			      
			      
			      _SkimMax(1)
{

  cout << "AnalysisMain::AnalysisMain " << "Constructor called" << endl;
} // AnalysisMain::AnalysisMain()


/******************************************************************************
 * AnalysisMain::~AnalysisMain()                                              *
 *                                                                            *
 * Destructor                                                                 *
 * Do nothing except destruct                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
AnalysisMain::~AnalysisMain()
{
  //if(m_tool) delete m_tool;

} // AnalysisMain::~AnalysisMain()


/******************************************************************************
 * AnalysisMain::ParseCmdLine()                                               *
 *                                                                            *
 * Parse the Command Line                                                     *
 *                                                                            *
 * Input:  Command Line Arguements and Chains of Events                       *
 * Output: Return 0 if successful 1 if unsuccessful (and print error)         *        
 ******************************************************************************/
Int_t AnalysisMain::ParseCmdLine(int argc, char **argv, TChain *chainEV0, TChain *chainTruth, TChain *chainTrigger, 
				 TChain *chainInfo, TChain *chainDecision, TChain *chainFastSim0, TChain *chainBkgd)
{
  /////////////////////////////////////////////////////////////////////////////////
  // Declare and Initialize Variables
  /////////////////////////////////////////////////////////////////////////////////
  // Offset from begining of chain
  Int_t evtOffSet = 0;

  chainConfig    = new TChain("ConfigTreeMeta/BunchConfTree"); //special tree for bunch crossing program

  // Variables for determining whether a list or single file is read in
  Int_t  fileListIndex   = 0;      // Parameter number of file-list name in cmd line arguments
  Int_t  fileSingleIndex = 0;      // Parameter number of file name in cmd line arguements
  Bool_t fileListFlag    = kFALSE; // Flag signifying use of file-list
  Bool_t fileSingleFlag  = kFALSE; // Flag signifying use of a single root file
  Bool_t fileConfigSpecified = kFALSE; //Flag signifying if config file read in yet and if commandlines can override it.
  SetSourceName("NONE");
  SetSourceNumber(0);

  if(NULL == chainEV0) {
    // no event tree? Check if this is supposed to be fast sim
    if(NULL == chainFastSim0) {
      cout << "<AnalysisMain::ParseCmdLine> " << "No valid input Event or FastSim Chain." << endl;
      return -1;
    } else {
      SetDoFastSim(true);
    }
  }//if

  // Truth Tree
  if(NULL == chainTruth) {
    cout << "<AnalysisMain::ParseCmdLine> " << "No valid input Truth Chain." << endl;
    //cout << "<AnalysisMain::ParseCmdLine> " << "This is ok for now..." << endl;
    return -1;
  }//if

  // Trigger Tree
  if(NULL == chainTrigger) {
    cout << "<AnalysisMain::ParseCmdLine> " << "No valid input Trigger Chain." << endl;
    //cout << "<AnalysisMain::ParseCmdLine> " << "This is ok for now..." << endl;
    return -1;
  }//if

  // Info Tree
  if(NULL == chainInfo) {
    cout << "<AnalysisMain::ParseCmdLine> " << "No valid input Info Chain." << endl;
    //cout << "<AnalysisMain::ParseCmdLine> " << "This is ok for now..." << endl;
    return -1;
  }//if

  // Decision Tree
  if(NULL == chainDecision) {
    cout << "<AnalysisMain::ParseCmdLine> " << "No valid input Decision Chain." << endl;
    //cout << "<AnalysisMain::ParseCmdLine> " << "This is ok for now..." << endl;
    return -1;
  }//if

// Bkgd Tree
  if(NULL == chainBkgd) {
    cout << "<AnalysisMain::ParseCmdLine> " << "No valid input Bkgd Chain." << endl;
    //cout << "<AnalysisMain::ParseCmdLine> " << "This is ok for now..." << endl;
    return -1;
  }//if
// Config Tree
  if(NULL == chainConfig) {
    cout << "<AnalysisMain::ParseCmdLine> " << "No valid input Config Chain." << endl;
    //cout << "<AnalysisMain::ParseCmdLine> " << "This is ok for now..." << endl;
    return -1;
  }//if

 
  /////////////////////////////////////////////////////////////////////////////////
  // Loop over arguements in command line 
  /////////////////////////////////////////////////////////////////////////////////
  cout << "<AnalysisMain::ParseCmdLine> " << "Parsing command line." << endl;

  // Print message if not enough arguements on command line
  if(argc < 3) {
    cout << "<AnalysisMain::ParseCmdLine> "<<"****************************************************************************" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"* Program: AnalysisMain                                                    *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                                                                          *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"* Purpose: analyze a parton level single top root file.                    *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                                                                          *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"* Usage: AnalysisMain -inlist filename                                     *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  You need to specify a root file or list of files                        *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  Note that commandline inputs will overide those in a config file        *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  BUT the config file must be specified FIRST                             *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  Input/Output files must be specified from the commandline               *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                                                                          *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -inlist file    - open a file with a list of root files.                *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -infile file    - load a single root file.                              *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                                                                          *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -config cfgfile - use this configuration file.                          *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                                                                          *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -source name    - name of the source we are processing.                 *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                    Can also be passed through the list file.             *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -sourcenum #    - number for the source we are processing.              *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                    Can also be passed through the list file.             *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -num #          - number of events to process.                          *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -offSet #       - start processing not at 0 but here.                   *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -skip #         - same as -offSet.                                      *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -hfile file     - filename for output histograms.                       *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -skimfile file  - filename to write selected events to.                 *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -JESError #     - JES uncertainty                                       *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                    Can also be passed through config file.               *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -JESShift #     - shift jet energy up or down by this amount            *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                    Can also be passed through config file.               *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -bTagAlgo name  - use this algorithm for -btagging (TRFIP2D, TRFSV1)    *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                    Can also be passed through config file.               *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -bTagCut #      - use this to specify cut on b-tagging weight (4.5, etc)*" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                    used for b-tagging uncertainty.               *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                    Can also be passed through config file.               *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -MisTagCut #    - use this to specify cut on b-tagging weight (4.5, etc)*" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                    used for mis-tag rate uncertainty.               *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                    Can also be passed through config file.               *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*  -evtList #       A comma separated list of events to run on             *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"*                                                                          *" << endl;
    cout << "<AnalysisMain::ParseCmdLine> "<<"****************************************************************************" << endl;
    return 1;
  } //if

  cout <<"in"<<endl;
  for (Int_t i = 1; i < argc; ++i) {
    //    cout << "Command line parameter " << i << " is " << argv[i] << endl;
    // FLAG: List of root files
    if( !strcmp(argv[i],"-inlist") ) {
      // Check if file name parameter is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing file-list name in cmd line" << endl;
	return 1;
      } //if argc < i+1
      fileListIndex = i+1;
      ++i;
      fileListFlag = kTRUE;
    } //if list of root files

    // FLAG: single root file
    else if( !strcmp(argv[i],"-infile") ) { 
      // Check if file name parameter is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing input-file name in cmd line" << endl;
	return 1;
      }//if 
      fileSingleIndex = i + 1;
      ++i;
      fileSingleFlag = kTRUE;
    }//if single root file

    else if( !strcmp(argv[i],"-evtList") ) { 
      // Check if file name parameter is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing evt-list-file name in cmd line" << endl;
	return 1;
      }//if 
      std::string listInStr = argv[i+1];
      ++i;
      std::stringstream ss(listInStr);
      
      int i;

      while (ss.good()) {
	std::string substr;
	getline(ss,substr,',');
	_eventsToRunOn.push_back(std::stoi(substr));
      }
      cout << "<AnalysisMain::ParseCmdLine> " << "Running over only events: ";
      for (auto const evtID : _eventsToRunOn) cout << evtID << ",";
      cout << endl;
    }//if single root file

    // FLAG: config file
    else if (!strcmp(argv[i],"-config")) {  
      // Check if file name parameter is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing configuration file name in cmd line" << endl;
	return 1;
      }//if
      TString confFileName = argv[i+1];
      // Check if this input file exists!
      FILE *confFile;
      if(!(confFile = fopen(confFileName.Data(), "r"))) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Input configuration file " << confFileName.Data() 
	     << " does not exist!" << endl;
	return -1;
      } //if
      fclose(confFile);
      _config.ReadFile(confFileName.Data(), kEnvChange);
      ++i;
      ////////////////////////////////////////////////////////////////////////////////
      // Parse the configuration file (part 1)- does involve options commandline inputs can override
      /////////////////////////////////////////////////////////////////////////////////
	
	// Debug level
	if(_config.Defined("DebugLevel")) {
	  SetDebugLevel(_config.GetValue("DebugLevel",0));
	  if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Parsing configuration file " << endl;
	  if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Setting debug level to " <<GetDebugLevel() << endl;
	} //if
	if(_config.Defined("Include.ObjectIDFile")) {
	  TString objectFile=_config.GetValue("Include.ObjectIDFile","NONE");
	  // check if this input file exists!
	  FILE *objectidFile;
	  if(!(objectidFile = fopen(objectFile.Data(), "r"))) {
	    cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Object ID file " << objectFile.Data()<<" does not exist!" << endl;
	    return -1;
	  } //if
	  fclose(objectidFile);
	  // The file exists, read it in:
	  if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Reading object ID config file: " << objectFile.Data() << endl;    
	  _config.ReadFile(objectFile.Data(), kEnvChange);

	  SetJESError(_config.GetValue("ObjectID.JES.Error",0.)); //percent uncertainty, so 0.1 is 10% JES uncertainty
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting jet energy scale error to "<<JESError()<<endl;
	  SetJESShift(_config.GetValue("ObjectID.JES.Shift",0)); // either +1 (for +1 sigma, -1 (for -1 sigma) or 0 (no shift, also default)
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting jet energy scale shift to "<<JESShift()<<endl;

	  SetJESUShift(_config.GetValue("ObjectID.JES.UniformShift",0)); // either +1 (for +1 sigma, -1 (for -1 sigma) or 0 (no shift, also default)
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting jet energy scale uniform shift to "<<JESUShift()<<endl;
	  SetJESPtShift(_config.GetValue("ObjectID.JES.PtShift",0)); // either +1 (for +1 sigma, -1 (for -1 sigma) or 0 (no shift, also default)
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting jet energy scale shift to "<<JESPtShift()<<endl;
	  SetJESLowPtShift(_config.GetValue("ObjectID.JES.PtLowShift",0)); // 999 for on, 0 for off
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting jet energy scale low Pt shift to "<<JESLowPtShift()<<endl;
	  SetJESHighPtShift(_config.GetValue("ObjectID.JES.PtHighShift",0)); // 999 for on, 0 for off
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting jet energy scale high Pt shift to "<<JESHighPtShift()<<endl;
	  SetJESEtaShift(_config.GetValue("ObjectID.JES.EtaShift",0)); // either +1 (for +1 sigma, -1 (for -1 sigma) or 0 (no shift, also default)
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting jet energy scale shift to "<<JESShift()<<endl;

	  SetBTagAlgo(_config.GetValue("ObjectID.BTAG.BTagAlgo","default")); //algorithm used for b-tagging, default is tree default (IP3D+SV1)
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting b-tagging algorithm to "<<BTagAlgo()<<endl;
	  SetBTagCut(_config.GetValue("ObjectID.BTAG.BTagCut",-999.)); //cut on weight, default of -999 is probably NOT a good value for most analyses
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting b-tagging weight cut to "<<BTagCut()<<endl;
	  SetMisTagCut(_config.GetValue("ObjectID.BTAG.MisTagCut",-999.)); //cut on weight, default of -999 is probably NOT a good value for most analyses
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting mis-tagging weight cut to "<<MisTagCut()<<endl;
	  SetcelloutShift(_config.GetValue("ObjectID.cellout.Shift",0.));
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting cellout shift to "<<celloutShift()<<endl;
	  SetsoftjetShift(_config.GetValue("ObjectID.softjet.Shift",0.));
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting softjet shift to "<<softjetShift()<<endl;
	  SetpileupShift(_config.GetValue("ObjectID.pileup.Shift",0.));
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting pileup shift to "<<pileupShift()<<endl;
	  
	  TString confFileName =_config.GetValue("ObjectID.JESupconfig","");
	  if (confFileName!=""){
	    // Check if this input file exists!
	    FILE *confFile;
	    if(!(confFile = fopen(confFileName.Data(), "r"))) {
	      cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: JESup configuration file " << confFileName.Data() 
		   << " does not exist!" << endl;
	      return -1;
	    } //if
	    fclose(confFile);
	    _JESconfig.ReadFile(confFileName.Data(), kEnvChange);
	    _JESconfigread=true;
	  }
	} //if
	fileConfigSpecified = kTRUE;
    }//if config file

    // FLAG: source name
    else if (!strcmp(argv[i],"-source")) { 
      // Check if file name parameter is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for source name -source!" << endl;
	return 1;
      }
      SetSourceName(argv[i+1]);
      ++i;
    } //if source name
  
    // FLAG: source number
    else if (!strcmp(argv[i],"-sourcenum")) { 
      // Check if file name parameter is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for source number -sourcenum!" << endl;
	return 1;
      }
      SetSourceNumber(atoi(argv[i+1]));
      ++i;
    } //if source number

    // FLAG: Max Number of Events
    else if (!strcmp(argv[i],"-num")) {  
      // Check if max number is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for number of events -num!" << endl;
	return 1;
      } //if
      // Convert command line string into integer
      _eventLimit = atoi(argv[i+1]) + evtOffSet;
      ++i;
    } //if num events

    // FLAG: Event Offset
    else if (!strcmp(argv[i],"-offSet") || !strcmp(argv[i],"-skip")) {  
      // Check if offset number is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for number of events -offSet or -skip!" << endl;
	return 1;
      } //if
      // Convert command line string into integer
      evtOffSet=atoi(argv[i+1]);
      if(_eventLimit>0) _eventLimit+=evtOffSet;
      ++i;
    } //if event offset

    // FLAG: Output histogram file
    else if (!strcmp(argv[i],"-hfile")) {  
      // Check if file name parameter is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for histogram filename -hfile!" << endl;
	return 1;
      } //if
      _histogramFileName=argv[i+1];
      ++i;
    } //if output histogram file name

    // FLAG: Skim file
    else if (!strcmp(argv[i],"-skimfile")) {  
      // Check if file name parameter is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for skim filename -skimfile!" << endl;
	return 1;
      } //if
      _doSkim = kTRUE;
      SetDoSkim(kTRUE);
      _skimFileName = argv[i+1];
      ++i;
    } //if skim file
    else if (!strcmp(argv[i],"-SkimMax")) {  
      // Check if file name parameter is in command line
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for skim filename -SkimMax!" << endl;
	return 1;
      } //if
      _SkimMax=atoi(argv[i+1]);
    }

    // Jet Energy Scale
    else if (!strcmp(argv[i],"-JESError") || !strcmp(argv[i],"-jesError")) {  
      // Check if jes shift is specified
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for number of jet energy scale -JESError or -jesError!" << endl;
	return 1;
      } //if
      SetJESError(atof(argv[i+1]));
      if (fileConfigSpecified == kTRUE){
	cout << "<AnalysisMain::ParseCmdLine> " << "Setting jet energy scale error to "<<JESError()<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: If this value was specified in the config file, it is overwritten"<<endl;
      }else{
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: Config file not specified yet "<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "jet energy scale error will be set to config file value if it exists and "<<JESError()<<" otherwise"<<endl;
      }
      ++i;
    }// if JESError
    else if (!strcmp(argv[i],"-JESShift") || !strcmp(argv[i],"-jesShift")) {  
      // Check if jes shift is specified
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for number of jet energy shift -JESShift or -jesShift!" << endl;
	return 1;
      } //if
      SetJESShift(atoi(argv[i+1])); 
      if (fileConfigSpecified == kTRUE){
	cout << "<AnalysisMain::ParseCmdLine> " << "Setting jet energy scale shift to "<<JESShift()<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: If this value was specified in the config file, it is overwritten"<<endl;
      }else{
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: Config file not specified yet "<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "jet energy scale shift will be set to config file value if it exists and "<<JESShift()<<" otherwise"<<endl;
      }
      ++i;
    }// if JES
    else if (!strcmp(argv[i],"-BTagAlgo") || !strcmp(argv[i],"-bTagAlgo")) {  
      // Check if b-tagging algorithm has been specified
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for b-tag algorithm -BTagAlgo or -bTagAlgo!" << endl;
	return 1;
      } //if
      SetBTagAlgo(argv[i+1]);  
      if (fileConfigSpecified == kTRUE){
	cout << "<AnalysisMain::ParseCmdLine> " << "Setting b-tagging algorithm to "<<BTagAlgo()<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: If this value was specified in the config file, it is overwritten"<<endl;
      }else{
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: Config file not specified yet "<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "b-tagging algorithm will be set to config file value if it exists and "<<BTagAlgo()<<" otherwise"<<endl;
      }
      ++i;
    }// if JES
    else if (!strcmp(argv[i],"-BTagCut") || !strcmp(argv[i],"-bTagCut")) {  
      // Check if b-tagging algorithm has been specified
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for b-tag weight cut -BTagCut or -bTagCut!" << endl;
	return 1;
      } //if
      SetBTagCut(atof(argv[i+1]));     
      if (fileConfigSpecified == kTRUE){
	cout << "<AnalysisMain::ParseCmdLine> " << "Setting b-tagging weight cut to "<<BTagCut()<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: If this value was specified in the config file, it is overwritten"<<endl;
      }else{
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: Config file not specified yet "<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "b-tagging weight cut will be set to config file value if it exists and "<<BTagCut()<<" otherwise"<<endl;
      }
      ++i;
    }// if JES
    else if (!strcmp(argv[i],"-MisTagCut") || !strcmp(argv[i],"-misTagCut")) {  
      // Check if mis-tagging algorithm has been specified
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for mis-tag weight cut -MisTagCut or -misTagCut!" << endl;
	return 1;
      } //if
      SetMisTagCut(atof(argv[i+1]));     
      if (fileConfigSpecified == kTRUE){
	cout << "<AnalysisMain::ParseCmdLine> " << "Setting mis-tagging weight cut to "<<MisTagCut()<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: If this value was specified in the config file, it is overwritten"<<endl;
      }else{
	cout << "<AnalysisMain::ParseCmdLine> " << "NOTE: Config file not specified yet "<<endl;
	cout << "<AnalysisMain::ParseCmdLine> " << "mis-tagging weight cut will be set to config file value if it exists and "<<MisTagCut()<<" otherwise"<<endl;
      }
      ++i;
    }// if JES
  }//for loop over command line arguments
  
  /////////////////////////////////////////////////////////////////////////////////
  // Chain together root files
  /////////////////////////////////////////////////////////////////////////////////
  Int_t iChain;                        // Number of files in each line of chain file (after wildcard expansion)
  Bool_t doPathWarning   = kTRUE;     // Prints warning message about directories
  Bool_t checkTrigger    = kTRUE;     // if false, stop checking for existence of Trigger tree
  Bool_t checkInfo       = kTRUE;     // if false, stop checking for existence of Info tree
  Bool_t checkDecision       = kTRUE;     // if false, stop checking for existence of Decision tree
  Bool_t checkBkgd       = kTRUE;     // if false, stop checking for existence of Decision tree
  Bool_t checkConfig       = kTRUE;     // if false, stop checking for existence of Decision tree
  SetDoTruth(kTRUE);                  // set truth to true initially

  //std::cout << "<AnalysisMain::ParseCmdLine> " << "DoTruth: " << DoTruth() << std::endl;
  

  if (fileListFlag) {
    // Use fstream
    ifstream chainFile(argv[fileListIndex]);
    Int_t chainLineCounter = 0;  // Total number of lines the file list
    Int_t chainFileCounter = 0;  // Total number of files in the file list (after wildcard expansion)
    
    if (chainFile) { // file with list exists
      string fileInChainName;
      
      cout << "<AnalysisMain::ParseCmdLine> " << "Chaining files listed in: ";
      cout << argv[fileListIndex] << endl;
      while (getline(chainFile, fileInChainName)) {
	// Check for comment lines or empty lines
	if(fileInChainName.find("#") != string::npos || fileInChainName.size()<1) {
	  continue;
	}
	// the line starting with Name: contains the name for this source.
	// just for informational purposes:
	else if(fileInChainName.find("Name") != string::npos) {
	  string::size_type pos = fileInChainName.rfind(":");
	  //string sname = fileInChainName.substr(pos+1);
	  string sname = fileInChainName.substr(pos+2);
	  cout<<"<AnalysisMain::ParseCmdLine> Name of this source: "<< fileInChainName<<endl;
	  SetSourceName(sname.c_str());
	  continue;
	}
	// If the tree name is something unusual, change the chain here
	else if(fileInChainName.find("InTree") != string::npos) {
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting source tree to name to: " << fileInChainName.c_str() << endl;
	  string::size_type pos = fileInChainName.rfind(":");
	  string snu = fileInChainName.substr(pos+2);
	  delete chainEV0;
	  TChain *chainEV0 = new TChain(TString(snu));
	  continue;
	}
	// The line starting with Number: contains the source number
	else if(fileInChainName.find("Number") != string::npos) {
	  cout << "<AnalysisMain::ParseCmdLine> " << "Setting source number to: "<< fileInChainName.c_str() << endl;
	  string::size_type pos = fileInChainName.rfind(":");
	  string snu = fileInChainName.substr(pos+1);
	  SetSourceNumber(atoi(snu.c_str()));
	  continue;
	} //if
	else {
	  cout << "<AnalysisMain::ParseCmdLine> " << "Reading root tree from " << fileInChainName.c_str() << endl;
	  
	  // Check if the trees exist in the root files
	  // If the event tree doesn't exist, then we do fast simulation

	  // Pointer to input file
	  TFile *inFilePtr = new TFile(fileInChainName.c_str());

	  // *****************************
	  // Check Recon tree
	  // *****************************
	  if(! (inFilePtr->Get(chainEV0->GetName()))) {
	    cout<<"<AnalysisMain::ParseCmdLine> **WARNING** No event tree found: "<<chainEV0->GetName()<<". Checking for FastSim tree..."<<endl;
	    // check if the event tree and trigger trees exists in this file
	    // if the event tree doesn't exist, then we do fast simulation
	    if(inFilePtr->Get(chainFastSim0->GetName())) {
	      cout<<""<<endl;
	      SetDoFastSim(true);
	    }
	    else {
	      cout<<"<AnalysisMain::ParseCmdLine> **ERROR** FastSim Tree not found either! Skipping this file."<<endl;
	      inFilePtr->Close();
	      delete inFilePtr;
	      continue;
	    } // end of else: no fastsim tree found.
	  } // end of if not found event tree 

	  // *****************************
	  // Check Truth tree
	  // *****************************
          // First check that command line parameter was passed
          if(DoTruth()){
            if(inFilePtr->Get(chainTruth->GetName())) {
              //cout<<"Found Truth object in file!"<<endl;
              SetDoTruth(kTRUE);
              cout << "<AnalysisMain::ParseCmdLine> **INFO** Truth required and tree found." <<endl;
            } //if
            else {
              cout << "<AnalysisMain::ParseCmdLine> **WARNING** Truth required but tree not found. Running without it!" <<endl;
              SetDoTruth(kFALSE);
            } //else
          } //if

	  // *****************************
	  if(checkConfig) {
	    if(inFilePtr -> Get(chainConfig -> GetName())) {
	      //cout<<"Found Config object in file!"<<endl;
	      SetDoConfig(kTRUE);
	    } //if
	    else {
	      cout << "<AnalysisMain::ParseCmdLine> **WARNING** Config tree not found. Running without it!" <<endl;
	      SetDoConfig(kFALSE);
	      checkConfig = kFALSE;
	    } //else
	  } // end of: if CheckConfig
 
	  // ****************************************
	  // Add files to chain for the various trees
	  // ****************************************
	  // FastSim / Recon
	  if(DoFastSim()) iChain = chainFastSim0 -> Add(fileInChainName.c_str());
	  else            iChain = chainEV0      -> Add(fileInChainName.c_str());

	  // Truth
	  if(DoTruth()) {
            //cout<<"Found Truth object in file!"<<endl;
            chainTruth -> Add(fileInChainName.c_str());
            SetDoTruth(kTRUE);
          } //if

	  if(DoConfig()) {
	    chainConfig -> Add(fileInChainName.c_str());
	    SetDoConfig(true);
	  } //if
	  //chainMVA      -> Add(fileInChainName.c_str());
	  /////////////////////////////////////////////////////////////////////////////////////////
	  //Get total number of events from AOD
	  //////////////////////////////////////////////////////////////////////////////////////  
	  if(DoTruth()||DoSkim()) {
	    //inFilePtr->GetObject("/LooseTopInputsTreeMeta/Lumi",_skimLooseTopInputsTreeMeta);
	    //inFilePtr ->GetObject("/AODBookkeeping/AODEvents", _hevents);
	    //_hevents->SetDirectory(0);
	    //_totalEvents += _hevents->Integral(1,1);//This is raw event number
	    //_totalMCatNLOEvents += _hevents->Integral(2,2);//This includes the MCatNLO weighting, used in EventWeight
	  }else{
	    //this is probably a data file, total events not used--> rel 16 v1 bug, diff name for histo
	    _totalEvents = 999;
	    _totalMCatNLOEvents = 999;
	  }

	  //cout<<"TOTAL AOD EVENTS "<<_totalEvents<<endl;
	  //Note- will break if file doesn't exist- if see a break here, check that file doesn't exist error hasn't also been thrown

	  // Increment count of lines in list 
	  if(iChain) ++chainLineCounter;
	  // Update total number of files in list (accounting for wildcard expansion)
	  chainFileCounter = chainFileCounter + iChain;

	  inFilePtr->Close();
	  delete inFilePtr;
	} // else: read in file
      } // while(getline)
      //      _hevents=NULL;
  
      cout<<"<AnalysisMain::ParseCmdLine> TOTAL AOD EVENT COUNT "<<_totalEvents<<endl;
      cout<<"<AnalysisMain::ParseCmdLine> TOTAL AOD EVENT COUNT WITH MCATNLO "<<_totalMCatNLOEvents<<endl;
      cout<<"<AnalysisMain::ParseCmdLine> TOTAL EVENT COUNT IN "<<chainEV0->GetName()<<" = "<<chainEV0->GetEntries()<<endl;
      // Print number of lines and files read in from list
      cout << "<AnalysisMain::ParseCmdLine> " << chainLineCounter << " lines read from chain list." << endl;
      cout << "<AnalysisMain::ParseCmdLine> " << chainFileCounter << " files successfully chained from chain list." << endl;
      
      // Exit if no lines or files were read from list
      if ( (0 == chainLineCounter) || (0 == chainFileCounter) ) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Empty list in file " << argv[fileListIndex] << endl;
	return 1;
      } //if

    }//if(chainFile)
    else{
      cout << "<AnalysisMain::ParseCmdLine> ERROR: Cannot find list file "<< argv[fileListIndex] << endl;
      return 1;
    }
  } //if(fileListFlag) 


  /////////////////////////////////////////////////////////////////////////////////
  // Use Single Root File
  /////////////////////////////////////////////////////////////////////////////////

  else if (fileSingleFlag) { 
    cout << "<AnalysisMain::ParseCmdLine> " << "Reading in single file: " << argv[fileSingleIndex] << endl;
    
    //had trouble coding error warning here, but ROOT will also throw out error message if TFile doesn't exist
    TFile *inFilePtr = new TFile(argv[fileSingleIndex]);
    
    // FastSim
    // If the event tree doesn't exist, then we do fast simulation
    if(! (inFilePtr -> Get(chainEV0->GetName()))) {
      cout<<"<AnalysisMain::ParseCmdLine> **WARNING** No event tree found. I'm going to assume that this is a fast simulation file."<<endl;
      SetDoFastSim(kTRUE);
    } //if

    // Truth
    // First check that command line parameter was passed
    if(DoTruth()){
      // If command line parameter passed check that Truth tree exists
      if(inFilePtr->Get(chainTruth->GetName())) {
        //cout<<"Found Truth object in file!"<<endl;
        SetDoTruth(kTRUE);
      } //if
      else {
        cout << "<AnalysisMain::ParseCmdLine> **WARNING** Truth tree not found. Running without it!" <<endl;
        SetDoTruth(kFALSE);
      } //else
    } //if
    else{
      cout << "<AnalysisMain::ParseCmdLine> **INFO** Truth parameter not passed. Not running over Truth tree." <<endl;
    } //if

     // Config
    if(inFilePtr -> Get(chainConfig -> GetName())) { SetDoConfig(kTRUE); }
    else {
      cout << "<AnalysisMain::ParseCmdLine> **WARNING** Config tree not found. Running without it!" <<endl;
      SetDoConfig(kFALSE);
    } //else

    /////////////////////////////////////////////////////////////////////////////////////////
    //Get total number of events from AOD
    //////////////////////////////////////////////////////////////////////////////////////  
    if(DoTruth()||DoSkim()) {
      inFilePtr->GetObject("/LooseTopInputsTreeMeta/Lumi",_skimLooseTopInputsTreeMeta);
      inFilePtr ->GetObject("/AODBookkeeping/AODEvents", _hevents);
      _hevents->SetDirectory(0);
      _totalEvents += _hevents->Integral(1,1);//This is raw event number
      _totalMCatNLOEvents += _hevents->Integral(2,2);//This includes the MCatNLO weighting, used in EventWeight
      //cout<<"TOTAL AOD EVENTS "<<_totalEvents<<endl;
      cout<<"<AnalysisMain::ParseCmdLine> TOTAL AOD EVENT COUNT "<<_totalEvents<<endl;
      cout<<"<AnalysisMain::ParseCmdLine> TOTAL AOD EVENT COUNT WITH MCATNLO "<<_totalMCatNLOEvents<<endl;
      //Note- will break if file doesn't exist- if see a break here, check that file doesn't exist error hasn't also been thrown
      
    }else{
      //this is probably a data file, total events not used--> rel 16 v1 bug, diff name for histo
      _totalEvents = 999;
      _totalMCatNLOEvents = 999;
    }
    
    // Close and delete file
    inFilePtr -> Close();
    delete inFilePtr;
	  
    // Add files to chains
    //iChain = chainCollection -> Add(argv[fileSingleIndex]);
    if(DoFastSim()) iChain = chainFastSim0   -> Add(argv[fileSingleIndex]);
    else            iChain = chainEV0        -> Add(argv[fileSingleIndex]);

    if(DoTruth())     chainTruth      -> Add(argv[fileSingleIndex]);
    if(DoConfig())    chainConfig     -> Add(argv[fileSingleIndex]);
    if (!iChain) {
      cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Chaining " << argv[fileSingleIndex];
      return 1;
    } //if
  } //if(fileSingleFlag)

  /////////////////////////////////////////////////////////////////////////////////
  // No chain or single file to input
  /////////////////////////////////////////////////////////////////////////////////
  else {
    cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: No chain or single file to input." << endl;
    return 1;
  } //else

 
  ////////////////////////////////////////////////////////////////////////////////
  // ???????????
  ////////////////////////////////////////////////////////////////////////////////
  // Create the root trees
  //CollectionTree* treeCollection = new CollectionTree(chainCollection);
  EventTree*      treeEventTree      = NULL;
  FastSimTree*    treeFastSimTree    = NULL;

  TruthTree*       treeTruthTree     = NULL;
  if(DoTruth())    treeTruthTree     = new TruthTree(chainTruth);

  if(DoFastSim()) {
    treeFastSimTree  = new FastSimTree(chainFastSim0);
  } //if
  else treeEventTree = new EventTree(chainEV0);

  // Open the output histogram file
  OpenHistogramFile(_histogramFileName);

  /////////////////////////////////////////////////////////////////////////////////
  // Parse the configuration file (part 2)
  //- nothing specified here related to commandline overrides (no object cut file)
  //- contains config setting which have default values
  /////////////////////////////////////////////////////////////////////////////////

  // Debug level
  if(_config.Defined("DebugLevel")) {
    SetDebugLevel(_config.GetValue("DebugLevel",0));
    if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Parsing configuration file (part2) " << endl;
    if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Setting debug level to " <<GetDebugLevel() << endl;
  } //if

  // Set cellout, softjet, and pileup systematic uncertainties.
  if ((celloutShift()!=0 ||softjetShift()!=0) && !_JESconfigread){
    cout << "<AnalysisMain::ParseCmdLine> " << "You have specified a cellout or softjet shift but no JES config file!  Your results will be invalid! " << endl;
  }

  // The target top mass used in some optimization (default: 172.5GeV), changed from 175 starting with 10TeV, v14 MC
  SetTargetTopMass(_config.GetValue("MTop",172.5));
  if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Setting target top mass to " << GetTargetTopMass() << endl;

  // Include a file containing MC total numbers
  if(_config.Defined("Include.MCTotEventNumberFile")) {
    TString MCevtFile = _config.GetValue("Include.MCTotEventNumberFile","NONE");
    // Check if this input file exists!
    FILE *evtFile;
    if(!(evtFile = fopen(MCevtFile.Data(), "r"))) {
      cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: MCTotEventNumberFile file "<<MCevtFile.Data()<<" does not exist!"<<endl;
      return -1;
    } //if
    fclose(evtFile);
    // The file exists, read it in:
    if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Reading cuts config file: " << MCevtFile.Data() << endl;    
    _config.ReadFile(MCevtFile.Data(),kEnvChange);
  } //if include cuts file

  // Include a file containing pile-up vtx reweighting weights
  if(_config.Defined("Include.pileUpReWeightingFile")) {
    TString MCevtFile = _config.GetValue("Include.pileUpReWeightingFile","NONE");
    // Check if this input file exists!
    FILE *evtFile;
    if(!(evtFile = fopen(MCevtFile.Data(), "r"))) {
      cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: pileUpReWeightingFile file "<<MCevtFile.Data()<<" does not exist!"<<endl;
      return -1;
    } //if
    fclose(evtFile);
    // The file exists, read it in:
    if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Reading cuts config file: " << MCevtFile.Data() << endl;    
    _config.ReadFile(MCevtFile.Data(),kEnvChange);
  } //if include cuts file

  // Include a file containing cuts
  if(_config.Defined("Include.CutsFile")) {
    TString cutsFile = _config.GetValue("Include.CutsFile","NONE");
    // Check if this input file exists!
    FILE *confFile;
    if(!(confFile = fopen(cutsFile.Data(), "r"))) {
      cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Cuts file "<<cutsFile.Data()<<" does not exist!"<<endl;
      return -1;
    } //if
    fclose(confFile);
    // The file exists, read it in:
    if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Reading cuts config file: " << cutsFile.Data() << endl;    
    _config.ReadFile(cutsFile.Data(),kEnvChange);
  } //if include cuts file

  // Include a file containing weights for each source
  if(_config.Defined("Include.WeightsFile")) {
    TString weightsFile=_config.GetValue("Include.WeightsFile","NONE");
    // check if this input file exists!
    FILE *weightFile;
    if(!(weightFile = fopen(weightsFile.Data(), "r"))) {
      cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Weights file " << weightsFile.Data()<<" does not exist!" << endl;
      return -1;
    } //if
    fclose(weightFile);
    // The file exists, read it in:
    if(GetDebugLevel() > 0) cout << "<AnalysisMain::ParseCmdLine> " << "Reading weights config file: " << weightsFile.Data() << endl;    
    _config.ReadFile(weightsFile.Data(), kEnvChange);
  } //if include weights file
  
   
  /////////////////////////////////////////////////////////////////////////////////
  // If we want to skim
  /////////////////////////////////////////////////////////////////////////////////
  if(DoSkim() && _skimFileName!="NONE") {
    // Open Skim File
    OpenSkimFile(_skimFileName);
    TDirectory *BSMTreeDir = _skimFile->mkdir("TNT","TNTD");
    BSMTreeDir->cd();
    //_skimConfigTreeMeta         = (TTree*)chainConfig-> CloneTree(0);
    //TDirectory *ConfigTreeMetaDir = _skimFile->mkdir("ConfigTreeMeta","ConfigTreeMetaD");
    //ConfigTreeMetaDir->Add(_skimConfigTreeMeta);

    // Create clones of chains
    //_skimCollectionTree = (TTree*)chainCollection -> CloneTree(0);

    if(DoFastSim()){
      _skimEventTree    = (TTree*)chainFastSim0 -> CloneTree(0);
    }else{
      _skimEventTree    = (TTree*)chainEV0      -> CloneTree(0);
    }
    cout << "Skimed chain name = " << chainEV0->GetName() << endl;
    cout << "Skimed tree name = " << _skimEventTree->GetName() << endl;
  } //if

  ////////////////////////////////////////////////////////////////////////////////
  // The chains were read in successfully. Now initialize the analysis classes
  ////////////////////////////////////////////////////////////////////////////////
  //if(DoFastSim()) EventContainer::InitializeFastSim(m_tool,treeFastSimTree, treeTruthTree, treeTriggerTree, treeInfoTree, treeDecisionTree, treeElectronPreTagTree, treeDiElectronsPreTagTree, treeElectronMuonPreTagTree, treeMuonPreTagTree, treeDiMuonsPreTagTree, treeElectronLooseTree, treeDiElectronsLooseTree, treeElectronMuonLooseTree, treeMuonLooseTree, treeDiMuonsLooseTree, treeBunchConfTree);
  //else            EventContainer::Initialize(m_tool, treeEventTree, treeTruthTree, treeTriggerTree, treeInfoTree, treeDecisionTree, treeElectronPreTagTree, treeDiElectronsPreTagTree, treeElectronMuonPreTagTree, treeMuonPreTagTree, treeDiMuonsPreTagTree, treeElectronLooseTree, treeDiElectronsLooseTree, treeElectronMuonLooseTree, treeMuonLooseTree, treeDiMuonsLooseTree,treeBunchConfTree );
  EventContainer::Initialize(treeEventTree, treeTruthTree);

  EventContainer::SetupObjectDefinitions();

  EventContainer::SetEventCount(evtOffSet);

  return 0;
} //end AnalysisMain::ParseCmdLine

/******************************************************************************
 * void AnalysisMain::Loop()                                                  *
 *                                                                            *
 * Loops over Events                                                          *
 * Main Action of Program                                                     *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
 void AnalysisMain::Loop() {

  // Add histogram to _histogramFile directory in memory
  SetDirectory (_histogramFile);
  // Add this information to the additional variable processor
  SetVariableHistogramDirectory (_histogramFile);
  SetVariableEventContainer(GetEventContainer());
  //Set the number of cuts
  SetNCuts(GetNumberOfCuts());
  // Book Histograms
  cout << "<AnalysisMain::Loop::BookHistogram Start> " << endl;
  CutListProcessor::BookHistogram();
  cout << "<AnalysisMain::Loop::BookHistogram End> " << endl;

  // Set flag if an event limit is specified (_eventLimit read from cmd line)
  Bool_t eventLimitFlag = (_eventLimit > 0);

  // Start the stopwatch
  TStopwatch stopwatch;
  stopwatch.Start();

  // Get Event count which is determined from position in chain
  Int_t eventInChain = GetNextEvent();

  // Initialize Event count which is determined by looping over events
  Int_t eventCounter = 0;
  
  // Initialize flags for writing events
  Bool_t writeThisEvent    = kFALSE;
  Bool_t firstEventWritten = kFALSE;

  //  if (_skimEventTree != NULL){
  // AdditionalVarsProcessor::BookBranches(_skimEventTree);
  //}

  // Loop over events in chain until (1) Out of events or (2) Have reached limit specified on command line
  //cout << "<AnalysisMain::Loop> " <<  "Beginning loop over events, _numEvent is "<<_eventLimit << ", eventInChain is "<<eventInChain<<endl;
  cout << "<AnalysisMain::Loop> " << endl;

  while( eventInChain != -1 && (!eventLimitFlag || eventInChain < _eventLimit) ) {  
    eventCounter++;
    if( 0 == (eventInChain%10000) ) cout << "<AnalysisMain::Loop> " << "Processing event " << eventInChain << endl;

    //If we're only doing specific events, check them here.
    if (_eventsToRunOn.size() > 0){
      bool skipEvent = true;
      for (auto const eventID : _eventsToRunOn){
	if (eventID == this->eventNumber){
	  skipEvent = false;
	  break;
	}
      }
      if (skipEvent) {
	eventInChain = GetNextEvent(); 
	continue;
      }
      std::cout << "Running on event: " << this->eventNumber << std::endl;
    }

    // Fill the histograms
    writeThisEvent = CutListProcessor::Apply(this);
    // Fill Trees for Skim Events
    if(DoSkim() && writeThisEvent) {
      if (!firstEventWritten) {
	firstEventWritten = kTRUE;
	if( NULL != _skimEventTree) {
	  AdditionalVarsProcessor::BookBranches(_skimEventTree);
	  _HFORb = -999;
	  _newBranchHFORb = _skimEventTree->Branch("hfor_type", &_HFORb, "hfor_type/I");
	  _EventNBeforePreselb = -999;
	  _newBranchNBPb                 = _skimEventTree->Branch("EventNBeforePresel",   &_EventNBeforePreselb,   "EventNBeforePresel/F");

	}
      } //if
    
      if(NULL != _skimEventTree) {
	_EventNBeforePreselb   = _totalMCatNLOEvents;
	AdditionalVarsProcessor::ResetBranches();
	AdditionalVarsProcessor::FillBranches(this);
	//AdditionalVarsProcessor::OutputBranches();
	_skimEventTree          -> Fill();
      }
    } //if

    // Read the next event in the chain
    //cout << "<AnalysisMain::Loop> " << "Getting next event in chain " << eventInChain<<endl;
    // cout<<GetHFOR()<<endl;
    eventInChain = GetNextEvent();
  } // while
  

  /******************************************************************************
  // Print everything!
  ******************************************************************************/
  
  PrintCutDetails();

  /******************************************************************************
  // Stop the timer and print out timing statistics
  ******************************************************************************/
  stopwatch.Stop();
  Double_t timeReal = stopwatch.RealTime();
  Double_t timeEvent = timeReal / static_cast<Double_t>(eventCounter);

  cout << "<AnalysisMain::Loop> " << endl;
  cout << "<AnalysisMain::Loop> " <<"           ------------------------------- "  << endl;
  cout << "<AnalysisMain::Loop> " <<"          |      Processing Summary        |"  << endl;
  cout << "<AnalysisMain::Loop> " <<"           ------------------------------- "  << endl;
  cout << "<AnalysisMain::Loop> " <<"          |  # of events:    " << setw(10)    << setprecision(4) <<eventCounter   << "  " << " |"  << endl;

  cout << "<AnalysisMain::Loop> " <<"          |  Real time:      " << setw(10)    << setprecision(4) << timeReal      << " s" << " |"  << endl;
  cout << "<AnalysisMain::Loop> " <<"          |  Time per event: " << setw(10)    << setprecision(4) << timeEvent     << " s" << " |"  << endl;
  cout << "<AnalysisMain::Loop> " <<"          |  Events/second:  " << setw(10)    << setprecision(4) << 1.0/timeEvent << "  " << " |"  << endl;
  cout << "<AnalysisMain::Loop> " <<"           ------------------------------- "  << endl;
  cout << "<AnalysisMain::Loop> " << endl;

  // Write and Close histogram file
  // Write output to files
  cout << "<AnalysisMain::Loop> " << "Closing histogram file: " << GetHistogramFileName() << endl;
  WriteHistogramFile();
  CloseHistogramFile();
  
  // If applicable Write and Close skim file
  if(_doSkim) {
    cout << "<AnalysisMain::Loop> " << "Write skim file: " << GetSkimFileName() << endl;
    WriteSkimFile();
    cout << "<AnalysisMain::Loop> " << "Closing skim file: " << GetSkimFileName() << endl;
    CloseSkimFile();
  } //if

} // void AnalysisMain::Loop()

