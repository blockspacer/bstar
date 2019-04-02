/******************************************************************************
 * Wt.cpp                                                                     *
 *                                                                            *
 * Top Level file in analysis package                                         *
 * Drives the package                                                         *
 *                                                                            *
 * This file is the main file which will loop over all of the events, make    *
 * cuts and generate histograms.                                              *
 *                                                                            *
 * For general information on the analysis package as a whole, see            *
 * https://hep.pa.msu.edu/twiki/bin/view/AtlasSingleTop/HowToUseTheAnalysisFramework *
 *                                                                            *
 * Modifications:                                                             *
 *   23 Nov 2010 - H. ZHANG, modify first version from SandB analysis         *
 *****************************************************************************/
#define Glodfsdbal_cxx
#include "SingleTopRootAnalysis/Base/Dictionary/AnalysisMain.hpp"
#include <iostream>
#include <string>
#include "TSystem.h"

// Include histogramming classes
#include "SingleTopRootAnalysis/Cuts/Weights/EventWeight.hpp"
//#include "SingleTopRootAnalysis/Histogramming/Topological/HistogrammingWtDiLepTopology.hpp"
#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingMuon.hpp"
#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingElectron.hpp"
#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingMET.hpp"
#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingMtW.hpp"
#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingJetAngular.hpp"
#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingJet.hpp"
#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingNPvtx.hpp"
// Include cuts classes
//#include "SingleTopRootAnalysis/Cuts/Other/CutTriangularSumDeltaPhiLepMET.hpp"
//#include "SingleTopRootAnalysis/Cuts/Other/CutEMuOverlap.hpp"
#include "SingleTopRootAnalysis/Cuts/Jet/CutJetN.hpp"
#include "SingleTopRootAnalysis/Cuts/TaggedJet/CutTaggedJetN.hpp"
//#include "SingleTopRootAnalysis/Cuts/Jet/CutJetPt1.hpp"

//#include "SingleTopRootAnalysis/Cuts/Other/CutMissingEt.hpp"
#include "SingleTopRootAnalysis/Cuts/Muon/CutMuonN.hpp"
//#include "SingleTopRootAnalysis/Cuts/Muon/CutMuonTighterPt.hpp"
#include "SingleTopRootAnalysis/Cuts/Electron/CutElectronN.hpp"
//#include "SingleTopRootAnalysis/Cuts/Electron/CutElectronTighterPt.hpp"
//#include "SingleTopRootAnalysis/Cuts/Other/CutHTJET1.hpp"
//#include "SingleTopRootAnalysis/Cuts/Lepton/CutLeptonOppositeCharge.hpp"
//#include "SingleTopRootAnalysis/Cuts/Other/CutLarBurst.hpp"
#include "SingleTopRootAnalysis/Cuts/Other/CutPrimaryVertex.hpp"
//
//#include "SingleTopRootAnalysis/Cuts/Other/CutZveto.hpp"
#include "SingleTopRootAnalysis/Cuts/Other/CutTriggerSelection.hpp"
//#include "SingleTopRootAnalysis/Cuts/Other/CutHFOR.hpp"
//#include "SingleTopRootAnalysis/Cuts/Other/CutBadPileupEvent.hpp"

//Include the additional variables that we would want in the skim tree
#include "SingleTopRootAnalysis/Vars/TestVar.hpp"
#include "SingleTopRootAnalysis/Vars/BDTVars.hpp"
#include "SingleTopRootAnalysis/Vars/WeightVars.hpp"
#include "SingleTopRootAnalysis/Vars/ChannelFlag.hpp"

using std::cout;
using std::endl;

/******************************************************************************
 * int main(int argc, char **argv)                                            *
 *                                                                            *
 * Main Routine                                                               *
 *                                                                            *
 * Input:  Command line arguements                                            *
 * Output: Returns 0 of successful                                            *
 ******************************************************************************/
// Only declare main if not in CINT (not interactive?)
#ifndef __CINT__
int main(int argc, char **argv) 
#endif
{
  cout << "<driver> Starting parton_analysis" << endl;
  gROOT->ProcessLine("#include <vector>");

  /////////////////////////////////////////////////////////////////////////////////
  // Decleartions and Instantiations
  /////////////////////////////////////////////////////////////////////////////////
  // Instantiate the analysis class
  AnalysisMain mystudy;

  TChain *chainReco       = new TChain("TNT/BOOM");
  TChain *chainTruth      = new TChain("TruthTree");
  TChain *chainTrigger    = new TChain("TriggerTree");
  TChain *chainInfo       = new TChain("InfoTree");
  TChain *chainDecision   = new TChain("DecisionTree");
  TChain *chainFastSim    = new TChain("FastSimTree");

 
 // Check command line for Fast Running
  // In Fast Running only a few cuts are made and histograms filled
  Bool_t doFast = kFALSE;
  string mcStr="";
  Bool_t doMC = kFALSE;
  Bool_t doPileup = kFALSE;
  Bool_t dobWeight = kFALSE;
  Bool_t useInvertedIsolation = kFALSE;
  Bool_t useLeptonSFs = kFALSE;
  Bool_t usebTagReweight = kFALSE;
  TString leptonTypeToSelect = "Tight"; //This variable is altered to unisolated when running QCD estimation.
  string evtListFileName="";
  int whichtrig = -1;
  // A couple of jet selection overrides to limit the number of config faces.  
  // -1 means use the value from the config file
  Int_t nJets = -1;
  Int_t nbJets = -1;

  //  cout << "<Driver> Reading arguments" << end;
  for (int i = 1; i < argc; ++i) {
    cout<<"Command line parameter "<<i<<" is "<<argv[i]<<endl;
    if (!strcmp(argv[i], "-inlist")) {
	string ListName = (string)(argv[i+1]);
	int b = ListName.find("SingleTop.");
	//	evtListFileName="WtDiElectron_"+ListName.substr(b)+".lst.root";
	//cout << "evtListFileName : " << evtListFileName << endl;
    }
    if (!strcmp(argv[i], "-lepSFs")){
      useLeptonSFs = kTRUE;
      cout << "Driver: Using lepton ID/iso SFs" << endl;
    }
    if (!strcmp(argv[i], "-bTagReshape")){
      usebTagReweight = kTRUE;
      cout << "Driver: Enable b-tag discriminant reshaping" << endl;
    }
    if (!strcmp(argv[i], "-fast")) {
      doFast = kTRUE;
      cout << "Driver: Fast running " << endl;
    }//if fast
    if (!strcmp(argv[i], "-MC") ||!strcmp(argv[i], "-mc") ) {
      doMC = kTRUE;
      cout << "Driver: This is a MC sample" << endl;
    }//if MC
    if (!strcmp(argv[i], "-MCatNLO")) {
      mcStr=mcStr+"MCatNLO";
      cout << "Driver: This is MCatNLO " << endl;
      //Number of total events automatically includes the number with MCatNLO weights, but the events themselves will not be weighted in current code version unless this flag is called
    }//if MCatNLO
    if (!strcmp(argv[i], "-PileUpWgt")) {
      mcStr=mcStr+"PileUpWgt";
      doPileup = kTRUE;
      cout << "Driver: Use PileUpWgt " << endl;
    }//if PileUpWgt
    if (!strcmp(argv[i], "-BWgt")) {
      mcStr=mcStr+"BWgt";
      dobWeight = kTRUE;
      cout << "Driver: Use b weight " << endl;
    }//if PileUpWgt
    if (!strcmp(argv[i], "-UseTotalEvtFromFile")) {
      mcStr=mcStr+"UseTotalEvtFromFile";
      cout << "Driver: UseTotalEvtFromFile " << endl;
    }//if UseTotalEvtFromFile
    if(!strcmp(argv[i],"-InvertIsolation")) {
      useInvertedIsolation = kTRUE;
      leptonTypeToSelect = "UnIsolated";
      cout << "Driver: useInvertedIsolation " << endl;
    }
    if(!strcmp(argv[i],"-SelectTrigger")) {
      if (argc < i+1) {
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for SelectTrigger - must choose either Electron or Muon" << endl;
	return 1;
      }
 
      if(!strcmp(argv[i+1],"Electron")) whichtrig = 0;
      if(!strcmp(argv[i+1],"Muon")) whichtrig = 1;
      if( whichtrig == -1){
	cout << "<AnalysisMain::ParseCmdLine> " << "ERROR: Incorrect Value for SelectTrigger - must choose either Electron or Muon" << endl;
	return 1;
      }
    }

    if (!strcmp(argv[i],"-BkgdTreeName")){
      // Check if BkgdTree is specified.  
      if (argc < i+1) {
	cout << " <AnalysisMain::ParseCmdLine> " << "ERROR: Missing Value for tree name for multivariate variables -BkgdTree!" << endl;
	return 1;
      } //if
      mystudy.SetBkgdTreeName(argv[i+1]);
    }// if BkgdTreeName
    if (!strcmp(argv[i],"-nJets")){
      nJets = atoi(argv[i+1]);
    }// if nJets
    if (!strcmp(argv[i],"-nbJets")){
      nbJets = atoi(argv[i+1]);
    }// if nbJets
  } // for

  TChain *chainBkgd = new TChain(mystudy.GetBkgdTreeName());//if nothing specified on commandline, this is just a "" string, which is ok if you don't want to use this chain (don't want a tree of multivariate variables).

  //////////////////////////////////////////////////////////////////////////////////
  // Parse the Command Line (the rest of it)
  //////////////////////////////////////////////////////////////////////////////////
  cout<<"USING THE TOPINPUTS TREE"<<endl;
  Int_t parseReturnValue = mystudy.ParseCmdLine(argc,argv,chainReco,chainTruth,chainTrigger,chainInfo,chainDecision, chainFastSim, chainBkgd);
  cout<<"ParseCmdLine done"<<endl;
  if(parseReturnValue != 0) return parseReturnValue;

  /////////////////////////////////////////////////////////////////////////////////
  // Define a Particle class to be the same as the AnalysisMain class
  /////////////////////////////////////////////////////////////////////////////////
  EventContainer *particlesObj = &mystudy;
  particlesObj->SetIsSimulation(doMC);
  TString BkgdTreeName=mystudy.GetBkgdTreeName();
  bool isee   = (BkgdTreeName == "DiElectronsPreTagTree"||BkgdTreeName=="DiElectronsLooseTree");
  bool ismumu = (BkgdTreeName == "DiMuonsPreTagTree"||BkgdTreeName=="DiMuonsLooseTree");
  bool isemu  = (BkgdTreeName == "ElectronMuonPreTagTree" ||  BkgdTreeName == "ElectronMuonLooseTree");
  particlesObj->SetUseUnisolatedLeptons(useInvertedIsolation,whichtrig);

  /////////////////////////////////////////////////////////////////////////////////
  // Add Cuts and Histograms applicable to Fast and Full Analyses
  /////////////////////////////////////////////////////////////////////////////////
  // ******** Cuts and Histograms applied to all studies ********

  mystudy.AddCut(new EventWeight(particlesObj,mystudy.GetTotalMCatNLOEvents(), mcStr, doPileup, dobWeight, useLeptonSFs, usebTagReweight));

  mystudy.AddCut(new HistogrammingMuon(particlesObj,"All"));  // make the muon plots, hopefully.
  mystudy.AddCut(new HistogrammingMuon(particlesObj,"Tight"));  // make the muon plots, hopefully.
  mystudy.AddCut(new HistogrammingMuon(particlesObj,"Veto"));  // make the muon plots, hopefully.
  mystudy.AddCut(new HistogrammingMuon(particlesObj,"UnIsolated"));  // make the muon plots, hopefully.
  mystudy.AddCut(new CutPrimaryVertex(particlesObj));
  mystudy.AddCut(new CutTriggerSelection(particlesObj, whichtrig));

  mystudy.AddCut(new HistogrammingMET(particlesObj));
  //mystudy.AddCut(new CutElectronTighterPt(particlesObj, "Tight")); 
  mystudy.AddCut(new CutMuonN(particlesObj, leptonTypeToSelect));     //require that lepton to be isolated, central, high pt

  mystudy.AddCut(new CutMuonN(particlesObj, "Veto"));     //require that lepton to be isolated, central, high pt
 


  mystudy.AddCut(new CutElectronN(particlesObj, leptonTypeToSelect)); //require that lepton to be isolated, central, high pt
  mystudy.AddCut(new CutElectronN(particlesObj, "Veto")); //require that lepton to be isolated, central, high pt

  mystudy.AddCut(new HistogrammingElectron(particlesObj,leptonTypeToSelect));  // make the muon plots, hopefully.
  mystudy.AddCut(new HistogrammingElectron(particlesObj,"Veto"));  // make the muon plots, hopefully.

  mystudy.AddCut(new HistogrammingMET(particlesObj));
  mystudy.AddCut(new HistogrammingMtW(particlesObj,useInvertedIsolation));

  mystudy.AddCut(new HistogrammingMuon(particlesObj,leptonTypeToSelect));  // make the muon plots, hopefully.

  //mystudy.AddCut(new CutMuonTighterPt(particlesObj, "Tight")); //require that new Pt cut for leading and subleading muon  
  //if (isemu){
  //  mystudy.AddCut(new CutEMuOverlap(particlesObj));
  //}
  //mystudy.AddCut(new CutJetPt1(particlesObj));
  mystudy.AddCut(new CutJetN(particlesObj,nJets));
  
  mystudy.AddCut(new CutTaggedJetN(particlesObj,nbJets));

  mystudy.AddCut(new EventWeight(particlesObj,mystudy.GetTotalMCatNLOEvents(), mcStr, doPileup, dobWeight, useLeptonSFs, usebTagReweight));

  mystudy.AddCut(new HistogrammingMuon(particlesObj,"Tight"));  // make the muon plots, hopefully.

  mystudy.AddCut(new HistogrammingMET(particlesObj));
  mystudy.AddCut(new HistogrammingMtW(particlesObj,useInvertedIsolation));
  mystudy.AddCut(new HistogrammingJetAngular(particlesObj,useInvertedIsolation));
  mystudy.AddCut(new HistogrammingJet(particlesObj));
  mystudy.AddCut(new HistogrammingNPvtx(particlesObj));
  //mystudy.AddCut(new CutTriangularSumDeltaPhiLepMET(particlesObj));  
  //if (isemu){
  //  mystudy.AddCut(new CutHTJET1(particlesObj));
  //}
  //mystudy.AddCut(new CutMissingEt(particlesObj));
  //mystudy.AddCut(new CutLeptonOppositeCharge(particlesObj));
  //if (isee || ismumu){
  //  mystudy.AddCut(new CutZveto(particlesObj, "Tight"));
  //}

  //Add in any variables to the skim tree that you want here
  //  mystudy.AddVars(new TestVar());
  //if (whichtrig) mystudy.AddVars(new BDTVars(true));
  mystudy.AddVars(new BDTVars(true));
  mystudy.AddVars(new WeightVars());
  mystudy.AddVars(new ChannelFlag());

  TFile *_skimBDTFile;
  TString NNname = mystudy.GetHistogramFileName() + "skimBDT.root" ;
  _skimBDTFile = new TFile(NNname,"RECREATE"); 
  //     mystudy.AddCut(new HistogrammingWtDiLepTopology(particlesObj, _skimBDTFile));
  //mystudy.AddCut(new HistogrammingWtDiLepTopology(particlesObj, _skimBDTFile));  // used to select a Zee candidate sample
  
  cout << "<driver> Starting Loop over events in chain" << endl << endl;
  // Loop over events in chain
  //Here is where the events are actually looped and the cuts and histograms above are made
  mystudy.Loop();
  cout << "<driver> Finished Loop over events in chain" << endl;

  // Write output to files
  cout << "<driver> " << "Writing histograms to file: " << mystudy.GetHistogramFileName() << endl;
  
  _skimBDTFile->Write(); 
  _skimBDTFile->Close(); 

  mystudy.Finish();
  
  // Delete chains and set pointers to NULL
  cout << "<driver> Deleting chains" << endl;

  if(chainReco)       delete chainReco;
  if(chainTruth)      delete chainTruth;
  if(chainFastSim)    delete chainFastSim;

  //if(chainCollection) delete chainCollection;

  chainReco       = NULL;
  chainTruth      = NULL;
  chainFastSim    = NULL;
  
  cout << "<driver> Sucessful Completion." << endl;
  return 0;
}
