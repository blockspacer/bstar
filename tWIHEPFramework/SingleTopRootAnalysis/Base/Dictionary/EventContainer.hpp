/******************************************************************************
 * EventContainer.hpp                                                         *
 *                                                                            *
 * Read in variables from EventTree and put them into an Object               *
 *                                                                            *
 * Public Members of EventContainer    class                                  *
 *    EventContainer()                  -- Default Constructor                *
 *    ~EventContainer()                 -- Destructor                         *
 *    GetSourceName()                   -- Get name of source being run on    *
 *    GetSourceNumber()                 -- Get number of source being run on  *
 *    GetTargetTopMass()                -- Get Top Quark Mass                 *
 *    SetSourceName()                   -- Set name of source being run on    *
 *    SetSourceNumber()                 -- Set number of source being run on  *
 *    SetTargetTopMass()                -- Set Top Quark Mass                 *
 *    Initialize()                      -- Initialize class                   *
 *    MakeTopQuarks()                   -- Reconstruct Top in several ways    *
 *    GetNextEvent()                    -- Set up particles for next event    *
 *    SetEventCount()                   -- Set counter for specific event     *
 *    GetCollectionTree()               -- Get the collection tree            *
 *    GetEventTree()                    -- Get the event tree                 *
 *    GetTruthTree()                    -- Get the truth tree                 *
 *    GetTruthAllTree()                    -- Get the truth tree              *
 *    GetTriggerTree()                  -- Get the trigger tree               *
 *    GetFastSimTree()                  -- Get the fastsim tree               *
 *    GetDebugLevel()                   -- Get debug level                    *
 *    SetDebugLevel()                   -- Set debug level                    *
 *    SetDoTrigger()                    -- Set if trigger tree                *
 *    SetDoTruthAll()                   -- Set if TruthAll tree               *
 *    SetDoTruth()                      -- Set if Truth tree                  *
 *    DoTrigger()                       -- If trigger tree                    *
 *    DoTruthAll()                      -- If TruthAll tree                   *
 *    DoTruth()                         -- If Truth tree                      *
 *    GetEventWeight()                  -- Get Event Weight                   *
 *    GetOutputEventWeight()            -- Get Output EventWeight             *
 *    GetTreeEventWeight()              -- Get Tree Event Weight              *
 *    GetGlobalEventWeight()            -- Get Global EventWeight             *
 *    SetGlobalEventWeight()            -- Set global event weight            *
 *    SetOutputEventWeight()            -- Set output event weight            *
 *    GetConfig()                       -- Get Configuration                  *
 *                                                                            *
 *    Int_t ReadEvent()                 -- Create TLorentz Vectors            *
 *    Bool_t CreateFinalState()         -- Create final state particles       *
 *    Int_t runNumber                   -- Run number from the root tree      *
 *    Int_t eventNumber                 -- Event number from the root tree    *
 *    vector<Electron> electrons        -- Electron Vector                    *
 *    vector<Muon>     muons            -- Muon Vector                        *
 *    vector<Muon>     isolatedmuons    -- Isolated Muon Vector               *
 *    vector<Muon>     unisolatedmuons  -- UnIsolated Muon Vector             *
 *    vector<Tau>      taus             -- Tau Vector                         *
 *    vector<JetOR>      jetors         -- JetOR Vector (special overlap) vec *
 *    vector<Jet>      jets             -- Jet Vector                         *
 *    vector<Jet>      taggedJets       -- taggedJet Vector                   *
 *    vector<Jet>      untaggedJets     -- untaggedJet Vector                 *
 *    vector<MCParticle> MCParticles    -- MC Particle Vector                 *
 *    vector<MCMuon>   MCmuons          -- MC Muon Vector                     *
 *    vector<MCElectron>   MCElectrons  -- MC Electron Vector                 *
 *    vector<MCTau>   MCTaus            -- MC Tau Vector                      *
 *    vector<MCJet>   MCJets            -- MC Jet Vector                      *
 *    vector<MCJet>   MCBJets           -- MC BJet Vector                     *
 *    vector<MCJet>   MCCJets           -- MC CJet Vector                     *
 *    vector<MCJet>   MCLightQuarkJets  -- MC LightJetJet Vector              *
 *    vector<MCTop>   MCTops            -- MC Top Vector                      *
 *    vector<MCW>     MCWs              -- MC W Vector                        * 
 *    vector<MCAllParticle> MCAllParticles -- MCAllParticle Vector            *
 *    vector<MCAllMuon> MCAllMuons      -- MCAllMuon Vector                   *
 *    vector<MCAllElectron> MCAllElectrons -- MCAll Electron Vector           *
 *    vector<MCAllJet> MCAllJets         -- MCAll Jet Vector                  *
 *    vector<MCAllTau> MCAllTaus         -- MCAll Tau Vector                  *
 *    vector<MCAllChargedPion> MCAllChargedPions -- MCAll Charged Pion Vector *
 *    vector<MCAllChargedK> MCAllChargedKs -- MCAll Charged K Vector          *
 *                                                                            *
 * Protected Members of EventContainer class                                  *
 *    Int_t debugLevel                  -- Debug level                        *
 *    TEnv _config                      -- Used to read configuration         *
 *                                                                            *
 * Private Members of EventContainer class                                    *
 *    Double_t _targetTopMass           -- Target top mass                    *
 *    CollectionTree *_collectionTree   -- Collection tree                    *
 *    EventTree *_eventTree             -- Event tree                         *
 *    TruthTree *_thruthTree            -- Truth tree                         *
 *    TruthTree *_thruthAllTree         -- TruthAll tree                      *
 *    Int_t _eventCount                 -- Location in event chain            *
 *    TString _sourceName               -- Name of source being read          *
 *    TString _sourceNumber             -- Number of source being read        *
 *    Double_t _globalEventWeight       -- Global Event Weight                *
 *    Double_t _treeEventWeight         -- Tree Event Weight                  *
 *    Double_t _outputEventWeight       -- Output Event Weight                *
 *                                                                            *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutJetN"                  *
 *                                                                            *
 * Private Data Members of CutJetN                                            *
 *    myTH1F* _hNJetsBefore              -- Hist of Num of jets before cuts   *
 *    myTH1F* _hNJetsAfter               -- Hist of Num of jets before cuts   *
 *    Int_t _nJetsMin;                   -- Minimum number of jets to require *
 *    Int_t _nJetsMax;                   -- Maximum number of jets allowed    *
 *                                                                            *
 * History                                                                    *
 *      17 June 2015   -   Huaqiao ZHANG Port from STAR                       *
 *****************************************************************************/

#ifndef particles_h
#define particles_h

#include <vector>
#include <string>
#include <iostream>

// MakeClass output:
#include "SingleTopRootAnalysis/Trees/EventTree.hpp"
#include "SingleTopRootAnalysis/Trees/FastSimTree.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

// Our own particles
#include "SingleTopRootAnalysis/Particles/Recon/Electron.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Muon.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Tau.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Jet.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/BoostedJet.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Neutrino.hpp"

// MC particles
#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCMuon.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCElectron.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCTau.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCJet.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCTop.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCW.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCNeutrino.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCBQuark.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCCQuark.hpp"

// A global CutFlowTable to keep track of individual cuts
#include "SingleTopRootAnalysis/Base/Dictionary/CutFlowTable.hpp"

#include "SingleTopRootAnalysis/Base/Dictionary/Top_MET.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/MET_Cleaning_Utils.hpp"


#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TEnv.h>
#include <TMath.h>
#include <TString.h>

//using namespace Analysis;

class EventContainer
{
 public:
  //Default Constructor
  EventContainer();

  //Destructor
  ~EventContainer();

  // Get and Set the NAME of the source being run over
  inline void SetSourceName(const TString& name) { _sourceName = name; };
  inline TString GetSourceName() const { return _sourceName; };

  // Get and Set the NUMBER for the source being run over
  inline void SetSourceNumber(Int_t number) { _sourceNumber = number; };
  inline Int_t GetSourceNumber() const { return _sourceNumber; };

  // Get and Set the target top mass for top quark constraint in best jet algorithm
  inline void SetTargetTopMass(const Double_t& mass) { _targetTopMass = mass; };
  inline Double_t GetTargetTopMass() const { return _targetTopMass; };

  // Initialize this class 
  void Initialize(EventTree* eventTree, TruthTree* truthTree);
 
  //Set up the object definitions
  void SetupObjectDefinitions();

  //void InitializeFastSim(Trig::D3PDBunchCrossingToolSA* m_tool,  FastSimTree* fsTree, TruthTree* truthTree, TriggerTree* triggerTree, InfoTree* infoTree, DecisionTree* decisionTree, ElectronPreTagTree* ElectronPreTagTree, DiElectronsPreTagTree* DiElectronsPreTagTree, ElectronMuonPreTagTree* ElectronMuonPreTagTree, MuonPreTagTree* MuonPreTagTree, DiMuonsPreTagTree* DiMuonsPreTagTree, ElectronLooseTree* ElectronLooseTree, DiElectronsLooseTree* DiElectronsLooseTree, ElectronMuonLooseTree* ElectronMuonLooseTree, MuonLooseTree* MuonLooseTree, DiMuonsLooseTree* DiMuonsLooseTree, BunchConfTree* treeBunchConfTree );

  // Create TLorentz Vectors for the particles in the interaction
  Int_t ReadEvent();

  // Create final state particles and make cuts if necessary
  Bool_t CreateFinalState();

  // Reconstruct the top in several different ways from the final state objects
  void MakeTopQuarks();

  // Set up the particles for the next event
  Int_t GetNextEvent();

  // Initialize the event counter to get a specific event in the chain
  inline void SetEventCount(Int_t i) { _eventCount = i; };

  // Get the collection tree
  //inline CollectionTree* GetCollectionTree() const { return _collectionTree; };

  // Get the event tree
  inline EventTree* GetEventTree() const { return _eventTree; };

  // Get the FastSim tree
  inline FastSimTree* GetFastSimTree() const { return _fastSimTree; };

  // Get the truth tree
  inline TruthTree* GetTruthTree() const { return _truthTree; };

    // Get and Set Debug Level
  inline Int_t GetDebugLevel() const { return _debugLevel; };
  inline void SetDebugLevel(Int_t i) { _debugLevel=i; };

  // Set whether we are processing a fastsim file or a full sim tree
  inline void SetDoFastSim(const Bool_t& a=true) { _doFastSim=a; };
  // get whether we should do the fastsim tree or not
  inline Bool_t DoFastSim(void) const { return _doFastSim; };

  // Set whether we should do the Truth
  inline void SetDoTruth(const Bool_t& a=kTRUE) { _doTruth=a; };
  // get whether we should do the Truth tree or not
  inline Bool_t DoTruth(void) const { return _doTruth; };

  // get whether we should do the skim tree or not
  inline void SetDoSkim(const Bool_t& a=true) { _doSkim=a; };
  // set whether we should do the skim tree or not
  inline Bool_t DoSkim(void) const { return _doSkim; };

 // Set whether we should do the Config tree
  inline void SetDoConfig(const Bool_t& aa=true) { _doConfig=aa; };
  // get whether we should do the Config tree or not
  inline Bool_t DoConfig(void) const { return _doConfig; };

  void SetUseUnisolatedLeptons(const Bool_t& a=true, int whichtrig=1);
  inline Bool_t UseUnisolatedLeptons(void) const {return _useUnisolatedLeptons; };

  // Set the global event weight
  inline void SetGlobalEventWeight(const Double_t& wt) { _globalEventWeight=wt; };
  // set the output event weight.
  inline void SetOutputEventWeight(const Double_t& wt) { _outputEventWeight=wt; };

  // Set the jet energy scale error
  inline void SetJESError(const Double_t& er) { _jesError=er; };
  // Get the jet energy scale error
  inline Double_t GetJESError(void) const { return _jesError; };
  inline Double_t JESError(void) const { return GetJESError(); };
  // 
  // set by how much we should shift the jet energy on this run.
  // this is either +1 or -1 or 0 (no shift)
  inline void SetJESShift(Int_t is) { _jesShift=is; };
  inline Int_t GetJESShift() const { return _jesShift; };
  inline Int_t JESShift() const { return GetJESShift(); };

  inline void SetJESUShift(Int_t is) { _jesUShift=is; };
  inline Int_t GetJESUShift() const { return _jesUShift; };
  inline Int_t JESUShift() const { return GetJESUShift(); };

  inline void SetJESPtShift(Int_t is) { _jesPtShift=is; };
  inline Int_t GetJESPtShift() const { return _jesPtShift; };
  inline Int_t JESPtShift() const { return GetJESPtShift(); };

  inline void SetJESLowPtShift(Int_t is) { _jesLowPtShift=is; };
  inline Int_t GetJESLowPtShift() const { return _jesLowPtShift; };
  inline Int_t JESLowPtShift() const { return GetJESLowPtShift(); };

  inline void SetJESHighPtShift(Int_t is) { _jesHighPtShift=is; };
  inline Int_t GetJESHighPtShift() const { return _jesHighPtShift; };
  inline Int_t JESHighPtShift() const { return GetJESHighPtShift(); };

  inline void SetJESEtaShift(Int_t is) { _jesEtaShift=is; };
  inline Int_t GetJESEtaShift() const { return _jesEtaShift; };
  inline Int_t JESEtaShift() const { return GetJESEtaShift(); };

  inline void SetcelloutShift(Double_t is) { _celloutShift=is; };
  inline Double_t GetcelloutShift() const { return _celloutShift; };
  inline Double_t celloutShift() const { return GetcelloutShift(); };

  inline void SetsoftjetShift(Double_t is) { _softjetShift=is; };
  inline Double_t GetsoftjetShift() const { return _softjetShift; };
  inline Double_t softjetShift() const { return GetsoftjetShift(); };

  inline void SetpileupShift(Double_t is) { _pileupShift=is; };
  inline Double_t GetpileupShift() const { return _pileupShift; };
  inline Double_t pileupShift() const { return GetpileupShift(); };

 // Bkgd Tree name (if set, is also a sign to set DoBkgdTree true
  inline void SetBkgdTreeName(const TString& tna) { _BkgdTreeName=tna; };
  inline TString GetBkgdTreeName() const { return _BkgdTreeName; };
  inline TString BkgdTreeName() const { return GetBkgdTreeName(); };
  //
  // b-tag algorithm
  inline void SetBTagAlgo(const TString& algo) { _bTagAlgo=algo; };
  inline void SetBTagAlgorithm(const TString& algo) { SetBTagAlgo(algo); };
  inline TString GetBTagAlgo() const { return _bTagAlgo; };
  inline TString GetBTagAlgorithm() const { return GetBTagAlgo(); };
  inline TString BTagAlgo() const { return GetBTagAlgo(); };

 // b-tag weight cut
  inline void SetBTagCut(Double_t cut) { _bTagCut=cut; };
  inline void SetBTagWeightCut(Double_t cut) { SetBTagCut(cut); };
  inline Double_t GetBTagCut() const { return _bTagCut; };
  inline Double_t GetBTagWeightCut() const { return GetBTagCut(); };
  inline Double_t BTagCut() const { return GetBTagCut(); };

 // mis-tag weight cut
  inline void SetMisTagCut(Double_t cutm) { _misTagCut=cutm; };
  inline void SetMisTagWeightCut(Double_t cutm) { SetMisTagCut(cutm); };
  inline Double_t GetMisTagCut() const { return _misTagCut; };
  inline Double_t GetMisTagWeightCut() const { return GetMisTagCut(); };
  inline Double_t MisTagCut() const { return GetMisTagCut(); };


  // Get the event weights for filling histograms
  inline Double_t GetEventWeight()       const { return _outputEventWeight; };
  inline Double_t GetOutputEventWeight() const { return _outputEventWeight; };
  inline Double_t GetTreeEventWeight()   const { return _treeEventWeight;   };
  inline Double_t GetGlobalEventWeight() const { return _globalEventWeight; };

  //pileup weight per event
  inline void SetEventPileupWeight(const Double_t& pileupweight) {_EventPileupWeight = pileupweight;}; //tagging prob per event
  inline Double_t GetEventPileupWeight() const {return _EventPileupWeight;}; 
  inline Double_t EventPileupWeight() const {return GetEventPileupWeight();}; 

  //pileup weight per event with min bias x-sec pushed up
  inline void SetEventPileupMinBiasUpWeight(const Double_t& pileupweight) {_EventPileupMinBiasUpWeight = pileupweight;}; 
  inline Double_t GetEventPileupMinBiasUpWeight() const {return _EventPileupMinBiasUpWeight;}; 
  inline Double_t EventPileupMinBiasUpWeight() const {return GetEventPileupMinBiasUpWeight();}; 

  //pileup weight per event with min bias x-sec pushed down
  inline void SetEventPileupMinBiasDownWeight(const Double_t& pileupweight) {_EventPileupMinBiasDownWeight = pileupweight;}; 
  inline Double_t GetEventPileupMinBiasDownWeight() const {return _EventPileupMinBiasDownWeight;}; 
  inline Double_t EventPileupMinBiasDownWeight() const {return GetEventPileupMinBiasDownWeight();}; 

 //b weight per event
  inline void SetEventbWeight(const Double_t& bweight) {_EventbWeight = bweight;}; //tagging prob per event
  inline Double_t GetEventbWeight() const {return _EventbWeight;}; 
  inline Double_t EventbWeight() const {return GetEventbWeight();}; 
  
 //W_tagging weight per event
  inline void SetEventWSF(const Double_t& w_WJet_) {_EventWSF = w_WJet_;}; //tagging prob per event
  inline Double_t GetEventWSF() const {return _EventWSF;};
  inline Double_t EventWSF() const {return GetEventWSF();};
// W_tagging weight up per event 
    inline void SetEventWSFUp(const Double_t& w_WJetUp_) {_EventWSFUp = w_WJetUp_;}; //tagging prob per event
    inline Double_t GetEventWSFUp() const {return _EventWSFUp;};
    inline Double_t EventWSFUp() const {return GetEventWSFUp();};
// W_tagging weight up per event 
     inline void SetEventWSFDown(const Double_t& w_WJetDown_) {_EventWSFDown = w_WJetDown_;}; //tagging prob per event
         inline Double_t GetEventWSFDown() const {return _EventWSFDown;};
         inline Double_t EventWSFDown() const {return GetEventWSFDown();};
             
  //lepton SF weight per event
  inline void SetEventLepSFWeight(const Double_t& lepSFweight) {_EventLepSFWeight = lepSFweight;};
  inline Double_t GetEventLepSFWeight() const {return _EventLepSFWeight;};
  inline Double_t EventLepSFWeight() const {return GetEventLepSFWeight();};

  //lepton SF weight up variation per event
  inline void SetEventLepSFWeightUp(const Double_t& lepSFweight) {_EventLepSFWeightUp = lepSFweight;};
  inline Double_t GetEventLepSFWeightUp() const {return _EventLepSFWeightUp;};
  inline Double_t EventLepSFWeightUp() const {return GetEventLepSFWeightUp();};

  //lepton SF weight down variation per event
  inline void SetEventLepSFWeightDown(const Double_t& lepSFweight) {_EventLepSFWeightDown = lepSFweight;};
  inline Double_t GetEventLepSFWeightDown() const {return _EventLepSFWeightDown;};
  inline Double_t EventLepSFWeightDown() const {return GetEventLepSFWeightDown();};

  //lepton SF weight per event
  inline void SetEventTrigSFWeight(const Double_t& lepSFweight) {_EventTrigSFWeight = lepSFweight;};
  inline Double_t GetEventTrigSFWeight() const {return _EventTrigSFWeight;};
  inline Double_t EventTrigSFWeight() const {return GetEventTrigSFWeight();};

  //lepton SF weight up variation per event
  inline void SetEventTrigSFWeightUp(const Double_t& lepSFweight) {_EventTrigSFWeightUp = lepSFweight;};
  inline Double_t GetEventTrigSFWeightUp() const {return _EventTrigSFWeightUp;};
  inline Double_t EventTrigSFWeightUp() const {return GetEventTrigSFWeightUp();};

  //lepton SF weight down variation per event
  inline void SetEventTrigSFWeightDown(const Double_t& lepSFweight) {_EventTrigSFWeightDown = lepSFweight;};
  inline Double_t GetEventTrigSFWeightDown() const {return _EventTrigSFWeightDown;};
  inline Double_t EventTrigSFWeightDown() const {return GetEventTrigSFWeightDown();};


  //set and get btag weights for central and systematics
  inline void SetEventbTagReshape(const Double_t& lepSFweight, std::string systName = "central" ) {_EventbTagReshape[systName] = lepSFweight;};
  inline Double_t GetEventbTagReshape(std::string systName = "central") {return _EventbTagReshape[systName];};
  inline Double_t EventbTagReshape() {return GetEventbTagReshape();};

  //set and get btag weights for central and systematics
  inline void SetEventMisTagReshape(const Double_t& lepSFweight, std::string systName = "central" ) {_EventMisTagReshape[systName] = lepSFweight;};
  inline Double_t GetEventMisTagReshape(std::string systName = "central") {return _EventMisTagReshape[systName];};
  inline Double_t EventMisTagReshape() {return GetEventMisTagReshape();};

  //Gen weight for event (actually just up or down)
  inline void SetGenWeight(const Double_t& genWeight) {_EventGenWeight = genWeight;};
  inline Double_t GetGenWeight() const {return _EventGenWeight;};
  inline Double_t GenWeight() const {return GetGenWeight();};

 //Tagging probabiliy weight per event
  inline void SetEventTagWeight(const Double_t& taggingweight) {_EventTagWeight = taggingweight;}; //tagging prob per event
  inline Double_t GetEventTagWeight() const {return _EventTagWeight;}; 
  inline Double_t EventTagWeight() const {return GetEventTagWeight();}; 

  //EventTagging probabiliy weight per event with a shift up on b SF
  inline void SetEventTagWeight_Bup(const Double_t& taggingweight_Bup) {_EventTagWeight_Bup = taggingweight_Bup;}; //tagging prob per event
  inline Double_t GetEventTagWeight_Bup() const {return _EventTagWeight_Bup;}; 
  inline Double_t EventTagWeight_Bup() const {return GetEventTagWeight_Bup();}; 

  //EventTagging probabiliy weight per event with a shift up on b SF
  inline void SetEventTagWeight_Bdown(const Double_t& taggingweight_Bdown) {_EventTagWeight_Bdown = taggingweight_Bdown;}; //tagging prob per event
  inline Double_t GetEventTagWeight_Bdown() const {return _EventTagWeight_Bdown;}; 
  inline Double_t EventTagWeight_Bdown() const {return GetEventTagWeight_Bdown();}; 

  //EventTagging probabiliy weight per event with a shift up on b SF
  inline void SetEventTagWeight_Lqup(const Double_t& taggingweight_Lqup) {_EventTagWeight_Lqup = taggingweight_Lqup;}; //tagging prob per event
  inline Double_t GetEventTagWeight_Lqup() const {return _EventTagWeight_Lqup;}; 
  inline Double_t EventTagWeight_Lqup() const {return GetEventTagWeight_Lqup();}; 

  //EventTagging probabiliy weight per event with a shift up on b SF
  inline void SetEventTagWeight_Lqdown(const Double_t& taggingweight_Lqdown) {_EventTagWeight_Lqdown = taggingweight_Lqdown;}; //tagging prob per event
  inline Double_t GetEventTagWeight_Lqdown() const {return _EventTagWeight_Lqdown;}; 
  inline Double_t EventTagWeight_Lqdown() const {return GetEventTagWeight_Lqdown();}; 

    // Get and set heavy flavor tool flag
  inline Int_t GetHFOR() const { return _HFOR; };
  inline void SetHFOR(Int_t i) { _HFOR=i; };

   // default to 1, set to 0 if event has be evaluated at end of EventContainer
  inline Bool_t GetIsFirstEvent() const { return _isFirstEvent; };
  inline void SetIsFirstEvent(Bool_t i) { _isFirstEvent=i; };

  inline Bool_t GetBadJetEvent() const { return _badJetEvent; };
  inline void SetBadJetEvent(Bool_t bb) { _badJetEvent=bb; };

  inline void SetnPvt(const int& nPvt) { nPvtx = nPvt; };
  inline Int_t GetnPvt() const { return nPvtx; };
  inline Int_t nPvt() const {return nPvtx;};
  
  inline TString GetChannelName() const {return _channelName;};

  inline void SetIsSimulation(bool isSim){isSimulation = isSim;};
  inline bool GetIsSimulation(){return isSimulation;};
  //inline bool isSimulation() { return _isSimulation;};


   // Get the configuration 
  TEnv * GetConfig() {return &_config; };

  ////////////////////////////////////
  // Object Helpers:
  ////////////////////////////////////

  TString _BkgdTreeName;
  Bool_t _doOnlyMVATree;
  Int_t _HFOR;
  Bool_t _isFirstEvent;

  Float_t _EventTagWeight; //tagging weight per event
  Float_t _EventTagWeight_Bup;//tagging weight per event applying the b SF up shift
  Float_t _EventTagWeight_Bdown;//tagging weight per event applying the b SF down shift
  Float_t _EventTagWeight_Lqup;//tagging weight per event applying the Lq SF up shift
  Float_t _EventTagWeight_Lqdown;//tagging weight per event applying the Lq SF down shift
  Float_t _EventPileupWeight;
  Float_t _EventPileupMinBiasUpWeight;
  Float_t _EventPileupMinBiasDownWeight;
  Float_t _EventbWeight; // this may be the same as the tagging weight, but I'm making it different anyway
  Float_t _EventWSF; // this may be the same as the tagging weight, but I'm making it different anyway
  Float_t _EventLepSFWeight;
  Float_t _EventTrigSFWeight;
  Float_t _EventGenWeight;
  std::map<std::string,Float_t> _EventbTagReshape;
  std::map<std::string,Float_t> _EventMisTagReshape;

  //Add in the systematic variations to the SFs
  Float_t _EventLepSFWeightUp;
  Float_t _EventLepSFWeightDown;
  Float_t _EventTrigSFWeightUp;
  Float_t _EventTrigSFWeightDown;
  Float_t _EventWSFUp;
  Float_t _EventWSFDown;
  //MultijetJESUncertaintyProvider myJES;
  // CalibrationDataVariables CalibVar;
  std::string taggername;
  
  Top_MET top_met;
  MET_Cleaning_Utils MET_Util;
  //CalibrationDataInterfaceROOT CalibROOT("SV0", "config/btag/BTagCalibration.env");
  //CalibrationDataInterfaceROOT CalibROOT;

  // Information for partons
  Int_t runNumber;              // run number from the root tree
  Int_t eventNumber;            // event number from the root tree
  Bool_t larError;
  Float_t actualIntPerXing;
  Float_t averageIntPerXing;
  Float_t distns;
  Float_t distbunch;
  Float_t bcid;
  Int_t Pvtxall_n;

  Bool_t _badJetEvent;

  Bool_t isSimulation;
  Int_t safejetevent;
  Int_t safejeteventup;
  Int_t safejeteventdown;

  //top MET Vectors
  std::vector<TLorentzVector>  Electrons_tlv;
  std::vector<TLorentzVector>  Electrons_tlv_scaled;
  std::vector<TLorentzVector>  Jets_tlv_scaled;
  std::vector<TLorentzVector>  Jets_tlv;
  std::vector<TLorentzVector>  Jets_JESup_tlv;
  std::vector<TLorentzVector>  Jets_emscale_tlv;
  std::vector<TLorentzVector>  Muons_tlv_scaled;
  std::vector<TLorentzVector>  Muons_tlv;
  std::vector<TLorentzVector>  Muons_ms_tlv_scaled;
  std::vector<TLorentzVector>  Muons_ms_tlv;
  std::vector<TLorentzVector>  Muons_track_tlv_scaled;
  std::vector<TLorentzVector>  Muons_track_tlv;

  // Vectors of Reconstructed particles
  std::vector<Electron>   electrons;
  std::vector<Electron>   tightElectrons;
  std::vector<Electron>   vetoElectrons;
  std::vector<Electron>   ptetaElectrons;
  std::vector<Electron>   isolatedElectrons;
  std::vector<Electron>   unIsolatedElectrons;
  std::vector<Muon>       muons;
  std::vector<Muon>       tightMuons;
  std::vector<Muon>       vetoMuons;
  std::vector<Muon>       ptetaMuons;
  std::vector<Muon>       isolatedMuons;
  std::vector<Muon>       unIsolatedMuons;
  std::vector<Tau>        taus;
  std::vector<Jet>        jets;
  std::vector<BoostedJet> boostedjets;
  std::vector<Jet>        alljets;//no jet checks
  std::vector<Jet>        jetors;//overlap removal vector for electron/jets
  std::vector<Jet>        jetjetors;//overlap removal vector for jet/jets
  std::vector<Jet>        jetms;//overlap removal vector for muons
  std::vector<Jet>        taggedJets;
  std::vector<Jet>        unTaggedJets;
  std::vector<Jet>        bLabeledJets;
  std::vector<Jet>        cLabeledJets;
  std::vector<Jet>        tauLabeledJets;
  std::vector<Jet>        lightQuarkLabeledJets;
  std::vector<Neutrino>   neutrinos;

  //Pointers to collections
  std::vector<Electron> * electronsToUsePtr;
  std::vector<Muon>     * muonsToUsePtr;
  std::vector<BoostedJet> * boostedjetsToUsePtr;
  std::vector<Jet> * jetsToUsePtr;
  std::vector<Jet> * taggedjetsToUsePtr;

  std::vector<TLorentzVector>      jetmsSpecial;//MET tool
  
  // also store MET
  Double_t missingEt;
  Double_t missingEx;
  Double_t missingPhi;
  Double_t missingEy;
  TLorentzVector missingEtVec;
  // and the sum of scalar ET
  Double_t sumEt;

  Double_t missingEt_xy;
  Double_t missingEx_xy;
  Double_t missingEy_xy;
  Double_t missingPhi_xy;
  TLorentzVector missingEtVec_xy;

  //Track whterh it passes the MET filters
  Int_t passesMETFilters;

  //Various variables for plotting primary vertex information
  Int_t nPvtx;
  Int_t trueInteractions;
  Int_t npuVertices;

  Double_t particleHT;
  Double_t particleMTotal;

  int closeindex;
  double ejordr;
  double bestjetdr;
  Bool_t jeteoverlap;

  // MC particles:
  std::vector<MCParticle> MCParticles;    // The full list of all MC particles
  std::vector<MCMuon> MCMuons;            // MC Muons
  std::vector<MCElectron> MCElectrons;    // MC Electrons
  std::vector<MCTau> MCTaus;              // MC Taus
  std::vector<MCJet> MCJets;              // MC Jets (UDSCB)
  std::vector<MCJet> MCBJets;             // MC B Jets
  std::vector<MCJet> MCCJets;             // MC C Jets
  std::vector<MCJet> MCLightQuarkJets;    // U, D, S Jets
  std::vector<MCTop> MCTops;              // MC Tops
  std::vector<MCW> MCWs;                  // MC Ws
  std::vector<MCNeutrino> MCNeutrinos;    // MC Neutrinos
  std::vector<MCBQuark>  MCBQuarks;       // MC B Quarks
  std::vector<MCCQuark>  MCCQuarks;       // MC C Quarks

  // Iterators
  std::vector<MCParticle>::iterator MCParticleIterator;
  std::vector<MCMuon>::iterator MCMuonIterator;
  std::vector<MCMuon>::iterator MCMuonNewEnd;
  std::vector<MCElectron>::iterator MCElectronIterator;
  std::vector<MCElectron>::iterator MCElectronNewEnd;
  std::vector<MCTau>::iterator MCTauIterator;
  std::vector<MCTau>::iterator MCTauNewEnd;
  std::vector<MCJet>::iterator MCJetIterator;
  std::vector<MCJet>::iterator MCJetNewEnd;
  std::vector<MCJet>::iterator MCBJetIterator;
  std::vector<MCJet>::iterator MCBJetNewEnd;
  std::vector<MCJet>::iterator MCCJetIterator;
  std::vector<MCJet>::iterator MCCJetNewEnd;
  std::vector<MCJet>::iterator MCLightQuarkJetIterator;
  std::vector<MCJet>::iterator MCLightQuarkJetNewEnd;
  std::vector<MCTop>::iterator MCTopIterator;
  std::vector<MCTop>::iterator MCTopNewEnd;
  std::vector<MCW>::iterator MCWIterator;
  std::vector<MCW>::iterator MCWNewEnd;
  std::vector<MCNeutrino>::iterator MCNeutrinoIterator;
  std::vector<MCNeutrino>::iterator MCNeutrinoNewEnd;
  std::vector<MCBQuark>::iterator MCBQuarkIterator;
  std::vector<MCBQuark>::iterator MCBQuarkNewEnd;

  std::string jesInputFile; 
  // Final State Variables
  int finalLeptonPDGID;            //Leptons pdgId for HistogrammingTV_float
  std::vector<Jet> finalJetVector;  // Jets
  Particle finalLepton;             // Leptons
  Neutrino finalNeutrino;           // Neutrino
  Particle finalW;                  // W
  Particle finalCM;                 // Final CM systrem (Jets + Leptons + Neutrino)

  // Variables for making top quarks
  Jet finalTop_leadingJet;
  Jet finalTop_bestJet;
  Jet finalTop_bestJet12;
  Jet finalTop_bgJet;
  Jet bgJet;
  Jet finalTop_taggedJet;
  Jet finalNeutrinoTopC;
  Jet finalWTopC;
  Jet finalTop_taggedJetTopC;

  // Variables for making top quarks
  Jet *bestJet12Ptr;   
  Jet *nonBestJet12Ptr;
  Jet *bestJetPtr;     
  Jet *bgJetbJetPtr;   
  Jet *nonBgJetPtr;    
  Jet *bgJetGluonPtr;  
  Jet *taggedJetPtr;   
  Jet *untaggedPtJetPtr; 
  Jet *untaggedEtaJetPtr;
  Jet *forwardJetPtr;  
  Jet *leadingJetPtr;

  Int_t iBestJet; 
  int useElecEnergyScale;
  int useElecEnergySmear;
  bool useMCWithConstantTerm;
  //bool isMCFlag;

protected:
  
  // Debugging level
  Int_t _debugLevel;

  // Use TEnv to read configuration parameters
  TEnv _config;
  TEnv _JESconfig;
  bool _JESconfigread;
private:

  //The name of the channel being used. this is currently used for synch, but might be used for trigger or something.
  TString _channelName;
  
  // Target top mass used in the best jet algorithm
  Double_t _targetTopMass;

  // jet energy scale relative uncertainty
  Double_t _jesError;
  // jet energy scale shift, either +1 (for +1 sigma, -1 (for -1 sigma) or 0 (no shift)
  Int_t _jesShift;
  Int_t _jesUShift;
  Int_t _jesPtShift;
  Int_t _jesLowPtShift;
  Int_t _jesHighPtShift;
  Int_t _jesEtaShift;

  Double_t _celloutShift;
  Double_t _softjetShift;
  Double_t _pileupShift;
  Int_t _larShift;

  //met uncertainty shift, +1,-1 or 0 as above.
  Int_t _metShift;

  // which b-tag algorithm to use?
  TString _bTagAlgo;
 // which b-tag weight cut to use?
  Double_t _bTagCut;
// which mis-tag weight cut to use?
  Double_t _misTagCut;

  // The classes that contain the information from the trees:
  //CollectionTree *_collectionTree;
  EventTree      *_eventTree;
  FastSimTree    *_fastSimTree;
  TruthTree      *_truthTree;
  // Counters to keep track of location in event chain
  Int_t _eventCount;

  Bool_t _doFastSim;    // should we run over fast simulation or full simulation?
  Bool_t _doTruth;      // should we include the truth tree at all?
  Bool_t _doSkim;       // should we skim?
  Bool_t _doConfig;     // should we include the Bkgd tree at all?
  Bool_t _useUnisolatedLeptons;   // should we use unisolated leptons? (for QCD estimation)
  Int_t _trigID;

  // Info for source that we are reading (i.e. tb, tq, tW, tt_lepjets, tt_dilep, etc.
  TString _sourceName;    // Source Name
  Int_t   _sourceNumber;  // Source Number

  // Event Weight
  Double_t _globalEventWeight;         // Global event weight, unchanging event to event
  Double_t _treeEventWeight;           // Event weight read in from the tree
  Double_t _outputEventWeight;         // Output event weight passed to the histograms
 
  // Event objects. By putting these here we can define them once along with their cuts instead of creating them 
  // once for each event and accessing the config for each particle.
  Muon newMuon;
  Electron newElectron;
  Jet newJet;
  BoostedJet newBoostedJet;

  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(EventContainer,0)
  ////////////////////////////////////////////////////////////////////////////////
};


// MCParticle - Sort by PT - decending order
inline Bool_t operator <(const Particle& left, const Particle& right)
{
  if (left.Pt() <= right.Pt()) return kFALSE;
  else return kTRUE;
} //<


inline Bool_t RemoveDuplicate(const MCParticle& one, const MCParticle& two)
{
  return (TMath::Abs(one.Pt() - two.Pt()) <= 0.001);
} //RemoveEntry

inline Bool_t RemoveEntry(const Particle& particle) {return (particle.Pt() == 0);} 

#endif
