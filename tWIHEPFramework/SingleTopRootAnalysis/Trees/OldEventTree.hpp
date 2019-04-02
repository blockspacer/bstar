//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 16 21:43:44 2015 by ROOT version 6.02/05
// from TTree BOOM/BOOM
// found on file: OutTree.root
//////////////////////////////////////////////////////////

#ifndef EventTree_h
#define EventTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
using namespace std;
#include <vector>

class EventTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<int>     *Trigger_decision;
   vector<string>  *Trigger_names;
   vector<double>  *Muon_pt;
   vector<double>  *Muon_eta;
   vector<double>  *Muon_phi;
   vector<double>  *Muon_p;
   vector<double>  *Muon_energy;
   vector<double>  *Muon_charge;
   vector<bool>    *Muon_tight;
   vector<bool>    *Muon_loose;
   vector<bool>    *Muon_soft;
   vector<bool>    *Muon_isHighPt;
   vector<bool>    *Muon_pf;
   vector<double>  *Muon_isoCharged;
   vector<double>  *Muon_isoSum;
   vector<double>  *Muon_isoCharParPt;
   vector<double>  *Muon_chi2;
   vector<double>  *Muon_validHits;
   vector<double>  *Muon_validHitsInner;
   vector<double>  *Muon_matchedStat;
   vector<double>  *Muon_dxy;
   vector<double>  *Muon_TLayers;
   vector<double>  *Muon_dz;
   vector<double>  *Muon_isoNeutralHadron;
   vector<double>  *Muon_isoPhoton;
   vector<double>  *Muon_isoPU;
   vector<double>  *patElectron_pt;
   vector<double>  *patElectron_eta;
   vector<double>  *patElectron_phi;
   vector<double>  *patElectron_energy;
   vector<double>  *patElectron_charge;
   vector<int>     *patElectron_isPassVeto;
   vector<int>     *patElectron_isPassLoose;
   vector<int>     *patElectron_isPassMedium;
   vector<int>     *patElectron_isPassTight;
   vector<int>     *patElectron_isPassHEEPId;
   vector<double>  *patElectron_d0;
   vector<double>  *patElectron_dz;
   vector<int>     *patElectron_expectedMissingInnerHits;
   vector<int>     *patElectron_passConversionVeto;
   vector<double>  *patElectron_isoChargedHadrons;
   vector<double>  *patElectron_isoNeutralHadrons;
   vector<double>  *patElectron_isoPhotons;
   vector<double>  *patElectron_isoPU;
   vector<double>  *Tau_eta;
   vector<double>  *Tau_phi;
   vector<double>  *Tau_pt;
   vector<double>  *Tau_energy;
   vector<double>  *Tau_charge;
   vector<int>     *Tau_decayModeFinding;
   vector<int>     *Tau_decayModeFindingNewDMs;
   vector<double>  *Tau_chargedIsoPtSum;
   vector<double>  *Tau_neutralIsoPtSum;
   vector<int>     *Tau_againstMuonTight3;
   vector<int>     *Tau_againstElectronMVATightMVA5;
   vector<double>  *Tau_nProngs;
   vector<double>  *Tau_puCorrPtSum;
   vector<int>     *Tau_byLooseCombinedIsolationDeltaBetaCorr;
   vector<int>     *Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *Tau_byMediumCombinedIsolationDeltaBetaCorr;
   vector<int>     *Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *Tau_byTightCombinedIsolationDeltaBetaCorr;
   vector<int>     *Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;
   vector<int>     *Tau_byLooseIsolationMVA3newDMwLT;
   vector<int>     *Tau_byLooseIsolationMVA3newDMwoLT;
   vector<int>     *Tau_byLooseIsolationMva3oldDMwLT;
   vector<int>     *Tau_byLooseIsolationMVA3oldDMwoLT;
   vector<int>     *Tau_byMediumIsolationMVA3newDMwLT;
   vector<int>     *Tau_byMediumIsolationMVA3newDMwoLT;
   vector<int>     *Tau_byMediumIsolationMva3oldDMwLT;
   vector<int>     *Tau_byMediumIsolationMVA3oldDMwoLT;
   vector<int>     *Tau_byTightIsolationMVA3newDMwLT;
   vector<int>     *Tau_byTightIsolationMVA3newDMwoLT;
   vector<int>     *Tau_byTightIsolationMva3oldDMwLT;
   vector<int>     *Tau_byTightIsolationMVA3oldDMwoLT;
   vector<int>     *Tau_againstMuonLoose2;
   vector<int>     *Tau_againstMuonLoose3;
   vector<int>     *Tau_againstMuonTight2;
   vector<int>     *Tau_againstElectronMVALooseMVA5;
   vector<int>     *Tau_againstElectronMVAMediumMVA5;
   vector<int>     *Tau_byVLooseCombinedIsolationDeltaBetaCorr;
   vector<double>  *Tau_leadChargedCandPt;
   vector<double>  *Tau_leadChargedCandCharge;
   vector<double>  *Tau_leadChargedCandEta;
   vector<double>  *Tau_leadChargedCandPhi;
   vector<double>  *Jet_pt;
   vector<double>  *Jet_eta;
   vector<double>  *Jet_phi;
   vector<double>  *Jet_energy;
   vector<double>  *Jet_bDiscriminator;
   vector<double>  *Jet_pileupId;
   vector<double>  *Jet_mass;
   vector<double>  *Jet_neutralHadEnergyFraction;
   vector<double>  *Jet_neutralEmEmEnergyFraction;
   vector<double>  *Jet_chargedHadronEnergyFraction;
   vector<double>  *Jet_chargedEmEnergyFraction;
   vector<double>  *Jet_muonEnergyFraction;
   vector<int>     *Jet_numberOfConstituents;
   vector<int>     *Jet_chargedMultiplicity;
   vector<double>  *Jet_electronEnergy;
   vector<double>  *Jet_photonEnergy;
   vector<double>  *UncorrJet_pt;
   vector<float>   *Photon_pt;
   vector<float>   *Photon_eta;
   vector<float>   *Photon_phi;
   vector<float>   *Photon_energy;
   vector<float>   *Photon_et;
   vector<float>   *Photon_HoverE;
   vector<float>   *Photon_phoR9;
   vector<float>   *Photon_SigmaIEtaIEta;
   vector<float>   *Photon_SigmaIPhiIPhi;
   vector<float>   *Photon_PFChIso;
   vector<float>   *Photon_PFPhoIso;
   vector<float>   *Photon_PFNeuIso;
   vector<int>     *Photon_EleVeto;
   vector<int>     *Photon_hasPixelSeed;

   // List of branches
   TBranch        *b_Trigger_decision;   //!
   TBranch        *b_Trigger_names;   //!
   TBranch        *b_Muon_pt;   //!
   TBranch        *b_Muon_eta;   //!
   TBranch        *b_Muon_phi;   //!
   TBranch        *b_Muon_p;   //!
   TBranch        *b_Muon_energy;   //!
   TBranch        *b_Muon_charge;   //!
   TBranch        *b_Muon_tight;   //!
   TBranch        *b_Muon_loose;   //!
   TBranch        *b_Muon_soft;   //!
   TBranch        *b_Muon_isHighPt;   //!
   TBranch        *b_Muon_pf;   //!
   TBranch        *b_Muon_isoCharged;   //!
   TBranch        *b_Muon_isoSum;   //!
   TBranch        *b_Muon_isoCharParPt;   //!
   TBranch        *b_Muon_chi2;   //!
   TBranch        *b_Muon_validHits;   //!
   TBranch        *b_Muon_validHitsInner;   //!
   TBranch        *b_Muon_matchedStat;   //!
   TBranch        *b_Muon_dxy;   //!
   TBranch        *b_Muon_TLayers;   //!
   TBranch        *b_Muon_dz;   //!
   TBranch        *b_Muon_isoNeutralHadron;   //!
   TBranch        *b_Muon_isoPhoton;   //!
   TBranch        *b_Muon_isoPU;   //!
   TBranch        *b_patElectron_pt;   //!
   TBranch        *b_patElectron_eta;   //!
   TBranch        *b_patElectron_phi;   //!
   TBranch        *b_patElectron_energy;   //!
   TBranch        *b_patElectron_charge;   //!
   TBranch        *b_patElectron_isPassVeto;   //!
   TBranch        *b_patElectron_isPassLoose;   //!
   TBranch        *b_patElectron_isPassMedium;   //!
   TBranch        *b_patElectron_isPassTight;   //!
   TBranch        *b_patElectron_isPassHEEPId;   //!
   TBranch        *b_patElectron_d0;   //!
   TBranch        *b_patElectron_dz;   //!
   TBranch        *b_patElectron_expectedMissingInnerHits;   //!
   TBranch        *b_patElectron_passConversionVeto;   //!
   TBranch        *b_patElectron_isoChargedHadrons;   //!
   TBranch        *b_patElectron_isoNeutralHadrons;   //!
   TBranch        *b_patElectron_isoPhotons;   //!
   TBranch        *b_patElectron_isoPU;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_pt;   //!
   TBranch        *b_Tau_energy;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_decayModeFinding;   //!
   TBranch        *b_Tau_decayModeFindingNewDMs;   //!
   TBranch        *b_Tau_chargedIsoPtSum;   //!
   TBranch        *b_Tau_neutralIsoPtSum;   //!
   TBranch        *b_Tau_againstMuonTight3;   //!
   TBranch        *b_Tau_againstElectronMVATightMVA5;   //!
   TBranch        *b_Tau_nProngs;   //!
   TBranch        *b_Tau_puCorrPtSum;   //!
   TBranch        *b_Tau_byLooseCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byMediumCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byTightCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_Tau_byLooseIsolationMVA3newDMwLT;   //!
   TBranch        *b_Tau_byLooseIsolationMVA3newDMwoLT;   //!
   TBranch        *b_Tau_byLooseIsolationMva3oldDMwLT;   //!
   TBranch        *b_Tau_byLooseIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVA3newDMwLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVA3newDMwoLT;   //!
   TBranch        *b_Tau_byMediumIsolationMva3oldDMwLT;   //!
   TBranch        *b_Tau_byMediumIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_byTightIsolationMVA3newDMwLT;   //!
   TBranch        *b_Tau_byTightIsolationMVA3newDMwoLT;   //!
   TBranch        *b_Tau_byTightIsolationMva3oldDMwLT;   //!
   TBranch        *b_Tau_byTightIsolationMVA3oldDMwoLT;   //!
   TBranch        *b_Tau_againstMuonLoose2;   //!
   TBranch        *b_Tau_againstMuonLoose3;   //!
   TBranch        *b_Tau_againstMuonTight2;   //!
   TBranch        *b_Tau_againstElectronMVALooseMVA5;   //!
   TBranch        *b_Tau_againstElectronMVAMediumMVA5;   //!
   TBranch        *b_Tau_byVLooseCombinedIsolationDeltaBetaCorr;   //!
   TBranch        *b_Tau_leadChargedCandPt;   //!
   TBranch        *b_Tau_leadChargedCandCharge;   //!
   TBranch        *b_Tau_leadChargedCandEta;   //!
   TBranch        *b_Tau_leadChargedCandPhi;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_energy;   //!
   TBranch        *b_Jet_bDiscriminator;   //!
   TBranch        *b_Jet_pileupId;   //!
   TBranch        *b_Jet_mass;   //!
   TBranch        *b_Jet_neutralHadEnergyFraction;   //!
   TBranch        *b_Jet_neutralEmEmEnergyFraction;   //!
   TBranch        *b_Jet_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jet_chargedEmEnergyFraction;   //!
   TBranch        *b_Jet_muonEnergyFraction;   //!
   TBranch        *b_Jet_numberOfConstituents;   //!
   TBranch        *b_Jet_chargedMultiplicity;   //!
   TBranch        *b_Jet_electronEnergy;   //!
   TBranch        *b_Jet_photonEnergy;   //!
   TBranch        *b_UncorrJet_pt;   //!
   TBranch        *b_Photon_pt;   //!
   TBranch        *b_Photon_eta;   //!
   TBranch        *b_Photon_phi;   //!
   TBranch        *b_Photon_energy;   //!
   TBranch        *b_Photon_et;   //!
   TBranch        *b_Photon_HoverE;   //!
   TBranch        *b_Photon_phoR9;   //!
   TBranch        *b_Photon_SigmaIEtaIEta;   //!
   TBranch        *b_Photon_SigmaIPhiIPhi;   //!
   TBranch        *b_Photon_PFChIso;   //!
   TBranch        *b_Photon_PFPhoIso;   //!
   TBranch        *b_Photon_PFNeuIso;   //!
   TBranch        *b_Photon_EleVeto;   //!
   TBranch        *b_Photon_hasPixelSeed;   //!

   EventTree(TTree *tree=0);
   virtual ~EventTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EventTree_cxx
EventTree::EventTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("OutTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("OutTree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("OutTree.root:/TNT");
      dir->GetObject("BOOM",tree);

   }
   Init(tree);
}

EventTree::~EventTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Trigger_decision = 0;
   Trigger_names = 0;
   Muon_pt = 0;
   Muon_eta = 0;
   Muon_phi = 0;
   Muon_p = 0;
   Muon_energy = 0;
   Muon_charge = 0;
   Muon_tight = 0;
   Muon_loose = 0;
   Muon_soft = 0;
   Muon_isHighPt = 0;
   Muon_pf = 0;
   Muon_isoCharged = 0;
   Muon_isoSum = 0;
   Muon_isoCharParPt = 0;
   Muon_chi2 = 0;
   Muon_validHits = 0;
   Muon_validHitsInner = 0;
   Muon_matchedStat = 0;
   Muon_dxy = 0;
   Muon_TLayers = 0;
   Muon_dz = 0;
   Muon_isoNeutralHadron = 0;
   Muon_isoPhoton = 0;
   Muon_isoPU = 0;
   patElectron_pt = 0;
   patElectron_eta = 0;
   patElectron_phi = 0;
   patElectron_energy = 0;
   patElectron_charge = 0;
   patElectron_isPassVeto = 0;
   patElectron_isPassLoose = 0;
   patElectron_isPassMedium = 0;
   patElectron_isPassTight = 0;
   patElectron_isPassHEEPId = 0;
   patElectron_d0 = 0;
   patElectron_dz = 0;
   patElectron_expectedMissingInnerHits = 0;
   patElectron_passConversionVeto = 0;
   patElectron_isoChargedHadrons = 0;
   patElectron_isoNeutralHadrons = 0;
   patElectron_isoPhotons = 0;
   patElectron_isoPU = 0;
   Tau_eta = 0;
   Tau_phi = 0;
   Tau_pt = 0;
   Tau_energy = 0;
   Tau_charge = 0;
   Tau_decayModeFinding = 0;
   Tau_decayModeFindingNewDMs = 0;
   Tau_chargedIsoPtSum = 0;
   Tau_neutralIsoPtSum = 0;
   Tau_againstMuonTight3 = 0;
   Tau_againstElectronMVATightMVA5 = 0;
   Tau_nProngs = 0;
   Tau_puCorrPtSum = 0;
   Tau_byLooseCombinedIsolationDeltaBetaCorr = 0;
   Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byMediumCombinedIsolationDeltaBetaCorr = 0;
   Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byTightCombinedIsolationDeltaBetaCorr = 0;
   Tau_byTightCombinedIsolationDeltaBetaCorr3Hits = 0;
   Tau_byLooseIsolationMVA3newDMwLT = 0;
   Tau_byLooseIsolationMVA3newDMwoLT = 0;
   Tau_byLooseIsolationMva3oldDMwLT = 0;
   Tau_byLooseIsolationMVA3oldDMwoLT = 0;
   Tau_byMediumIsolationMVA3newDMwLT = 0;
   Tau_byMediumIsolationMVA3newDMwoLT = 0;
   Tau_byMediumIsolationMva3oldDMwLT = 0;
   Tau_byMediumIsolationMVA3oldDMwoLT = 0;
   Tau_byTightIsolationMVA3newDMwLT = 0;
   Tau_byTightIsolationMVA3newDMwoLT = 0;
   Tau_byTightIsolationMva3oldDMwLT = 0;
   Tau_byTightIsolationMVA3oldDMwoLT = 0;
   Tau_againstMuonLoose2 = 0;
   Tau_againstMuonLoose3 = 0;
   Tau_againstMuonTight2 = 0;
   Tau_againstElectronMVALooseMVA5 = 0;
   Tau_againstElectronMVAMediumMVA5 = 0;
   Tau_byVLooseCombinedIsolationDeltaBetaCorr = 0;
   Tau_leadChargedCandPt = 0;
   Tau_leadChargedCandCharge = 0;
   Tau_leadChargedCandEta = 0;
   Tau_leadChargedCandPhi = 0;
   Jet_pt = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_energy = 0;
   Jet_bDiscriminator = 0;
   Jet_pileupId = 0;
   Jet_mass = 0;
   Jet_neutralHadEnergyFraction = 0;
   Jet_neutralEmEmEnergyFraction = 0;
   Jet_chargedHadronEnergyFraction = 0;
   Jet_chargedEmEnergyFraction = 0;
   Jet_muonEnergyFraction = 0;
   Jet_numberOfConstituents = 0;
   Jet_chargedMultiplicity = 0;
   Jet_electronEnergy = 0;
   Jet_photonEnergy = 0;
   UncorrJet_pt = 0;
   Photon_pt = 0;
   Photon_eta = 0;
   Photon_phi = 0;
   Photon_energy = 0;
   Photon_et = 0;
   Photon_HoverE = 0;
   Photon_phoR9 = 0;
   Photon_SigmaIEtaIEta = 0;
   Photon_SigmaIPhiIPhi = 0;
   Photon_PFChIso = 0;
   Photon_PFPhoIso = 0;
   Photon_PFNeuIso = 0;
   Photon_EleVeto = 0;
   Photon_hasPixelSeed = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Trigger_decision", &Trigger_decision, &b_Trigger_decision);
   fChain->SetBranchAddress("Trigger_names", &Trigger_names, &b_Trigger_names);
   fChain->SetBranchAddress("Muon_pt", &Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_eta", &Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_phi", &Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_p", &Muon_p, &b_Muon_p);
   fChain->SetBranchAddress("Muon_energy", &Muon_energy, &b_Muon_energy);
   fChain->SetBranchAddress("Muon_charge", &Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_tight", &Muon_tight, &b_Muon_tight);
   fChain->SetBranchAddress("Muon_loose", &Muon_loose, &b_Muon_loose);
   fChain->SetBranchAddress("Muon_soft", &Muon_soft, &b_Muon_soft);
   fChain->SetBranchAddress("Muon_isHighPt", &Muon_isHighPt, &b_Muon_isHighPt);
   fChain->SetBranchAddress("Muon_pf", &Muon_pf, &b_Muon_pf);
   fChain->SetBranchAddress("Muon_isoCharged", &Muon_isoCharged, &b_Muon_isoCharged);
   fChain->SetBranchAddress("Muon_isoSum", &Muon_isoSum, &b_Muon_isoSum);
   fChain->SetBranchAddress("Muon_isoCharParPt", &Muon_isoCharParPt, &b_Muon_isoCharParPt);
   fChain->SetBranchAddress("Muon_chi2", &Muon_chi2, &b_Muon_chi2);
   fChain->SetBranchAddress("Muon_validHits", &Muon_validHits, &b_Muon_validHits);
   fChain->SetBranchAddress("Muon_validHitsInner", &Muon_validHitsInner, &b_Muon_validHitsInner);
   fChain->SetBranchAddress("Muon_matchedStat", &Muon_matchedStat, &b_Muon_matchedStat);
   fChain->SetBranchAddress("Muon_dxy", &Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_TLayers", &Muon_TLayers, &b_Muon_TLayers);
   fChain->SetBranchAddress("Muon_dz", &Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_isoNeutralHadron", &Muon_isoNeutralHadron, &b_Muon_isoNeutralHadron);
   fChain->SetBranchAddress("Muon_isoPhoton", &Muon_isoPhoton, &b_Muon_isoPhoton);
   fChain->SetBranchAddress("Muon_isoPU", &Muon_isoPU, &b_Muon_isoPU);
   fChain->SetBranchAddress("patElectron_pt", &patElectron_pt, &b_patElectron_pt);
   fChain->SetBranchAddress("patElectron_eta", &patElectron_eta, &b_patElectron_eta);
   fChain->SetBranchAddress("patElectron_phi", &patElectron_phi, &b_patElectron_phi);
   fChain->SetBranchAddress("patElectron_energy", &patElectron_energy, &b_patElectron_energy);
   fChain->SetBranchAddress("patElectron_charge", &patElectron_charge, &b_patElectron_charge);
   fChain->SetBranchAddress("patElectron_isPassVeto", &patElectron_isPassVeto, &b_patElectron_isPassVeto);
   fChain->SetBranchAddress("patElectron_isPassLoose", &patElectron_isPassLoose, &b_patElectron_isPassLoose);
   fChain->SetBranchAddress("patElectron_isPassMedium", &patElectron_isPassMedium, &b_patElectron_isPassMedium);
   fChain->SetBranchAddress("patElectron_isPassTight", &patElectron_isPassTight, &b_patElectron_isPassTight);
   fChain->SetBranchAddress("patElectron_isPassHEEPId", &patElectron_isPassHEEPId, &b_patElectron_isPassHEEPId);
   fChain->SetBranchAddress("patElectron_d0", &patElectron_d0, &b_patElectron_d0);
   fChain->SetBranchAddress("patElectron_dz", &patElectron_dz, &b_patElectron_dz);
   fChain->SetBranchAddress("patElectron_expectedMissingInnerHits", &patElectron_expectedMissingInnerHits, &b_patElectron_expectedMissingInnerHits);
   fChain->SetBranchAddress("patElectron_passConversionVeto", &patElectron_passConversionVeto, &b_patElectron_passConversionVeto);
   fChain->SetBranchAddress("patElectron_isoChargedHadrons", &patElectron_isoChargedHadrons, &b_patElectron_isoChargedHadrons);
   fChain->SetBranchAddress("patElectron_isoNeutralHadrons", &patElectron_isoNeutralHadrons, &b_patElectron_isoNeutralHadrons);
   fChain->SetBranchAddress("patElectron_isoPhotons", &patElectron_isoPhotons, &b_patElectron_isoPhotons);
   fChain->SetBranchAddress("patElectron_isoPU", &patElectron_isoPU, &b_patElectron_isoPU);
   fChain->SetBranchAddress("Tau_eta", &Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_phi", &Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_pt", &Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_energy", &Tau_energy, &b_Tau_energy);
   fChain->SetBranchAddress("Tau_charge", &Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_decayModeFinding", &Tau_decayModeFinding, &b_Tau_decayModeFinding);
   fChain->SetBranchAddress("Tau_decayModeFindingNewDMs", &Tau_decayModeFindingNewDMs, &b_Tau_decayModeFindingNewDMs);
   fChain->SetBranchAddress("Tau_chargedIsoPtSum", &Tau_chargedIsoPtSum, &b_Tau_chargedIsoPtSum);
   fChain->SetBranchAddress("Tau_neutralIsoPtSum", &Tau_neutralIsoPtSum, &b_Tau_neutralIsoPtSum);
   fChain->SetBranchAddress("Tau_againstMuonTight3", &Tau_againstMuonTight3, &b_Tau_againstMuonTight3);
   fChain->SetBranchAddress("Tau_againstElectronMVATightMVA5", &Tau_againstElectronMVATightMVA5, &b_Tau_againstElectronMVATightMVA5);
   fChain->SetBranchAddress("Tau_nProngs", &Tau_nProngs, &b_Tau_nProngs);
   fChain->SetBranchAddress("Tau_puCorrPtSum", &Tau_puCorrPtSum, &b_Tau_puCorrPtSum);
   fChain->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr", &Tau_byLooseCombinedIsolationDeltaBetaCorr, &b_Tau_byLooseCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits", &Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byLooseCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr", &Tau_byMediumCombinedIsolationDeltaBetaCorr, &b_Tau_byMediumCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits", &Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byMediumCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr", &Tau_byTightCombinedIsolationDeltaBetaCorr, &b_Tau_byTightCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("Tau_byTightCombinedIsolationDeltaBetaCorr3Hits", &Tau_byTightCombinedIsolationDeltaBetaCorr3Hits, &b_Tau_byTightCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("Tau_byLooseIsolationMVA3newDMwLT", &Tau_byLooseIsolationMVA3newDMwLT, &b_Tau_byLooseIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("Tau_byLooseIsolationMVA3newDMwoLT", &Tau_byLooseIsolationMVA3newDMwoLT, &b_Tau_byLooseIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("Tau_byLooseIsolationMva3oldDMwLT", &Tau_byLooseIsolationMva3oldDMwLT, &b_Tau_byLooseIsolationMva3oldDMwLT);
   fChain->SetBranchAddress("Tau_byLooseIsolationMVA3oldDMwoLT", &Tau_byLooseIsolationMVA3oldDMwoLT, &b_Tau_byLooseIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("Tau_byMediumIsolationMVA3newDMwLT", &Tau_byMediumIsolationMVA3newDMwLT, &b_Tau_byMediumIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("Tau_byMediumIsolationMVA3newDMwoLT", &Tau_byMediumIsolationMVA3newDMwoLT, &b_Tau_byMediumIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("Tau_byMediumIsolationMva3oldDMwLT", &Tau_byMediumIsolationMva3oldDMwLT, &b_Tau_byMediumIsolationMva3oldDMwLT);
   fChain->SetBranchAddress("Tau_byMediumIsolationMVA3oldDMwoLT", &Tau_byMediumIsolationMVA3oldDMwoLT, &b_Tau_byMediumIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("Tau_byTightIsolationMVA3newDMwLT", &Tau_byTightIsolationMVA3newDMwLT, &b_Tau_byTightIsolationMVA3newDMwLT);
   fChain->SetBranchAddress("Tau_byTightIsolationMVA3newDMwoLT", &Tau_byTightIsolationMVA3newDMwoLT, &b_Tau_byTightIsolationMVA3newDMwoLT);
   fChain->SetBranchAddress("Tau_byTightIsolationMva3oldDMwLT", &Tau_byTightIsolationMva3oldDMwLT, &b_Tau_byTightIsolationMva3oldDMwLT);
   fChain->SetBranchAddress("Tau_byTightIsolationMVA3oldDMwoLT", &Tau_byTightIsolationMVA3oldDMwoLT, &b_Tau_byTightIsolationMVA3oldDMwoLT);
   fChain->SetBranchAddress("Tau_againstMuonLoose2", &Tau_againstMuonLoose2, &b_Tau_againstMuonLoose2);
   fChain->SetBranchAddress("Tau_againstMuonLoose3", &Tau_againstMuonLoose3, &b_Tau_againstMuonLoose3);
   fChain->SetBranchAddress("Tau_againstMuonTight2", &Tau_againstMuonTight2, &b_Tau_againstMuonTight2);
   fChain->SetBranchAddress("Tau_againstElectronMVALooseMVA5", &Tau_againstElectronMVALooseMVA5, &b_Tau_againstElectronMVALooseMVA5);
   fChain->SetBranchAddress("Tau_againstElectronMVAMediumMVA5", &Tau_againstElectronMVAMediumMVA5, &b_Tau_againstElectronMVAMediumMVA5);
   fChain->SetBranchAddress("Tau_byVLooseCombinedIsolationDeltaBetaCorr", &Tau_byVLooseCombinedIsolationDeltaBetaCorr, &b_Tau_byVLooseCombinedIsolationDeltaBetaCorr);
   fChain->SetBranchAddress("Tau_leadChargedCandPt", &Tau_leadChargedCandPt, &b_Tau_leadChargedCandPt);
   fChain->SetBranchAddress("Tau_leadChargedCandCharge", &Tau_leadChargedCandCharge, &b_Tau_leadChargedCandCharge);
   fChain->SetBranchAddress("Tau_leadChargedCandEta", &Tau_leadChargedCandEta, &b_Tau_leadChargedCandEta);
   fChain->SetBranchAddress("Tau_leadChargedCandPhi", &Tau_leadChargedCandPhi, &b_Tau_leadChargedCandPhi);
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_energy", &Jet_energy, &b_Jet_energy);
   fChain->SetBranchAddress("Jet_bDiscriminator", &Jet_bDiscriminator, &b_Jet_bDiscriminator);
   fChain->SetBranchAddress("Jet_pileupId", &Jet_pileupId, &b_Jet_pileupId);
   fChain->SetBranchAddress("Jet_mass", &Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_neutralHadEnergyFraction", &Jet_neutralHadEnergyFraction, &b_Jet_neutralHadEnergyFraction);
   fChain->SetBranchAddress("Jet_neutralEmEmEnergyFraction", &Jet_neutralEmEmEnergyFraction, &b_Jet_neutralEmEmEnergyFraction);
   fChain->SetBranchAddress("Jet_chargedHadronEnergyFraction", &Jet_chargedHadronEnergyFraction, &b_Jet_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jet_chargedEmEnergyFraction", &Jet_chargedEmEnergyFraction, &b_Jet_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jet_muonEnergyFraction", &Jet_muonEnergyFraction, &b_Jet_muonEnergyFraction);
   fChain->SetBranchAddress("Jet_numberOfConstituents", &Jet_numberOfConstituents, &b_Jet_numberOfConstituents);
   fChain->SetBranchAddress("Jet_chargedMultiplicity", &Jet_chargedMultiplicity, &b_Jet_chargedMultiplicity);
   fChain->SetBranchAddress("Jet_electronEnergy", &Jet_electronEnergy, &b_Jet_electronEnergy);
   fChain->SetBranchAddress("Jet_photonEnergy", &Jet_photonEnergy, &b_Jet_photonEnergy);
   fChain->SetBranchAddress("UncorrJet_pt", &UncorrJet_pt, &b_UncorrJet_pt);
   fChain->SetBranchAddress("Photon_pt", &Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_eta", &Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_phi", &Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_energy", &Photon_energy, &b_Photon_energy);
   fChain->SetBranchAddress("Photon_et", &Photon_et, &b_Photon_et);
   fChain->SetBranchAddress("Photon_HoverE", &Photon_HoverE, &b_Photon_HoverE);
   fChain->SetBranchAddress("Photon_phoR9", &Photon_phoR9, &b_Photon_phoR9);
   fChain->SetBranchAddress("Photon_SigmaIEtaIEta", &Photon_SigmaIEtaIEta, &b_Photon_SigmaIEtaIEta);
   fChain->SetBranchAddress("Photon_SigmaIPhiIPhi", &Photon_SigmaIPhiIPhi, &b_Photon_SigmaIPhiIPhi);
   fChain->SetBranchAddress("Photon_PFChIso", &Photon_PFChIso, &b_Photon_PFChIso);
   fChain->SetBranchAddress("Photon_PFPhoIso", &Photon_PFPhoIso, &b_Photon_PFPhoIso);
   fChain->SetBranchAddress("Photon_PFNeuIso", &Photon_PFNeuIso, &b_Photon_PFNeuIso);
   fChain->SetBranchAddress("Photon_EleVeto", &Photon_EleVeto, &b_Photon_EleVeto);
   fChain->SetBranchAddress("Photon_hasPixelSeed", &Photon_hasPixelSeed, &b_Photon_hasPixelSeed);
   Notify();
}

Bool_t EventTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EventTree_cxx
