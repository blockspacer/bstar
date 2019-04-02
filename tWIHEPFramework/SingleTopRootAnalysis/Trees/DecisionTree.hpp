//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 23 15:27:29 2011 by ROOT version 5.28/00a
// from TTree DecisionTree/DecisionTree
// found on file: /msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root
//////////////////////////////////////////////////////////

#ifndef DecisionTree_h
#define DecisionTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class DecisionTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          RunNumber;
   UInt_t          EventNumber;
   UInt_t          lbn;
   Double_t        weight;
   Bool_t          ElectronPreTag;
   Bool_t          MuonPreTag;
   Bool_t          ElectronLoose;
   Bool_t          MuonLoose;
   Bool_t          DiElectronsPreTag;
   Bool_t          DiMuonsPreTag;
   Bool_t          ElectronMuonPreTag;
   Bool_t          DiElectronsLoose;
   Bool_t          DiMuonsLoose;
   Bool_t          ElectronMuonLoose;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_ElectronPreTag;   //!
   TBranch        *b_MuonPreTag;   //!
   TBranch        *b_ElectronLoose;   //!
   TBranch        *b_MuonLoose;   //!
   TBranch        *b_DiElectronsPreTag;   //!
   TBranch        *b_DiMuonsPreTag;   //!
   TBranch        *b_ElectronMuonPreTag;   //!
   TBranch        *b_DiElectronsLoose;   //!
   TBranch        *b_DiMuonsLoose;   //!
   TBranch        *b_ElectronMuonLoose;   //!

   DecisionTree(TTree *tree=0);
   virtual ~DecisionTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DecisionTree_cxx
DecisionTree::DecisionTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root");
      if (!f) {
         f = new TFile("/msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root");
      }
      tree = (TTree*)gDirectory->Get("DecisionTree");

   }
   Init(tree);
}

DecisionTree::~DecisionTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DecisionTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DecisionTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DecisionTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("ElectronPreTag", &ElectronPreTag, &b_ElectronPreTag);
   fChain->SetBranchAddress("MuonPreTag", &MuonPreTag, &b_MuonPreTag);
   fChain->SetBranchAddress("ElectronLoose", &ElectronLoose, &b_ElectronLoose);
   fChain->SetBranchAddress("MuonLoose", &MuonLoose, &b_MuonLoose);
   fChain->SetBranchAddress("DiElectronsPreTag", &DiElectronsPreTag, &b_DiElectronsPreTag);
   fChain->SetBranchAddress("DiMuonsPreTag", &DiMuonsPreTag, &b_DiMuonsPreTag);
   fChain->SetBranchAddress("ElectronMuonPreTag", &ElectronMuonPreTag, &b_ElectronMuonPreTag);
   fChain->SetBranchAddress("DiElectronsLoose", &DiElectronsLoose, &b_DiElectronsLoose);
   fChain->SetBranchAddress("DiMuonsLoose", &DiMuonsLoose, &b_DiMuonsLoose);
   fChain->SetBranchAddress("ElectronMuonLoose", &ElectronMuonLoose, &b_ElectronMuonLoose);
   Notify();
}

Bool_t DecisionTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DecisionTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DecisionTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DecisionTree_cxx
