//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 23 15:27:41 2011 by ROOT version 5.28/00a
// from TTree TriggerTree/TriggerTree
// found on file: /msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root
//////////////////////////////////////////////////////////

#ifndef TriggerTree_h
#define TriggerTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class TriggerTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Bool_t          EF_g17_etcut;
   Bool_t          EF_e10_medium;
   Bool_t          EF_e15_medium;
   Bool_t          EF_e18_medium;
   Bool_t          EF_e20_medium;
   Bool_t          EF_e20_loose;
   Bool_t          EF_mu10_MSonly;
   Bool_t          EF_mu13;
   Bool_t          EF_mu13_tight;
   Bool_t          EF_mu15_tight;
   Bool_t          EF_mu18;
   Bool_t          EF_j20;
   Bool_t          EF_j20_jetNoCut;
   Bool_t          EF_j20_jetNoEF;
   Bool_t          EF_j30_jetNoCut;
   Bool_t          EF_j30_jetNoEF;
   Bool_t          EF_j35_jetNoCut;
   Bool_t          EF_j35_jetNoEF;
   Bool_t          EF_j50_jetNoCut;
   Bool_t          EF_j50_jetNoEF;
   UInt_t          RunNumber;
   UInt_t          EventNumber;
   UInt_t          lbn;
   Double_t        weight;

   // List of branches
   TBranch        *b_EF_g17_etcut;   //!
   TBranch        *b_EF_e10_medium;   //!
   TBranch        *b_EF_e15_medium;   //!
   TBranch        *b_EF_e18_medium;   //!
   TBranch        *b_EF_e20_medium;   //!
   TBranch        *b_EF_e20_loose;   //!
   TBranch        *b_EF_mu10_MSonly;   //!
   TBranch        *b_EF_mu13;   //!
   TBranch        *b_EF_mu13_tight;   //!
   TBranch        *b_EF_mu15_tight;   //!
   TBranch        *b_EF_mu18;   //!
   TBranch        *b_EF_j20;   //!
   TBranch        *b_EF_j20_jetNoCut;   //!
   TBranch        *b_EF_j20_jetNoEF;   //!
   TBranch        *b_EF_j30_jetNoCut;   //!
   TBranch        *b_EF_j30_jetNoEF;   //!
   TBranch        *b_EF_j35_jetNoCut;   //!
   TBranch        *b_EF_j35_jetNoEF;   //!
   TBranch        *b_EF_j50_jetNoCut;   //!
   TBranch        *b_EF_j50_jetNoEF;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_weight;   //!

   TriggerTree(TTree *tree=0);
   virtual ~TriggerTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TriggerTree_cxx
TriggerTree::TriggerTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root");
      if (!f) {
         f = new TFile("/msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root");
      }
      tree = (TTree*)gDirectory->Get("TriggerTree");

   }
   Init(tree);
}

TriggerTree::~TriggerTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TriggerTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TriggerTree::LoadTree(Long64_t entry)
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

void TriggerTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("EF_g17_etcut", &EF_g17_etcut, &b_EF_g17_etcut);
   fChain->SetBranchAddress("EF_e10_medium", &EF_e10_medium, &b_EF_e10_medium);
   fChain->SetBranchAddress("EF_e15_medium", &EF_e15_medium, &b_EF_e15_medium);
   fChain->SetBranchAddress("EF_e18_medium", &EF_e18_medium, &b_EF_e18_medium);
   fChain->SetBranchAddress("EF_e20_medium", &EF_e20_medium, &b_EF_e20_medium);
   fChain->SetBranchAddress("EF_e20_loose", &EF_e20_loose, &b_EF_e20_loose);
   fChain->SetBranchAddress("EF_mu10_MSonly", &EF_mu10_MSonly, &b_EF_mu10_MSonly);
   fChain->SetBranchAddress("EF_mu13", &EF_mu13, &b_EF_mu13);
   fChain->SetBranchAddress("EF_mu13_tight", &EF_mu13_tight, &b_EF_mu13_tight);
   fChain->SetBranchAddress("EF_mu15_tight", &EF_mu15_tight, &b_EF_mu15_tight);
   fChain->SetBranchAddress("EF_mu18", &EF_mu18, &b_EF_mu18);
   fChain->SetBranchAddress("EF_j20", &EF_j20, &b_EF_j20);
   fChain->SetBranchAddress("EF_j20_jetNoCut", &EF_j20_jetNoCut, &b_EF_j20_jetNoCut);
   fChain->SetBranchAddress("EF_j20_jetNoEF", &EF_j20_jetNoEF, &b_EF_j20_jetNoEF);
   fChain->SetBranchAddress("EF_j30_jetNoCut", &EF_j30_jetNoCut, &b_EF_j30_jetNoCut);
   fChain->SetBranchAddress("EF_j30_jetNoEF", &EF_j30_jetNoEF, &b_EF_j30_jetNoEF);
   fChain->SetBranchAddress("EF_j35_jetNoCut", &EF_j35_jetNoCut, &b_EF_j35_jetNoCut);
   fChain->SetBranchAddress("EF_j35_jetNoEF", &EF_j35_jetNoEF, &b_EF_j35_jetNoEF);
   fChain->SetBranchAddress("EF_j50_jetNoCut", &EF_j50_jetNoCut, &b_EF_j50_jetNoCut);
   fChain->SetBranchAddress("EF_j50_jetNoEF", &EF_j50_jetNoEF, &b_EF_j50_jetNoEF);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   Notify();
}

Bool_t TriggerTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TriggerTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TriggerTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TriggerTree_cxx
