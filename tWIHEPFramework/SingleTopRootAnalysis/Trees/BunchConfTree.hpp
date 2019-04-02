//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep  1 11:27:48 2011 by ROOT version 5.28/00a
// from TTree BunchConfTree/BunchConfTree
// found on file: /msu/data/t3fast4/single_top/TopPhysDPDMaker/Data10bLP//user.csuhr.v16.6.6.2.2.data11_7TeV.00182886.physics_Egamma.SGTOP_NTUP_TOPEL.f382_m866_p568_p569_v1.110708170256/user.csuhr.004844.NTUP_SGTOP._00001.root
//////////////////////////////////////////////////////////

#ifndef BunchConfTree_h
#define BunchConfTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class BunchConfTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          ConfigID;
   std::vector<int>     *FilledBunches;
   std::vector<float>   *BunchIntensities;

   // List of branches
   TBranch        *b_ConfigID;   //!
   TBranch        *b_FilledBunches;   //!
   TBranch        *b_BunchIntensities;   //!

   BunchConfTree(TTree *tree=0);
   virtual ~BunchConfTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BunchConfTree_cxx
BunchConfTree::BunchConfTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/msu/data/t3fast4/single_top/TopPhysDPDMaker/Data10bLP//user.csuhr.v16.6.6.2.2.data11_7TeV.00182886.physics_Egamma.SGTOP_NTUP_TOPEL.f382_m866_p568_p569_v1.110708170256/user.csuhr.004844.NTUP_SGTOP._00001.root");
      if (!f) {
         f = new TFile("/msu/data/t3fast4/single_top/TopPhysDPDMaker/Data10bLP//user.csuhr.v16.6.6.2.2.data11_7TeV.00182886.physics_Egamma.SGTOP_NTUP_TOPEL.f382_m866_p568_p569_v1.110708170256/user.csuhr.004844.NTUP_SGTOP._00001.root");
      }
      tree = (TTree*)gDirectory->Get("BunchConfTree");

   }
   Init(tree);
}

BunchConfTree::~BunchConfTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BunchConfTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BunchConfTree::LoadTree(Long64_t entry)
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

void BunchConfTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   FilledBunches = 0;
   BunchIntensities = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ConfigID", &ConfigID, &b_ConfigID);
   fChain->SetBranchAddress("FilledBunches", &FilledBunches, &b_FilledBunches);
   fChain->SetBranchAddress("BunchIntensities", &BunchIntensities, &b_BunchIntensities);
   Notify();
}

Bool_t BunchConfTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BunchConfTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BunchConfTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BunchConfTree_cxx
