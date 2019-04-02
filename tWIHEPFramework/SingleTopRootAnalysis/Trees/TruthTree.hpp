//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 14 14:03:31 2011 by ROOT version 5.28/00a
// from TTree TruthTree/TruthTree
// found on file: /msu/data/t3fast4/single_top/TopPhysDPDMaker/MC10bLP/user.xiaohu.mc10_7TeV.117367.st_Wtchan_dil_AcerMC.merge.AOD.e832_s933_s946_r2302_r2300_MC10b_2011JUL08.110712123253/user.xiaohu.003059.NTUP_SGTOP._00001.root
//////////////////////////////////////////////////////////

#ifndef TruthTree_h
#define TruthTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class TruthTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Tchannel_n;
   std::vector<int>     *Tchannel_pdgId;
   std::vector<int>     *Tchannel_barcode;
   std::vector<int>     *Tchannel_status;
   std::vector<int>     *Tchannel_nParents;
   std::vector<int>     *Tchannel_nDecay;
   std::vector<int>     *Tchannel_Type;
   std::vector<int>     *Tchannel_Origin;
   std::vector<double>  *Tchannel_m;
   std::vector<double>  *Tchannel_charge;
   std::vector<double>  *Tchannel_e;
   std::vector<double>  *Tchannel_et;
   std::vector<double>  *Tchannel_pt;
   std::vector<double>  *Tchannel_eta;
   std::vector<double>  *Tchannel_phi;
   std::vector<std::vector<double> > *Tchannel_Parent_pdgId;
   std::vector<std::vector<double> > *Tchannel_Parent_barcode;
   std::vector<std::vector<double> > *Tchannel_Parent_status;
   std::vector<std::vector<double> > *Tchannel_Parent_pt;
   std::vector<std::vector<double> > *Tchannel_Parent_eta;
   std::vector<std::vector<double> > *Tchannel_Parent_phi;
   std::vector<std::vector<double> > *Tchannel_Decay_pdgId;
   std::vector<std::vector<double> > *Tchannel_Decay_barcode;
   std::vector<std::vector<double> > *Tchannel_Decay_status;
   std::vector<std::vector<double> > *Tchannel_Decay_pt;
   std::vector<std::vector<double> > *Tchannel_Decay_eta;
   std::vector<std::vector<double> > *Tchannel_Decay_phi;
   Int_t           Wtchannel_n;
   std::vector<int>     *Wtchannel_pdgId;
   std::vector<int>     *Wtchannel_barcode;
   std::vector<int>     *Wtchannel_status;
   std::vector<int>     *Wtchannel_nParents;
   std::vector<int>     *Wtchannel_nDecay;
   std::vector<int>     *Wtchannel_Type;
   std::vector<int>     *Wtchannel_Origin;
   std::vector<double>  *Wtchannel_m;
   std::vector<double>  *Wtchannel_charge;
   std::vector<double>  *Wtchannel_e;
   std::vector<double>  *Wtchannel_et;
   std::vector<double>  *Wtchannel_pt;
   std::vector<double>  *Wtchannel_eta;
   std::vector<double>  *Wtchannel_phi;
   std::vector<std::vector<double> > *Wtchannel_Parent_pdgId;
   std::vector<std::vector<double> > *Wtchannel_Parent_barcode;
   std::vector<std::vector<double> > *Wtchannel_Parent_status;
   std::vector<std::vector<double> > *Wtchannel_Parent_pt;
   std::vector<std::vector<double> > *Wtchannel_Parent_eta;
   std::vector<std::vector<double> > *Wtchannel_Parent_phi;
   std::vector<std::vector<double> > *Wtchannel_Decay_pdgId;
   std::vector<std::vector<double> > *Wtchannel_Decay_barcode;
   std::vector<std::vector<double> > *Wtchannel_Decay_status;
   std::vector<std::vector<double> > *Wtchannel_Decay_pt;
   std::vector<std::vector<double> > *Wtchannel_Decay_eta;
   std::vector<std::vector<double> > *Wtchannel_Decay_phi;
   Int_t           hfor_type;
   UInt_t          RunNumber;
   UInt_t          EventNumber;
   UInt_t          lbn;
   Double_t        weight;

   // List of branches
   TBranch        *b_Tchannel_n;   //!
   TBranch        *b_Tchannel_pdgId;   //!
   TBranch        *b_Tchannel_barcode;   //!
   TBranch        *b_Tchannel_status;   //!
   TBranch        *b_Tchannel_nParents;   //!
   TBranch        *b_Tchannel_nDecay;   //!
   TBranch        *b_Tchannel_Type;   //!
   TBranch        *b_Tchannel_Origin;   //!
   TBranch        *b_Tchannel_m;   //!
   TBranch        *b_Tchannel_charge;   //!
   TBranch        *b_Tchannel_e;   //!
   TBranch        *b_Tchannel_et;   //!
   TBranch        *b_Tchannel_pt;   //!
   TBranch        *b_Tchannel_eta;   //!
   TBranch        *b_Tchannel_phi;   //!
   TBranch        *b_Tchannel_Parent_pdgId;   //!
   TBranch        *b_Tchannel_Parent_barcode;   //!
   TBranch        *b_Tchannel_Parent_status;   //!
   TBranch        *b_Tchannel_Parent_pt;   //!
   TBranch        *b_Tchannel_Parent_eta;   //!
   TBranch        *b_Tchannel_Parent_phi;   //!
   TBranch        *b_Tchannel_Decay_pdgId;   //!
   TBranch        *b_Tchannel_Decay_barcode;   //!
   TBranch        *b_Tchannel_Decay_status;   //!
   TBranch        *b_Tchannel_Decay_pt;   //!
   TBranch        *b_Tchannel_Decay_eta;   //!
   TBranch        *b_Tchannel_Decay_phi;   //!
   TBranch        *b_Wtchannel_n;   //!
   TBranch        *b_Wtchannel_pdgId;   //!
   TBranch        *b_Wtchannel_barcode;   //!
   TBranch        *b_Wtchannel_status;   //!
   TBranch        *b_Wtchannel_nParents;   //!
   TBranch        *b_Wtchannel_nDecay;   //!
   TBranch        *b_Wtchannel_Type;   //!
   TBranch        *b_Wtchannel_Origin;   //!
   TBranch        *b_Wtchannel_m;   //!
   TBranch        *b_Wtchannel_charge;   //!
   TBranch        *b_Wtchannel_e;   //!
   TBranch        *b_Wtchannel_et;   //!
   TBranch        *b_Wtchannel_pt;   //!
   TBranch        *b_Wtchannel_eta;   //!
   TBranch        *b_Wtchannel_phi;   //!
   TBranch        *b_Wtchannel_Parent_pdgId;   //!
   TBranch        *b_Wtchannel_Parent_barcode;   //!
   TBranch        *b_Wtchannel_Parent_status;   //!
   TBranch        *b_Wtchannel_Parent_pt;   //!
   TBranch        *b_Wtchannel_Parent_eta;   //!
   TBranch        *b_Wtchannel_Parent_phi;   //!
   TBranch        *b_Wtchannel_Decay_pdgId;   //!
   TBranch        *b_Wtchannel_Decay_barcode;   //!
   TBranch        *b_Wtchannel_Decay_status;   //!
   TBranch        *b_Wtchannel_Decay_pt;   //!
   TBranch        *b_Wtchannel_Decay_eta;   //!
   TBranch        *b_Wtchannel_Decay_phi;   //!
   TBranch        *b_hfor_type;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_weight;   //!

   TruthTree(TTree *tree=0);
   virtual ~TruthTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TruthTree_cxx
TruthTree::TruthTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/msu/data/t3fast4/single_top/TopPhysDPDMaker/MC10bLP/user.xiaohu.mc10_7TeV.117367.st_Wtchan_dil_AcerMC.merge.AOD.e832_s933_s946_r2302_r2300_MC10b_2011JUL08.110712123253/user.xiaohu.003059.NTUP_SGTOP._00001.root");
      if (!f) {
         f = new TFile("/msu/data/t3fast4/single_top/TopPhysDPDMaker/MC10bLP/user.xiaohu.mc10_7TeV.117367.st_Wtchan_dil_AcerMC.merge.AOD.e832_s933_s946_r2302_r2300_MC10b_2011JUL08.110712123253/user.xiaohu.003059.NTUP_SGTOP._00001.root");
      }
      tree = (TTree*)gDirectory->Get("TruthTree");

   }
   Init(tree);
}

TruthTree::~TruthTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TruthTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TruthTree::LoadTree(Long64_t entry)
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

void TruthTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Tchannel_pdgId = 0;
   Tchannel_barcode = 0;
   Tchannel_status = 0;
   Tchannel_nParents = 0;
   Tchannel_nDecay = 0;
   Tchannel_Type = 0;
   Tchannel_Origin = 0;
   Tchannel_m = 0;
   Tchannel_charge = 0;
   Tchannel_e = 0;
   Tchannel_et = 0;
   Tchannel_pt = 0;
   Tchannel_eta = 0;
   Tchannel_phi = 0;
   Tchannel_Parent_pdgId = 0;
   Tchannel_Parent_barcode = 0;
   Tchannel_Parent_status = 0;
   Tchannel_Parent_pt = 0;
   Tchannel_Parent_eta = 0;
   Tchannel_Parent_phi = 0;
   Tchannel_Decay_pdgId = 0;
   Tchannel_Decay_barcode = 0;
   Tchannel_Decay_status = 0;
   Tchannel_Decay_pt = 0;
   Tchannel_Decay_eta = 0;
   Tchannel_Decay_phi = 0;
   Wtchannel_pdgId = 0;
   Wtchannel_barcode = 0;
   Wtchannel_status = 0;
   Wtchannel_nParents = 0;
   Wtchannel_nDecay = 0;
   Wtchannel_Type = 0;
   Wtchannel_Origin = 0;
   Wtchannel_m = 0;
   Wtchannel_charge = 0;
   Wtchannel_e = 0;
   Wtchannel_et = 0;
   Wtchannel_pt = 0;
   Wtchannel_eta = 0;
   Wtchannel_phi = 0;
   Wtchannel_Parent_pdgId = 0;
   Wtchannel_Parent_barcode = 0;
   Wtchannel_Parent_status = 0;
   Wtchannel_Parent_pt = 0;
   Wtchannel_Parent_eta = 0;
   Wtchannel_Parent_phi = 0;
   Wtchannel_Decay_pdgId = 0;
   Wtchannel_Decay_barcode = 0;
   Wtchannel_Decay_status = 0;
   Wtchannel_Decay_pt = 0;
   Wtchannel_Decay_eta = 0;
   Wtchannel_Decay_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Tchannel_n", &Tchannel_n, &b_Tchannel_n);
   fChain->SetBranchAddress("Tchannel_pdgId", &Tchannel_pdgId, &b_Tchannel_pdgId);
   fChain->SetBranchAddress("Tchannel_barcode", &Tchannel_barcode, &b_Tchannel_barcode);
   fChain->SetBranchAddress("Tchannel_status", &Tchannel_status, &b_Tchannel_status);
   fChain->SetBranchAddress("Tchannel_nParents", &Tchannel_nParents, &b_Tchannel_nParents);
   fChain->SetBranchAddress("Tchannel_nDecay", &Tchannel_nDecay, &b_Tchannel_nDecay);
   fChain->SetBranchAddress("Tchannel_Type", &Tchannel_Type, &b_Tchannel_Type);
   fChain->SetBranchAddress("Tchannel_Origin", &Tchannel_Origin, &b_Tchannel_Origin);
   fChain->SetBranchAddress("Tchannel_m", &Tchannel_m, &b_Tchannel_m);
   fChain->SetBranchAddress("Tchannel_charge", &Tchannel_charge, &b_Tchannel_charge);
   fChain->SetBranchAddress("Tchannel_e", &Tchannel_e, &b_Tchannel_e);
   fChain->SetBranchAddress("Tchannel_et", &Tchannel_et, &b_Tchannel_et);
   fChain->SetBranchAddress("Tchannel_pt", &Tchannel_pt, &b_Tchannel_pt);
   fChain->SetBranchAddress("Tchannel_eta", &Tchannel_eta, &b_Tchannel_eta);
   fChain->SetBranchAddress("Tchannel_phi", &Tchannel_phi, &b_Tchannel_phi);
   fChain->SetBranchAddress("Tchannel_Parent_pdgId", &Tchannel_Parent_pdgId, &b_Tchannel_Parent_pdgId);
   fChain->SetBranchAddress("Tchannel_Parent_barcode", &Tchannel_Parent_barcode, &b_Tchannel_Parent_barcode);
   fChain->SetBranchAddress("Tchannel_Parent_status", &Tchannel_Parent_status, &b_Tchannel_Parent_status);
   fChain->SetBranchAddress("Tchannel_Parent_pt", &Tchannel_Parent_pt, &b_Tchannel_Parent_pt);
   fChain->SetBranchAddress("Tchannel_Parent_eta", &Tchannel_Parent_eta, &b_Tchannel_Parent_eta);
   fChain->SetBranchAddress("Tchannel_Parent_phi", &Tchannel_Parent_phi, &b_Tchannel_Parent_phi);
   fChain->SetBranchAddress("Tchannel_Decay_pdgId", &Tchannel_Decay_pdgId, &b_Tchannel_Decay_pdgId);
   fChain->SetBranchAddress("Tchannel_Decay_barcode", &Tchannel_Decay_barcode, &b_Tchannel_Decay_barcode);
   fChain->SetBranchAddress("Tchannel_Decay_status", &Tchannel_Decay_status, &b_Tchannel_Decay_status);
   fChain->SetBranchAddress("Tchannel_Decay_pt", &Tchannel_Decay_pt, &b_Tchannel_Decay_pt);
   fChain->SetBranchAddress("Tchannel_Decay_eta", &Tchannel_Decay_eta, &b_Tchannel_Decay_eta);
   fChain->SetBranchAddress("Tchannel_Decay_phi", &Tchannel_Decay_phi, &b_Tchannel_Decay_phi);
   fChain->SetBranchAddress("Wtchannel_n", &Wtchannel_n, &b_Wtchannel_n);
   fChain->SetBranchAddress("Wtchannel_pdgId", &Wtchannel_pdgId, &b_Wtchannel_pdgId);
   fChain->SetBranchAddress("Wtchannel_barcode", &Wtchannel_barcode, &b_Wtchannel_barcode);
   fChain->SetBranchAddress("Wtchannel_status", &Wtchannel_status, &b_Wtchannel_status);
   fChain->SetBranchAddress("Wtchannel_nParents", &Wtchannel_nParents, &b_Wtchannel_nParents);
   fChain->SetBranchAddress("Wtchannel_nDecay", &Wtchannel_nDecay, &b_Wtchannel_nDecay);
   fChain->SetBranchAddress("Wtchannel_Type", &Wtchannel_Type, &b_Wtchannel_Type);
   fChain->SetBranchAddress("Wtchannel_Origin", &Wtchannel_Origin, &b_Wtchannel_Origin);
   fChain->SetBranchAddress("Wtchannel_m", &Wtchannel_m, &b_Wtchannel_m);
   fChain->SetBranchAddress("Wtchannel_charge", &Wtchannel_charge, &b_Wtchannel_charge);
   fChain->SetBranchAddress("Wtchannel_e", &Wtchannel_e, &b_Wtchannel_e);
   fChain->SetBranchAddress("Wtchannel_et", &Wtchannel_et, &b_Wtchannel_et);
   fChain->SetBranchAddress("Wtchannel_pt", &Wtchannel_pt, &b_Wtchannel_pt);
   fChain->SetBranchAddress("Wtchannel_eta", &Wtchannel_eta, &b_Wtchannel_eta);
   fChain->SetBranchAddress("Wtchannel_phi", &Wtchannel_phi, &b_Wtchannel_phi);
   fChain->SetBranchAddress("Wtchannel_Parent_pdgId", &Wtchannel_Parent_pdgId, &b_Wtchannel_Parent_pdgId);
   fChain->SetBranchAddress("Wtchannel_Parent_barcode", &Wtchannel_Parent_barcode, &b_Wtchannel_Parent_barcode);
   fChain->SetBranchAddress("Wtchannel_Parent_status", &Wtchannel_Parent_status, &b_Wtchannel_Parent_status);
   fChain->SetBranchAddress("Wtchannel_Parent_pt", &Wtchannel_Parent_pt, &b_Wtchannel_Parent_pt);
   fChain->SetBranchAddress("Wtchannel_Parent_eta", &Wtchannel_Parent_eta, &b_Wtchannel_Parent_eta);
   fChain->SetBranchAddress("Wtchannel_Parent_phi", &Wtchannel_Parent_phi, &b_Wtchannel_Parent_phi);
   fChain->SetBranchAddress("Wtchannel_Decay_pdgId", &Wtchannel_Decay_pdgId, &b_Wtchannel_Decay_pdgId);
   fChain->SetBranchAddress("Wtchannel_Decay_barcode", &Wtchannel_Decay_barcode, &b_Wtchannel_Decay_barcode);
   fChain->SetBranchAddress("Wtchannel_Decay_status", &Wtchannel_Decay_status, &b_Wtchannel_Decay_status);
   fChain->SetBranchAddress("Wtchannel_Decay_pt", &Wtchannel_Decay_pt, &b_Wtchannel_Decay_pt);
   fChain->SetBranchAddress("Wtchannel_Decay_eta", &Wtchannel_Decay_eta, &b_Wtchannel_Decay_eta);
   fChain->SetBranchAddress("Wtchannel_Decay_phi", &Wtchannel_Decay_phi, &b_Wtchannel_Decay_phi);
   fChain->SetBranchAddress("hfor_type", &hfor_type, &b_hfor_type);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   Notify();
}

Bool_t TruthTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TruthTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TruthTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TruthTree_cxx
