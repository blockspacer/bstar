//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 23 15:27:35 2011 by ROOT version 5.28/00a
// from TTree ElectronMuonPreTagTree/ElectronMuonPreTagTree
// found on file: /msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root
//////////////////////////////////////////////////////////

#ifndef ElectronMuonPreTagTree_h
#define ElectronMuonPreTagTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class ElectronMuonPreTagTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          Lepton1_loose;
   UInt_t          Lepton1_medium;
   UInt_t          Lepton1_mediumIso;
   UInt_t          Lepton1_tight;
   UInt_t          Lepton1_tightIso;
   UInt_t          Lepton1_trackBlayer;
   UInt_t          Lepton1_trackMatchEoverP;
   UInt_t          Lepton1_trackTRThits;
   UInt_t          Lepton1_trackTRTratio;
   UInt_t          Lepton1_trackA0Tight;
   UInt_t          Lepton1_trackMatchEtaTight;
   UInt_t          Lepton1_trackMatchPhi;
   UInt_t          Lepton1_conversionMatch;
   Int_t           Lepton1_Barcode;
   Int_t           Lepton1_pdgId;
   Float_t         Lepton1_charge;
   Float_t         Lepton1_Etcone45;
   Float_t         Lepton1_Etcone20;
   Float_t         Lepton1_Etcone30;
   Float_t         Lepton1_Etcone40;
   Float_t         Lepton1_ptcone20;
   Float_t         Lepton1_ptcone30;
   Float_t         Lepton1_ptcone40;
   Float_t         Lepton1_cluster_pt;
   Float_t         Lepton1_cluster_eta;
   Float_t         Lepton1_cluster_phi;
   Float_t         Lepton1_cluster_E;
   Float_t         Lepton1_cluster_Et;
   Float_t         Lepton1_cluster_track_Et;
   Float_t         Lepton1_track_d0;
   Float_t         Lepton1_track_z0;
   Float_t         Lepton1_track_phi;
   Float_t         Lepton1_track_theta;
   Float_t         Lepton1_track_qoverp;
   Float_t         Lepton1_track_pt;
   Float_t         Lepton1_track_eta;
   Int_t           Lepton2_Barcode;
   Int_t           Lepton2_pdgId;
   Int_t           Lepton2_isGoodIDTrack;
   Float_t         Lepton2_pt;
   Float_t         Lepton2_eta;
   Float_t         Lepton2_phi;
   Float_t         Lepton2_E;
   Float_t         Lepton2_Et;
   Float_t         Lepton2_charge;
   Float_t         Lepton2_etcone20;
   Float_t         Lepton2_etcone30;
   Float_t         Lepton2_etcone40;
   Float_t         Lepton2_ptcone20;
   Float_t         Lepton2_ptcone30;
   Float_t         Lepton2_ptcone40;
   Float_t         Lepton2_id_pt;
   Float_t         Lepton2_id_eta;
   Float_t         Lepton2_me_pt;
   Float_t         Lepton2_me_eta;
   Float_t         MissingEt_etx;
   Float_t         MissingEt_ety;
   Float_t         MissingEt_et;
   Float_t         MissingEt_sumet;
   Float_t         MissingEt_phi;
   Int_t           Jet_n;
   std::vector<int>     *Jet_Barcode;
   std::vector<int>     *Jet_flavor_truth_trueflav;
   std::vector<int>     *Jet_flavor_truth_BHadronpdg;
   std::vector<float>   *Jet_pt;
   std::vector<float>   *Jet_eta;
   std::vector<float>   *Jet_phi;
   std::vector<float>   *Jet_E;
   std::vector<float>   *Jet_Et;
   std::vector<float>   *Jet_emscale_pt;
   std::vector<float>   *Jet_emscale_eta;
   std::vector<float>   *Jet_emscale_phi;
   std::vector<float>   *Jet_emscale_E;
   std::vector<double>  *Jet_flavor_weight;
   std::vector<double>  *Jet_flavor_weight_TrackCounting2D;
   std::vector<double>  *Jet_flavor_weight_JetProb;
   std::vector<double>  *Jet_flavor_weight_IP1D;
   std::vector<double>  *Jet_flavor_weight_IP2D;
   std::vector<double>  *Jet_flavor_weight_IP3D;
   std::vector<double>  *Jet_flavor_weight_SV0;
   std::vector<double>  *Jet_flavor_weight_SV1;
   std::vector<double>  *Jet_flavor_weight_SV2;
   std::vector<double>  *Jet_flavor_weight_JetFitterTag;
   std::vector<double>  *Jet_flavor_weight_JetFitterCOMB;
   std::vector<double>  *Jet_flavor_weight_JetFitterCOMBNN;
   std::vector<double>  *Jet_flavor_weight_SoftMuonTag;
   std::vector<double>  *Jet_flavor_weight_SoftElectronTag;
   std::vector<double>  *Jet_flavor_truth_dRMinToB;
   std::vector<double>  *Jet_flavor_truth_dRMinToC;
   std::vector<double>  *Jet_flavor_truth_dRMinToT;
   Float_t         DiLeptons_m;
   Float_t         DiLeptons_E;
   Float_t         DiLeptons_Et;
   Float_t         DiLeptons_pt;
   Float_t         DiLeptons_eta;
   Float_t         DiLeptons_phi;
   Float_t         DiLeptons_charge;
   Int_t           BJet_n;
   std::vector<int>     *BJet_Barcode;
   std::vector<int>     *BJet_flavor_truth_trueflav;
   std::vector<int>     *BJet_flavor_truth_BHadronpdg;
   std::vector<float>   *BJet_pt;
   std::vector<float>   *BJet_eta;
   std::vector<float>   *BJet_phi;
   std::vector<float>   *BJet_E;
   std::vector<float>   *BJet_Et;
   std::vector<float>   *BJet_emscale_pt;
   std::vector<float>   *BJet_emscale_eta;
   std::vector<float>   *BJet_emscale_phi;
   std::vector<float>   *BJet_emscale_E;
   std::vector<double>  *BJet_flavor_weight;
   std::vector<double>  *BJet_flavor_weight_TrackCounting2D;
   std::vector<double>  *BJet_flavor_weight_JetProb;
   std::vector<double>  *BJet_flavor_weight_IP1D;
   std::vector<double>  *BJet_flavor_weight_IP2D;
   std::vector<double>  *BJet_flavor_weight_IP3D;
   std::vector<double>  *BJet_flavor_weight_SV0;
   std::vector<double>  *BJet_flavor_weight_SV1;
   std::vector<double>  *BJet_flavor_weight_SV2;
   std::vector<double>  *BJet_flavor_weight_JetFitterTag;
   std::vector<double>  *BJet_flavor_weight_JetFitterCOMB;
   std::vector<double>  *BJet_flavor_weight_JetFitterCOMBNN;
   std::vector<double>  *BJet_flavor_weight_SoftMuonTag;
   std::vector<double>  *BJet_flavor_weight_SoftElectronTag;
   std::vector<double>  *BJet_flavor_truth_dRMinToB;
   std::vector<double>  *BJet_flavor_truth_dRMinToC;
   std::vector<double>  *BJet_flavor_truth_dRMinToT;
   std::vector<int>     *BJet_index;
   UInt_t          RunNumber;
   UInt_t          EventNumber;
   UInt_t          lbn;
   Double_t        weight;
   Bool_t          Lepton1_inTrigger;
   Bool_t          Lepton2_inTrigger;
   Bool_t          isElectronTriggerPassed;
   Bool_t          isMuonTriggerPassed;
   Int_t           VetoJet_n;
   std::vector<int>     *VetoJet_Barcode;
   std::vector<int>     *VetoJet_flavor_truth_trueflav;
   std::vector<int>     *VetoJet_flavor_truth_BHadronpdg;
   std::vector<float>   *VetoJet_pt;
   std::vector<float>   *VetoJet_eta;
   std::vector<float>   *VetoJet_phi;
   std::vector<float>   *VetoJet_E;
   std::vector<float>   *VetoJet_Et;
   std::vector<float>   *VetoJet_emscale_pt;
   std::vector<float>   *VetoJet_emscale_eta;
   std::vector<float>   *VetoJet_emscale_phi;
   std::vector<float>   *VetoJet_emscale_E;
   std::vector<double>  *VetoJet_flavor_weight;
   std::vector<double>  *VetoJet_flavor_weight_TrackCounting2D;
   std::vector<double>  *VetoJet_flavor_weight_JetProb;
   std::vector<double>  *VetoJet_flavor_weight_IP1D;
   std::vector<double>  *VetoJet_flavor_weight_IP2D;
   std::vector<double>  *VetoJet_flavor_weight_IP3D;
   std::vector<double>  *VetoJet_flavor_weight_SV0;
   std::vector<double>  *VetoJet_flavor_weight_SV1;
   std::vector<double>  *VetoJet_flavor_weight_SV2;
   std::vector<double>  *VetoJet_flavor_weight_JetFitterTag;
   std::vector<double>  *VetoJet_flavor_weight_JetFitterCOMB;
   std::vector<double>  *VetoJet_flavor_weight_JetFitterCOMBNN;
   std::vector<double>  *VetoJet_flavor_weight_SoftMuonTag;
   std::vector<double>  *VetoJet_flavor_weight_SoftElectronTag;
   std::vector<double>  *VetoJet_flavor_truth_dRMinToB;
   std::vector<double>  *VetoJet_flavor_truth_dRMinToC;
   std::vector<double>  *VetoJet_flavor_truth_dRMinToT;
   Bool_t          isNoEleMuOverlap;

   // List of branches
   TBranch        *b_Lepton1_loose;   //!
   TBranch        *b_Lepton1_medium;   //!
   TBranch        *b_Lepton1_mediumIso;   //!
   TBranch        *b_Lepton1_tight;   //!
   TBranch        *b_Lepton1_tightIso;   //!
   TBranch        *b_Lepton1_trackBlayer;   //!
   TBranch        *b_Lepton1_trackMatchEoverP;   //!
   TBranch        *b_Lepton1_trackTRThits;   //!
   TBranch        *b_Lepton1_trackTRTratio;   //!
   TBranch        *b_Lepton1_trackA0Tight;   //!
   TBranch        *b_Lepton1_trackMatchEtaTight;   //!
   TBranch        *b_Lepton1_trackMatchPhi;   //!
   TBranch        *b_Lepton1_conversionMatch;   //!
   TBranch        *b_Lepton1_Barcode;   //!
   TBranch        *b_Lepton1_pdgId;   //!
   TBranch        *b_Lepton1_charge;   //!
   TBranch        *b_Lepton1_Etcone45;   //!
   TBranch        *b_Lepton1_Etcone20;   //!
   TBranch        *b_Lepton1_Etcone30;   //!
   TBranch        *b_Lepton1_Etcone40;   //!
   TBranch        *b_Lepton1_ptcone20;   //!
   TBranch        *b_Lepton1_ptcone30;   //!
   TBranch        *b_Lepton1_ptcone40;   //!
   TBranch        *b_Lepton1_cluster_pt;   //!
   TBranch        *b_Lepton1_cluster_eta;   //!
   TBranch        *b_Lepton1_cluster_phi;   //!
   TBranch        *b_Lepton1_cluster_E;   //!
   TBranch        *b_Lepton1_cluster_Et;   //!
   TBranch        *b_Lepton1_cluster_track_Et;   //!
   TBranch        *b_Lepton1_track_d0;   //!
   TBranch        *b_Lepton1_track_z0;   //!
   TBranch        *b_Lepton1_track_phi;   //!
   TBranch        *b_Lepton1_track_theta;   //!
   TBranch        *b_Lepton1_track_qoverp;   //!
   TBranch        *b_Lepton1_track_pt;   //!
   TBranch        *b_Lepton1_track_eta;   //!
   TBranch        *b_Lepton2_Barcode;   //!
   TBranch        *b_Lepton2_pdgId;   //!
   TBranch        *b_Lepton2_isGoodIDTrack;   //!
   TBranch        *b_Lepton2_pt;   //!
   TBranch        *b_Lepton2_eta;   //!
   TBranch        *b_Lepton2_phi;   //!
   TBranch        *b_Lepton2_E;   //!
   TBranch        *b_Lepton2_Et;   //!
   TBranch        *b_Lepton2_charge;   //!
   TBranch        *b_Lepton2_etcone20;   //!
   TBranch        *b_Lepton2_etcone30;   //!
   TBranch        *b_Lepton2_etcone40;   //!
   TBranch        *b_Lepton2_ptcone20;   //!
   TBranch        *b_Lepton2_ptcone30;   //!
   TBranch        *b_Lepton2_ptcone40;   //!
   TBranch        *b_Lepton2_id_pt;   //!
   TBranch        *b_Lepton2_id_eta;   //!
   TBranch        *b_Lepton2_me_pt;   //!
   TBranch        *b_Lepton2_me_eta;   //!
   TBranch        *b_MissingEt_etx;   //!
   TBranch        *b_MissingEt_ety;   //!
   TBranch        *b_MissingEt_et;   //!
   TBranch        *b_MissingEt_sumet;   //!
   TBranch        *b_MissingEt_phi;   //!
   TBranch        *b_Jet_n;   //!
   TBranch        *b_Jet_Barcode;   //!
   TBranch        *b_Jet_flavor_truth_trueflav;   //!
   TBranch        *b_Jet_flavor_truth_BHadronpdg;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_E;   //!
   TBranch        *b_Jet_Et;   //!
   TBranch        *b_Jet_emscale_pt;   //!
   TBranch        *b_Jet_emscale_eta;   //!
   TBranch        *b_Jet_emscale_phi;   //!
   TBranch        *b_Jet_emscale_E;   //!
   TBranch        *b_Jet_flavor_weight;   //!
   TBranch        *b_Jet_flavor_weight_TrackCounting2D;   //!
   TBranch        *b_Jet_flavor_weight_JetProb;   //!
   TBranch        *b_Jet_flavor_weight_IP1D;   //!
   TBranch        *b_Jet_flavor_weight_IP2D;   //!
   TBranch        *b_Jet_flavor_weight_IP3D;   //!
   TBranch        *b_Jet_flavor_weight_SV0;   //!
   TBranch        *b_Jet_flavor_weight_SV1;   //!
   TBranch        *b_Jet_flavor_weight_SV2;   //!
   TBranch        *b_Jet_flavor_weight_JetFitterTag;   //!
   TBranch        *b_Jet_flavor_weight_JetFitterCOMB;   //!
   TBranch        *b_Jet_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_Jet_flavor_weight_SoftMuonTag;   //!
   TBranch        *b_Jet_flavor_weight_SoftElectronTag;   //!
   TBranch        *b_Jet_flavor_truth_dRMinToB;   //!
   TBranch        *b_Jet_flavor_truth_dRMinToC;   //!
   TBranch        *b_Jet_flavor_truth_dRMinToT;   //!
   TBranch        *b_DiLeptons_m;   //!
   TBranch        *b_DiLeptons_E;   //!
   TBranch        *b_DiLeptons_Et;   //!
   TBranch        *b_DiLeptons_pt;   //!
   TBranch        *b_DiLeptons_eta;   //!
   TBranch        *b_DiLeptons_phi;   //!
   TBranch        *b_DiLeptons_charge;   //!
   TBranch        *b_BJet_n;   //!
   TBranch        *b_BJet_Barcode;   //!
   TBranch        *b_BJet_flavor_truth_trueflav;   //!
   TBranch        *b_BJet_flavor_truth_BHadronpdg;   //!
   TBranch        *b_BJet_pt;   //!
   TBranch        *b_BJet_eta;   //!
   TBranch        *b_BJet_phi;   //!
   TBranch        *b_BJet_E;   //!
   TBranch        *b_BJet_Et;   //!
   TBranch        *b_BJet_emscale_pt;   //!
   TBranch        *b_BJet_emscale_eta;   //!
   TBranch        *b_BJet_emscale_phi;   //!
   TBranch        *b_BJet_emscale_E;   //!
   TBranch        *b_BJet_flavor_weight;   //!
   TBranch        *b_BJet_flavor_weight_TrackCounting2D;   //!
   TBranch        *b_BJet_flavor_weight_JetProb;   //!
   TBranch        *b_BJet_flavor_weight_IP1D;   //!
   TBranch        *b_BJet_flavor_weight_IP2D;   //!
   TBranch        *b_BJet_flavor_weight_IP3D;   //!
   TBranch        *b_BJet_flavor_weight_SV0;   //!
   TBranch        *b_BJet_flavor_weight_SV1;   //!
   TBranch        *b_BJet_flavor_weight_SV2;   //!
   TBranch        *b_BJet_flavor_weight_JetFitterTag;   //!
   TBranch        *b_BJet_flavor_weight_JetFitterCOMB;   //!
   TBranch        *b_BJet_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_BJet_flavor_weight_SoftMuonTag;   //!
   TBranch        *b_BJet_flavor_weight_SoftElectronTag;   //!
   TBranch        *b_BJet_flavor_truth_dRMinToB;   //!
   TBranch        *b_BJet_flavor_truth_dRMinToC;   //!
   TBranch        *b_BJet_flavor_truth_dRMinToT;   //!
   TBranch        *b_BJet_index;   //!
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventNumber;   //!
   TBranch        *b_lbn;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_Lepton1_inTrigger;   //!
   TBranch        *b_Lepton2_inTrigger;   //!
   TBranch        *b_isElectronTriggerPassed;   //!
   TBranch        *b_isMuonTriggerPassed;   //!
   TBranch        *b_VetoJet_n;   //!
   TBranch        *b_VetoJet_Barcode;   //!
   TBranch        *b_VetoJet_flavor_truth_trueflav;   //!
   TBranch        *b_VetoJet_flavor_truth_BHadronpdg;   //!
   TBranch        *b_VetoJet_pt;   //!
   TBranch        *b_VetoJet_eta;   //!
   TBranch        *b_VetoJet_phi;   //!
   TBranch        *b_VetoJet_E;   //!
   TBranch        *b_VetoJet_Et;   //!
   TBranch        *b_VetoJet_emscale_pt;   //!
   TBranch        *b_VetoJet_emscale_eta;   //!
   TBranch        *b_VetoJet_emscale_phi;   //!
   TBranch        *b_VetoJet_emscale_E;   //!
   TBranch        *b_VetoJet_flavor_weight;   //!
   TBranch        *b_VetoJet_flavor_weight_TrackCounting2D;   //!
   TBranch        *b_VetoJet_flavor_weight_JetProb;   //!
   TBranch        *b_VetoJet_flavor_weight_IP1D;   //!
   TBranch        *b_VetoJet_flavor_weight_IP2D;   //!
   TBranch        *b_VetoJet_flavor_weight_IP3D;   //!
   TBranch        *b_VetoJet_flavor_weight_SV0;   //!
   TBranch        *b_VetoJet_flavor_weight_SV1;   //!
   TBranch        *b_VetoJet_flavor_weight_SV2;   //!
   TBranch        *b_VetoJet_flavor_weight_JetFitterTag;   //!
   TBranch        *b_VetoJet_flavor_weight_JetFitterCOMB;   //!
   TBranch        *b_VetoJet_flavor_weight_JetFitterCOMBNN;   //!
   TBranch        *b_VetoJet_flavor_weight_SoftMuonTag;   //!
   TBranch        *b_VetoJet_flavor_weight_SoftElectronTag;   //!
   TBranch        *b_VetoJet_flavor_truth_dRMinToB;   //!
   TBranch        *b_VetoJet_flavor_truth_dRMinToC;   //!
   TBranch        *b_VetoJet_flavor_truth_dRMinToT;   //!
   TBranch        *b_isNoEleMuOverlap;   //!

   ElectronMuonPreTagTree(TTree *tree=0);
   virtual ~ElectronMuonPreTagTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ElectronMuonPreTagTree_cxx
ElectronMuonPreTagTree::ElectronMuonPreTagTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root");
      if (!f) {
         f = new TFile("/msu/data/t3fast2/single_top/TopPhysDPDMaker/MC10b/test/user.xiaohu.mc10_7TeV.108341.st_tchan_munu_McAtNlo_Jimmy.merge.AOD.e598_s933_s946_r2302_r2300_MC10b_EPS2011JUN09.110609121014/user.xiaohu.002350.NTUP_SGTOP._00001.root");
      }
      tree = (TTree*)gDirectory->Get("ElectronMuonPreTagTree");

   }
   Init(tree);
}

ElectronMuonPreTagTree::~ElectronMuonPreTagTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ElectronMuonPreTagTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ElectronMuonPreTagTree::LoadTree(Long64_t entry)
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

void ElectronMuonPreTagTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Jet_Barcode = 0;
   Jet_flavor_truth_trueflav = 0;
   Jet_flavor_truth_BHadronpdg = 0;
   Jet_pt = 0;
   Jet_eta = 0;
   Jet_phi = 0;
   Jet_E = 0;
   Jet_Et = 0;
   Jet_emscale_pt = 0;
   Jet_emscale_eta = 0;
   Jet_emscale_phi = 0;
   Jet_emscale_E = 0;
   Jet_flavor_weight = 0;
   Jet_flavor_weight_TrackCounting2D = 0;
   Jet_flavor_weight_JetProb = 0;
   Jet_flavor_weight_IP1D = 0;
   Jet_flavor_weight_IP2D = 0;
   Jet_flavor_weight_IP3D = 0;
   Jet_flavor_weight_SV0 = 0;
   Jet_flavor_weight_SV1 = 0;
   Jet_flavor_weight_SV2 = 0;
   Jet_flavor_weight_JetFitterTag = 0;
   Jet_flavor_weight_JetFitterCOMB = 0;
   Jet_flavor_weight_JetFitterCOMBNN = 0;
   Jet_flavor_weight_SoftMuonTag = 0;
   Jet_flavor_weight_SoftElectronTag = 0;
   Jet_flavor_truth_dRMinToB = 0;
   Jet_flavor_truth_dRMinToC = 0;
   Jet_flavor_truth_dRMinToT = 0;
   BJet_Barcode = 0;
   BJet_flavor_truth_trueflav = 0;
   BJet_flavor_truth_BHadronpdg = 0;
   BJet_pt = 0;
   BJet_eta = 0;
   BJet_phi = 0;
   BJet_E = 0;
   BJet_Et = 0;
   BJet_emscale_pt = 0;
   BJet_emscale_eta = 0;
   BJet_emscale_phi = 0;
   BJet_emscale_E = 0;
   BJet_flavor_weight = 0;
   BJet_flavor_weight_TrackCounting2D = 0;
   BJet_flavor_weight_JetProb = 0;
   BJet_flavor_weight_IP1D = 0;
   BJet_flavor_weight_IP2D = 0;
   BJet_flavor_weight_IP3D = 0;
   BJet_flavor_weight_SV0 = 0;
   BJet_flavor_weight_SV1 = 0;
   BJet_flavor_weight_SV2 = 0;
   BJet_flavor_weight_JetFitterTag = 0;
   BJet_flavor_weight_JetFitterCOMB = 0;
   BJet_flavor_weight_JetFitterCOMBNN = 0;
   BJet_flavor_weight_SoftMuonTag = 0;
   BJet_flavor_weight_SoftElectronTag = 0;
   BJet_flavor_truth_dRMinToB = 0;
   BJet_flavor_truth_dRMinToC = 0;
   BJet_flavor_truth_dRMinToT = 0;
   BJet_index = 0;
   VetoJet_Barcode = 0;
   VetoJet_flavor_truth_trueflav = 0;
   VetoJet_flavor_truth_BHadronpdg = 0;
   VetoJet_pt = 0;
   VetoJet_eta = 0;
   VetoJet_phi = 0;
   VetoJet_E = 0;
   VetoJet_Et = 0;
   VetoJet_emscale_pt = 0;
   VetoJet_emscale_eta = 0;
   VetoJet_emscale_phi = 0;
   VetoJet_emscale_E = 0;
   VetoJet_flavor_weight = 0;
   VetoJet_flavor_weight_TrackCounting2D = 0;
   VetoJet_flavor_weight_JetProb = 0;
   VetoJet_flavor_weight_IP1D = 0;
   VetoJet_flavor_weight_IP2D = 0;
   VetoJet_flavor_weight_IP3D = 0;
   VetoJet_flavor_weight_SV0 = 0;
   VetoJet_flavor_weight_SV1 = 0;
   VetoJet_flavor_weight_SV2 = 0;
   VetoJet_flavor_weight_JetFitterTag = 0;
   VetoJet_flavor_weight_JetFitterCOMB = 0;
   VetoJet_flavor_weight_JetFitterCOMBNN = 0;
   VetoJet_flavor_weight_SoftMuonTag = 0;
   VetoJet_flavor_weight_SoftElectronTag = 0;
   VetoJet_flavor_truth_dRMinToB = 0;
   VetoJet_flavor_truth_dRMinToC = 0;
   VetoJet_flavor_truth_dRMinToT = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Lepton1_loose", &Lepton1_loose, &b_Lepton1_loose);
   fChain->SetBranchAddress("Lepton1_medium", &Lepton1_medium, &b_Lepton1_medium);
   fChain->SetBranchAddress("Lepton1_mediumIso", &Lepton1_mediumIso, &b_Lepton1_mediumIso);
   fChain->SetBranchAddress("Lepton1_tight", &Lepton1_tight, &b_Lepton1_tight);
   fChain->SetBranchAddress("Lepton1_tightIso", &Lepton1_tightIso, &b_Lepton1_tightIso);
   fChain->SetBranchAddress("Lepton1_trackBlayer", &Lepton1_trackBlayer, &b_Lepton1_trackBlayer);
   fChain->SetBranchAddress("Lepton1_trackMatchEoverP", &Lepton1_trackMatchEoverP, &b_Lepton1_trackMatchEoverP);
   fChain->SetBranchAddress("Lepton1_trackTRThits", &Lepton1_trackTRThits, &b_Lepton1_trackTRThits);
   fChain->SetBranchAddress("Lepton1_trackTRTratio", &Lepton1_trackTRTratio, &b_Lepton1_trackTRTratio);
   fChain->SetBranchAddress("Lepton1_trackA0Tight", &Lepton1_trackA0Tight, &b_Lepton1_trackA0Tight);
   fChain->SetBranchAddress("Lepton1_trackMatchEtaTight", &Lepton1_trackMatchEtaTight, &b_Lepton1_trackMatchEtaTight);
   fChain->SetBranchAddress("Lepton1_trackMatchPhi", &Lepton1_trackMatchPhi, &b_Lepton1_trackMatchPhi);
   fChain->SetBranchAddress("Lepton1_conversionMatch", &Lepton1_conversionMatch, &b_Lepton1_conversionMatch);
   fChain->SetBranchAddress("Lepton1_Barcode", &Lepton1_Barcode, &b_Lepton1_Barcode);
   fChain->SetBranchAddress("Lepton1_pdgId", &Lepton1_pdgId, &b_Lepton1_pdgId);
   fChain->SetBranchAddress("Lepton1_charge", &Lepton1_charge, &b_Lepton1_charge);
   fChain->SetBranchAddress("Lepton1_Etcone45", &Lepton1_Etcone45, &b_Lepton1_Etcone45);
   fChain->SetBranchAddress("Lepton1_Etcone20", &Lepton1_Etcone20, &b_Lepton1_Etcone20);
   fChain->SetBranchAddress("Lepton1_Etcone30", &Lepton1_Etcone30, &b_Lepton1_Etcone30);
   fChain->SetBranchAddress("Lepton1_Etcone40", &Lepton1_Etcone40, &b_Lepton1_Etcone40);
   fChain->SetBranchAddress("Lepton1_ptcone20", &Lepton1_ptcone20, &b_Lepton1_ptcone20);
   fChain->SetBranchAddress("Lepton1_ptcone30", &Lepton1_ptcone30, &b_Lepton1_ptcone30);
   fChain->SetBranchAddress("Lepton1_ptcone40", &Lepton1_ptcone40, &b_Lepton1_ptcone40);
   fChain->SetBranchAddress("Lepton1_cluster_pt", &Lepton1_cluster_pt, &b_Lepton1_cluster_pt);
   fChain->SetBranchAddress("Lepton1_cluster_eta", &Lepton1_cluster_eta, &b_Lepton1_cluster_eta);
   fChain->SetBranchAddress("Lepton1_cluster_phi", &Lepton1_cluster_phi, &b_Lepton1_cluster_phi);
   fChain->SetBranchAddress("Lepton1_cluster_E", &Lepton1_cluster_E, &b_Lepton1_cluster_E);
   fChain->SetBranchAddress("Lepton1_cluster_Et", &Lepton1_cluster_Et, &b_Lepton1_cluster_Et);
   fChain->SetBranchAddress("Lepton1_cluster_track_Et", &Lepton1_cluster_track_Et, &b_Lepton1_cluster_track_Et);
   fChain->SetBranchAddress("Lepton1_track_d0", &Lepton1_track_d0, &b_Lepton1_track_d0);
   fChain->SetBranchAddress("Lepton1_track_z0", &Lepton1_track_z0, &b_Lepton1_track_z0);
   fChain->SetBranchAddress("Lepton1_track_phi", &Lepton1_track_phi, &b_Lepton1_track_phi);
   fChain->SetBranchAddress("Lepton1_track_theta", &Lepton1_track_theta, &b_Lepton1_track_theta);
   fChain->SetBranchAddress("Lepton1_track_qoverp", &Lepton1_track_qoverp, &b_Lepton1_track_qoverp);
   fChain->SetBranchAddress("Lepton1_track_pt", &Lepton1_track_pt, &b_Lepton1_track_pt);
   fChain->SetBranchAddress("Lepton1_track_eta", &Lepton1_track_eta, &b_Lepton1_track_eta);
   fChain->SetBranchAddress("Lepton2_Barcode", &Lepton2_Barcode, &b_Lepton2_Barcode);
   fChain->SetBranchAddress("Lepton2_pdgId", &Lepton2_pdgId, &b_Lepton2_pdgId);
   fChain->SetBranchAddress("Lepton2_isGoodIDTrack", &Lepton2_isGoodIDTrack, &b_Lepton2_isGoodIDTrack);
   fChain->SetBranchAddress("Lepton2_pt", &Lepton2_pt, &b_Lepton2_pt);
   fChain->SetBranchAddress("Lepton2_eta", &Lepton2_eta, &b_Lepton2_eta);
   fChain->SetBranchAddress("Lepton2_phi", &Lepton2_phi, &b_Lepton2_phi);
   fChain->SetBranchAddress("Lepton2_E", &Lepton2_E, &b_Lepton2_E);
   fChain->SetBranchAddress("Lepton2_Et", &Lepton2_Et, &b_Lepton2_Et);
   fChain->SetBranchAddress("Lepton2_charge", &Lepton2_charge, &b_Lepton2_charge);
   fChain->SetBranchAddress("Lepton2_etcone20", &Lepton2_etcone20, &b_Lepton2_etcone20);
   fChain->SetBranchAddress("Lepton2_etcone30", &Lepton2_etcone30, &b_Lepton2_etcone30);
   fChain->SetBranchAddress("Lepton2_etcone40", &Lepton2_etcone40, &b_Lepton2_etcone40);
   fChain->SetBranchAddress("Lepton2_ptcone20", &Lepton2_ptcone20, &b_Lepton2_ptcone20);
   fChain->SetBranchAddress("Lepton2_ptcone30", &Lepton2_ptcone30, &b_Lepton2_ptcone30);
   fChain->SetBranchAddress("Lepton2_ptcone40", &Lepton2_ptcone40, &b_Lepton2_ptcone40);
   fChain->SetBranchAddress("Lepton2_id_pt", &Lepton2_id_pt, &b_Lepton2_id_pt);
   fChain->SetBranchAddress("Lepton2_id_eta", &Lepton2_id_eta, &b_Lepton2_id_eta);
   fChain->SetBranchAddress("Lepton2_me_pt", &Lepton2_me_pt, &b_Lepton2_me_pt);
   fChain->SetBranchAddress("Lepton2_me_eta", &Lepton2_me_eta, &b_Lepton2_me_eta);
   fChain->SetBranchAddress("MissingEt_etx", &MissingEt_etx, &b_MissingEt_etx);
   fChain->SetBranchAddress("MissingEt_ety", &MissingEt_ety, &b_MissingEt_ety);
   fChain->SetBranchAddress("MissingEt_et", &MissingEt_et, &b_MissingEt_et);
   fChain->SetBranchAddress("MissingEt_sumet", &MissingEt_sumet, &b_MissingEt_sumet);
   fChain->SetBranchAddress("MissingEt_phi", &MissingEt_phi, &b_MissingEt_phi);
   fChain->SetBranchAddress("Jet_n", &Jet_n, &b_Jet_n);
   fChain->SetBranchAddress("Jet_Barcode", &Jet_Barcode, &b_Jet_Barcode);
   fChain->SetBranchAddress("Jet_flavor_truth_trueflav", &Jet_flavor_truth_trueflav, &b_Jet_flavor_truth_trueflav);
   fChain->SetBranchAddress("Jet_flavor_truth_BHadronpdg", &Jet_flavor_truth_BHadronpdg, &b_Jet_flavor_truth_BHadronpdg);
   fChain->SetBranchAddress("Jet_pt", &Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", &Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", &Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_E", &Jet_E, &b_Jet_E);
   fChain->SetBranchAddress("Jet_Et", &Jet_Et, &b_Jet_Et);
   fChain->SetBranchAddress("Jet_emscale_pt", &Jet_emscale_pt, &b_Jet_emscale_pt);
   fChain->SetBranchAddress("Jet_emscale_eta", &Jet_emscale_eta, &b_Jet_emscale_eta);
   fChain->SetBranchAddress("Jet_emscale_phi", &Jet_emscale_phi, &b_Jet_emscale_phi);
   fChain->SetBranchAddress("Jet_emscale_E", &Jet_emscale_E, &b_Jet_emscale_E);
   fChain->SetBranchAddress("Jet_flavor_weight", &Jet_flavor_weight, &b_Jet_flavor_weight);
   fChain->SetBranchAddress("Jet_flavor_weight_TrackCounting2D", &Jet_flavor_weight_TrackCounting2D, &b_Jet_flavor_weight_TrackCounting2D);
   fChain->SetBranchAddress("Jet_flavor_weight_JetProb", &Jet_flavor_weight_JetProb, &b_Jet_flavor_weight_JetProb);
   fChain->SetBranchAddress("Jet_flavor_weight_IP1D", &Jet_flavor_weight_IP1D, &b_Jet_flavor_weight_IP1D);
   fChain->SetBranchAddress("Jet_flavor_weight_IP2D", &Jet_flavor_weight_IP2D, &b_Jet_flavor_weight_IP2D);
   fChain->SetBranchAddress("Jet_flavor_weight_IP3D", &Jet_flavor_weight_IP3D, &b_Jet_flavor_weight_IP3D);
   fChain->SetBranchAddress("Jet_flavor_weight_SV0", &Jet_flavor_weight_SV0, &b_Jet_flavor_weight_SV0);
   fChain->SetBranchAddress("Jet_flavor_weight_SV1", &Jet_flavor_weight_SV1, &b_Jet_flavor_weight_SV1);
   fChain->SetBranchAddress("Jet_flavor_weight_SV2", &Jet_flavor_weight_SV2, &b_Jet_flavor_weight_SV2);
   fChain->SetBranchAddress("Jet_flavor_weight_JetFitterTag", &Jet_flavor_weight_JetFitterTag, &b_Jet_flavor_weight_JetFitterTag);
   fChain->SetBranchAddress("Jet_flavor_weight_JetFitterCOMB", &Jet_flavor_weight_JetFitterCOMB, &b_Jet_flavor_weight_JetFitterCOMB);
   fChain->SetBranchAddress("Jet_flavor_weight_JetFitterCOMBNN", &Jet_flavor_weight_JetFitterCOMBNN, &b_Jet_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("Jet_flavor_weight_SoftMuonTag", &Jet_flavor_weight_SoftMuonTag, &b_Jet_flavor_weight_SoftMuonTag);
   fChain->SetBranchAddress("Jet_flavor_weight_SoftElectronTag", &Jet_flavor_weight_SoftElectronTag, &b_Jet_flavor_weight_SoftElectronTag);
   fChain->SetBranchAddress("Jet_flavor_truth_dRMinToB", &Jet_flavor_truth_dRMinToB, &b_Jet_flavor_truth_dRMinToB);
   fChain->SetBranchAddress("Jet_flavor_truth_dRMinToC", &Jet_flavor_truth_dRMinToC, &b_Jet_flavor_truth_dRMinToC);
   fChain->SetBranchAddress("Jet_flavor_truth_dRMinToT", &Jet_flavor_truth_dRMinToT, &b_Jet_flavor_truth_dRMinToT);
   fChain->SetBranchAddress("DiLeptons_m", &DiLeptons_m, &b_DiLeptons_m);
   fChain->SetBranchAddress("DiLeptons_E", &DiLeptons_E, &b_DiLeptons_E);
   fChain->SetBranchAddress("DiLeptons_Et", &DiLeptons_Et, &b_DiLeptons_Et);
   fChain->SetBranchAddress("DiLeptons_pt", &DiLeptons_pt, &b_DiLeptons_pt);
   fChain->SetBranchAddress("DiLeptons_eta", &DiLeptons_eta, &b_DiLeptons_eta);
   fChain->SetBranchAddress("DiLeptons_phi", &DiLeptons_phi, &b_DiLeptons_phi);
   fChain->SetBranchAddress("DiLeptons_charge", &DiLeptons_charge, &b_DiLeptons_charge);
   fChain->SetBranchAddress("BJet_n", &BJet_n, &b_BJet_n);
   fChain->SetBranchAddress("BJet_Barcode", &BJet_Barcode, &b_BJet_Barcode);
   fChain->SetBranchAddress("BJet_flavor_truth_trueflav", &BJet_flavor_truth_trueflav, &b_BJet_flavor_truth_trueflav);
   fChain->SetBranchAddress("BJet_flavor_truth_BHadronpdg", &BJet_flavor_truth_BHadronpdg, &b_BJet_flavor_truth_BHadronpdg);
   fChain->SetBranchAddress("BJet_pt", &BJet_pt, &b_BJet_pt);
   fChain->SetBranchAddress("BJet_eta", &BJet_eta, &b_BJet_eta);
   fChain->SetBranchAddress("BJet_phi", &BJet_phi, &b_BJet_phi);
   fChain->SetBranchAddress("BJet_E", &BJet_E, &b_BJet_E);
   fChain->SetBranchAddress("BJet_Et", &BJet_Et, &b_BJet_Et);
   fChain->SetBranchAddress("BJet_emscale_pt", &BJet_emscale_pt, &b_BJet_emscale_pt);
   fChain->SetBranchAddress("BJet_emscale_eta", &BJet_emscale_eta, &b_BJet_emscale_eta);
   fChain->SetBranchAddress("BJet_emscale_phi", &BJet_emscale_phi, &b_BJet_emscale_phi);
   fChain->SetBranchAddress("BJet_emscale_E", &BJet_emscale_E, &b_BJet_emscale_E);
   fChain->SetBranchAddress("BJet_flavor_weight", &BJet_flavor_weight, &b_BJet_flavor_weight);
   fChain->SetBranchAddress("BJet_flavor_weight_TrackCounting2D", &BJet_flavor_weight_TrackCounting2D, &b_BJet_flavor_weight_TrackCounting2D);
   fChain->SetBranchAddress("BJet_flavor_weight_JetProb", &BJet_flavor_weight_JetProb, &b_BJet_flavor_weight_JetProb);
   fChain->SetBranchAddress("BJet_flavor_weight_IP1D", &BJet_flavor_weight_IP1D, &b_BJet_flavor_weight_IP1D);
   fChain->SetBranchAddress("BJet_flavor_weight_IP2D", &BJet_flavor_weight_IP2D, &b_BJet_flavor_weight_IP2D);
   fChain->SetBranchAddress("BJet_flavor_weight_IP3D", &BJet_flavor_weight_IP3D, &b_BJet_flavor_weight_IP3D);
   fChain->SetBranchAddress("BJet_flavor_weight_SV0", &BJet_flavor_weight_SV0, &b_BJet_flavor_weight_SV0);
   fChain->SetBranchAddress("BJet_flavor_weight_SV1", &BJet_flavor_weight_SV1, &b_BJet_flavor_weight_SV1);
   fChain->SetBranchAddress("BJet_flavor_weight_SV2", &BJet_flavor_weight_SV2, &b_BJet_flavor_weight_SV2);
   fChain->SetBranchAddress("BJet_flavor_weight_JetFitterTag", &BJet_flavor_weight_JetFitterTag, &b_BJet_flavor_weight_JetFitterTag);
   fChain->SetBranchAddress("BJet_flavor_weight_JetFitterCOMB", &BJet_flavor_weight_JetFitterCOMB, &b_BJet_flavor_weight_JetFitterCOMB);
   fChain->SetBranchAddress("BJet_flavor_weight_JetFitterCOMBNN", &BJet_flavor_weight_JetFitterCOMBNN, &b_BJet_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("BJet_flavor_weight_SoftMuonTag", &BJet_flavor_weight_SoftMuonTag, &b_BJet_flavor_weight_SoftMuonTag);
   fChain->SetBranchAddress("BJet_flavor_weight_SoftElectronTag", &BJet_flavor_weight_SoftElectronTag, &b_BJet_flavor_weight_SoftElectronTag);
   fChain->SetBranchAddress("BJet_flavor_truth_dRMinToB", &BJet_flavor_truth_dRMinToB, &b_BJet_flavor_truth_dRMinToB);
   fChain->SetBranchAddress("BJet_flavor_truth_dRMinToC", &BJet_flavor_truth_dRMinToC, &b_BJet_flavor_truth_dRMinToC);
   fChain->SetBranchAddress("BJet_flavor_truth_dRMinToT", &BJet_flavor_truth_dRMinToT, &b_BJet_flavor_truth_dRMinToT);
   fChain->SetBranchAddress("BJet_index", &BJet_index, &b_BJet_index);
   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
   fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("Lepton1_inTrigger", &Lepton1_inTrigger, &b_Lepton1_inTrigger);
   fChain->SetBranchAddress("Lepton2_inTrigger", &Lepton2_inTrigger, &b_Lepton2_inTrigger);
   fChain->SetBranchAddress("isElectronTriggerPassed", &isElectronTriggerPassed, &b_isElectronTriggerPassed);
   fChain->SetBranchAddress("isMuonTriggerPassed", &isMuonTriggerPassed, &b_isMuonTriggerPassed);
   fChain->SetBranchAddress("VetoJet_n", &VetoJet_n, &b_VetoJet_n);
   fChain->SetBranchAddress("VetoJet_Barcode", &VetoJet_Barcode, &b_VetoJet_Barcode);
   fChain->SetBranchAddress("VetoJet_flavor_truth_trueflav", &VetoJet_flavor_truth_trueflav, &b_VetoJet_flavor_truth_trueflav);
   fChain->SetBranchAddress("VetoJet_flavor_truth_BHadronpdg", &VetoJet_flavor_truth_BHadronpdg, &b_VetoJet_flavor_truth_BHadronpdg);
   fChain->SetBranchAddress("VetoJet_pt", &VetoJet_pt, &b_VetoJet_pt);
   fChain->SetBranchAddress("VetoJet_eta", &VetoJet_eta, &b_VetoJet_eta);
   fChain->SetBranchAddress("VetoJet_phi", &VetoJet_phi, &b_VetoJet_phi);
   fChain->SetBranchAddress("VetoJet_E", &VetoJet_E, &b_VetoJet_E);
   fChain->SetBranchAddress("VetoJet_Et", &VetoJet_Et, &b_VetoJet_Et);
   fChain->SetBranchAddress("VetoJet_emscale_pt", &VetoJet_emscale_pt, &b_VetoJet_emscale_pt);
   fChain->SetBranchAddress("VetoJet_emscale_eta", &VetoJet_emscale_eta, &b_VetoJet_emscale_eta);
   fChain->SetBranchAddress("VetoJet_emscale_phi", &VetoJet_emscale_phi, &b_VetoJet_emscale_phi);
   fChain->SetBranchAddress("VetoJet_emscale_E", &VetoJet_emscale_E, &b_VetoJet_emscale_E);
   fChain->SetBranchAddress("VetoJet_flavor_weight", &VetoJet_flavor_weight, &b_VetoJet_flavor_weight);
   fChain->SetBranchAddress("VetoJet_flavor_weight_TrackCounting2D", &VetoJet_flavor_weight_TrackCounting2D, &b_VetoJet_flavor_weight_TrackCounting2D);
   fChain->SetBranchAddress("VetoJet_flavor_weight_JetProb", &VetoJet_flavor_weight_JetProb, &b_VetoJet_flavor_weight_JetProb);
   fChain->SetBranchAddress("VetoJet_flavor_weight_IP1D", &VetoJet_flavor_weight_IP1D, &b_VetoJet_flavor_weight_IP1D);
   fChain->SetBranchAddress("VetoJet_flavor_weight_IP2D", &VetoJet_flavor_weight_IP2D, &b_VetoJet_flavor_weight_IP2D);
   fChain->SetBranchAddress("VetoJet_flavor_weight_IP3D", &VetoJet_flavor_weight_IP3D, &b_VetoJet_flavor_weight_IP3D);
   fChain->SetBranchAddress("VetoJet_flavor_weight_SV0", &VetoJet_flavor_weight_SV0, &b_VetoJet_flavor_weight_SV0);
   fChain->SetBranchAddress("VetoJet_flavor_weight_SV1", &VetoJet_flavor_weight_SV1, &b_VetoJet_flavor_weight_SV1);
   fChain->SetBranchAddress("VetoJet_flavor_weight_SV2", &VetoJet_flavor_weight_SV2, &b_VetoJet_flavor_weight_SV2);
   fChain->SetBranchAddress("VetoJet_flavor_weight_JetFitterTag", &VetoJet_flavor_weight_JetFitterTag, &b_VetoJet_flavor_weight_JetFitterTag);
   fChain->SetBranchAddress("VetoJet_flavor_weight_JetFitterCOMB", &VetoJet_flavor_weight_JetFitterCOMB, &b_VetoJet_flavor_weight_JetFitterCOMB);
   fChain->SetBranchAddress("VetoJet_flavor_weight_JetFitterCOMBNN", &VetoJet_flavor_weight_JetFitterCOMBNN, &b_VetoJet_flavor_weight_JetFitterCOMBNN);
   fChain->SetBranchAddress("VetoJet_flavor_weight_SoftMuonTag", &VetoJet_flavor_weight_SoftMuonTag, &b_VetoJet_flavor_weight_SoftMuonTag);
   fChain->SetBranchAddress("VetoJet_flavor_weight_SoftElectronTag", &VetoJet_flavor_weight_SoftElectronTag, &b_VetoJet_flavor_weight_SoftElectronTag);
   fChain->SetBranchAddress("VetoJet_flavor_truth_dRMinToB", &VetoJet_flavor_truth_dRMinToB, &b_VetoJet_flavor_truth_dRMinToB);
   fChain->SetBranchAddress("VetoJet_flavor_truth_dRMinToC", &VetoJet_flavor_truth_dRMinToC, &b_VetoJet_flavor_truth_dRMinToC);
   fChain->SetBranchAddress("VetoJet_flavor_truth_dRMinToT", &VetoJet_flavor_truth_dRMinToT, &b_VetoJet_flavor_truth_dRMinToT);
   fChain->SetBranchAddress("isNoEleMuOverlap", &isNoEleMuOverlap, &b_isNoEleMuOverlap);
   Notify();
}

Bool_t ElectronMuonPreTagTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ElectronMuonPreTagTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ElectronMuonPreTagTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ElectronMuonPreTagTree_cxx
