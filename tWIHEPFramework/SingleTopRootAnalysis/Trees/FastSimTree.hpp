//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 11 13:16:41 2007 by ROOT version 5.14/00e
// from TTree FastSim0/FastSim0
// found on file: /work/hx6/atlas/ntuples/single_top/csc/TopView1214/TopView1214_MuidTauRec_Correction/005502.AcerMC_tchanv12000601/user.csct8.TopView121402_MuidTauRec_Correction.trig1_misal1_mc12.005502.AcerMC_tchanv12000601.002.AANT0._00001.root
// Modified by RS: - bomb if input file is not found.
//                 - replace "vector" by "std::vector"
//                 - remove duplicate occurrences of several variables (fBits, N, etc)
//
//////////////////////////////////////////////////////////

#ifndef FastSimTree_h
#define FastSimTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <vector>
#include <iostream>

class FastSimTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leave types
   Double_t        AtlfastMissingEt_ex;
   Double_t        AtlfastMissingEt_ey;
   Double_t        AtlfastMissingEt_et;
   Double_t        AtlfastMissingEt_sumet;
   Double_t        MissingEx;
   Double_t        MissingEy;
   Double_t        MissingEt;
   Double_t        SumEt;
   Double_t        eventWeight;
   Double_t        eventWeightMCatNLO;
   Double_t        eventWeightIP2DWp3_1incl;
   Double_t        eventWeightIP2DWp3_1excl;
   Double_t        eventWeightIP2DWp3_2incl;
   Double_t        eventWeightIP2DWp3_2excl;
   Double_t        eventWeightSV1Wp3_1incl;
   Double_t        eventWeightSV1Wp3_1excl;
   Double_t        eventWeightSV1Wp3_2incl;
   Double_t        eventWeightSV1Wp3_2excl;
   Double_t        eventWeightIP2DWp3_1incl_Plus;
   Double_t        eventWeightIP2DWp3_1excl_Plus;
   Double_t        eventWeightIP2DWp3_2incl_Plus;
   Double_t        eventWeightIP2DWp3_2excl_Plus;
   Double_t        eventWeightSV1Wp3_1incl_Plus;
   Double_t        eventWeightSV1Wp3_1excl_Plus;
   Double_t        eventWeightSV1Wp3_2incl_Plus;
   Double_t        eventWeightSV1Wp3_2excl_Plus;
   Double_t        eventWeightIP2DWp3_1incl_Minus;
   Double_t        eventWeightIP2DWp3_1excl_Minus;
   Double_t        eventWeightIP2DWp3_2incl_Minus;
   Double_t        eventWeightIP2DWp3_2excl_Minus;
   Double_t        eventWeightSV1Wp3_1incl_Minus;
   Double_t        eventWeightSV1Wp3_1excl_Minus;
   Double_t        eventWeightSV1Wp3_2incl_Minus;
   Double_t        eventWeightSV1Wp3_2excl_Minus;
   Double_t        MTotal;
   Double_t        HT;
   Double_t        Event_Sphericity;
   Double_t        Event_Aplanarity;
   Double_t        Event_Planarity;
   Int_t           N;
   Int_t           PJet_N;
   Int_t           runNumber;
   Int_t           eventNumber;
   //Int_t           N;
   Int_t           GoodLepton_N;
   Int_t           GoodJet_N;
   Int_t           NonLooseJet_N;
   Int_t           El_N;
   Int_t           Mu_N;
   Int_t           Ph_N;
   Int_t           Tau_N;
   Int_t           VectSumAll_N;
   Int_t           ForwardJet_N;
   std::vector<double>  *PJet_E;
   std::vector<double>  *PJet_p_T;
   std::vector<double>  *PJet_phi;
   std::vector<double>  *PJet_eta;
   std::vector<double>  *PJet_px;
   std::vector<double>  *PJet_py;
   std::vector<double>  *PJet_pz;
   std::vector<double>  *PJet_m;
   std::vector<double>  *PJet_charge;
   std::vector<double>  *PJet_pdgId;
   std::vector<double>  *PJet_MCCalibFn;
   std::vector<double>  *PJet_MCCalib;
   std::vector<double>  *PJet_MCCalibError;
   std::vector<double>  *PJet_BHypoCalFactor;
   std::vector<double>  *PJet_LightHypoCalFactor;
   std::vector<double>  *PJet_TauHypoCalFactor;
   std::vector<double>  *PJet_Tau1P3PHypoCalFactor;
   std::vector<double>  *PJet_DRTruB;
   std::vector<double>  *PJet_DRTruC;
   std::vector<double>  *PJet_DRTruT;
   std::vector<double>  *El_E;
   std::vector<double>  *El_p_T;
   std::vector<double>  *El_phi;
   std::vector<double>  *El_eta;
   std::vector<double>  *El_px;
   std::vector<double>  *El_py;
   std::vector<double>  *El_pz;
   std::vector<double>  *El_m;
   std::vector<double>  *El_charge;
   std::vector<double>  *El_pdgId;
   std::vector<double>  *Mu_E;
   std::vector<double>  *Mu_p_T;
   std::vector<double>  *Mu_phi;
   std::vector<double>  *Mu_eta;
   std::vector<double>  *Mu_px;
   std::vector<double>  *Mu_py;
   std::vector<double>  *Mu_pz;
   std::vector<double>  *Mu_m;
   std::vector<double>  *Mu_charge;
   std::vector<double>  *Mu_pdgId;
   std::vector<double>  *Ph_E;
   std::vector<double>  *Ph_p_T;
   std::vector<double>  *Ph_phi;
   std::vector<double>  *Ph_eta;
   std::vector<double>  *Ph_px;
   std::vector<double>  *Ph_py;
   std::vector<double>  *Ph_pz;
   std::vector<double>  *Ph_m;
   std::vector<double>  *Ph_charge;
   std::vector<double>  *Ph_pdgId;
   std::vector<double>  *Tau_E;
   std::vector<double>  *Tau_p_T;
   std::vector<double>  *Tau_phi;
   std::vector<double>  *Tau_eta;
   std::vector<double>  *Tau_px;
   std::vector<double>  *Tau_py;
   std::vector<double>  *Tau_pz;
   std::vector<double>  *Tau_m;
   std::vector<double>  *Tau_charge;
   std::vector<double>  *Tau_pdgId;
   std::vector<double>  *PJet_TRFIP2DWp3;
   std::vector<double>  *PJet_TRFSV1Wp3;
   std::vector<double>  *PJet_TRFIP2DWp3_Plus;
   std::vector<double>  *PJet_TRFSV1Wp3_Plus;
   std::vector<double>  *PJet_TRFIP2DWp3_Minus;
   std::vector<double>  *PJet_TRFSV1Wp3_Minus;
   std::vector<double>  *El_TransM;
   std::vector<double>  *Mu_TransM;
   std::vector<double>  *VectSumAll_E;
   std::vector<double>  *VectSumAll_p_T;
   std::vector<double>  *VectSumAll_phi;
   std::vector<double>  *VectSumAll_eta;
   std::vector<double>  *VectSumAll_px;
   std::vector<double>  *VectSumAll_py;
   std::vector<double>  *VectSumAll_pz;
   std::vector<double>  *VectSumAll_m;
   std::vector<double>  *VectSumAll_charge;
   std::vector<double>  *VectSumAll_pdgId;
   std::vector<double>  *PJet_cosThetaStar;
   std::vector<double>  *PJet_Tru_DeltaR;
   std::vector<double>  *PJet_Tru_E;
   std::vector<double>  *PJet_Tru_p_T;
   std::vector<double>  *PJet_Tru_phi;
   std::vector<double>  *PJet_Tru_eta;
   std::vector<double>  *PJet_Tru_px;
   std::vector<double>  *PJet_Tru_py;
   std::vector<double>  *PJet_Tru_pz;
   std::vector<double>  *PJet_Tru_m;
   std::vector<double>  *PJet_Tru_charge;
   std::vector<double>  *PJet_Tru_pdgId;
   std::vector<double>  *PJet_Tru_etcone;
   std::vector<double>  *PJet_Tru_etcone10;
   std::vector<double>  *PJet_Tru_etcone20;
   std::vector<double>  *PJet_Tru_etcone30;
   std::vector<double>  *PJet_Tru_etcone40;
   std::vector<double>  *PJet_Tru_etcone50;
   std::vector<double>  *PJet_Tru_etcone60;
   std::vector<double>  *PJet_Tru_etcone70;
   std::vector<double>  *PJet_Tru_status;
   std::vector<double>  *PJet_TruJet_DeltaR;
   std::vector<double>  *PJet_TruJet_E;
   std::vector<double>  *PJet_TruJet_p_T;
   std::vector<double>  *PJet_TruJet_phi;
   std::vector<double>  *PJet_TruJet_eta;
   std::vector<double>  *PJet_TruJet_px;
   std::vector<double>  *PJet_TruJet_py;
   std::vector<double>  *PJet_TruJet_pz;
   std::vector<double>  *PJet_TruJet_m;
   std::vector<double>  *PJet_TruJet_charge;
   std::vector<double>  *PJet_TruJet_pdgId;
   std::vector<double>  *El_Tru_DeltaR;
   std::vector<double>  *El_Tru_E;
   std::vector<double>  *El_Tru_p_T;
   std::vector<double>  *El_Tru_phi;
   std::vector<double>  *El_Tru_eta;
   std::vector<double>  *El_Tru_px;
   std::vector<double>  *El_Tru_py;
   std::vector<double>  *El_Tru_pz;
   std::vector<double>  *El_Tru_m;
   std::vector<double>  *El_Tru_charge;
   std::vector<double>  *El_Tru_pdgId;
   std::vector<double>  *El_Tru_etcone;
   std::vector<double>  *El_Tru_etcone10;
   std::vector<double>  *El_Tru_etcone20;
   std::vector<double>  *El_Tru_etcone30;
   std::vector<double>  *El_Tru_etcone40;
   std::vector<double>  *El_Tru_etcone50;
   std::vector<double>  *El_Tru_etcone60;
   std::vector<double>  *El_Tru_etcone70;
   std::vector<double>  *El_Tru_status;
   std::vector<double>  *El_TruJet_DeltaR;
   std::vector<double>  *El_TruTau_DeltaR;
   std::vector<double>  *Ph_Tru_DeltaR;
   std::vector<double>  *Ph_Tru_E;
   std::vector<double>  *Ph_Tru_p_T;
   std::vector<double>  *Ph_Tru_phi;
   std::vector<double>  *Ph_Tru_eta;
   std::vector<double>  *Ph_Tru_px;
   std::vector<double>  *Ph_Tru_py;
   std::vector<double>  *Ph_Tru_pz;
   std::vector<double>  *Ph_Tru_m;
   std::vector<double>  *Ph_Tru_charge;
   std::vector<double>  *Ph_Tru_pdgId;
   std::vector<double>  *Ph_Tru_etcone;
   std::vector<double>  *Ph_Tru_etcone10;
   std::vector<double>  *Ph_Tru_etcone20;
   std::vector<double>  *Ph_Tru_etcone30;
   std::vector<double>  *Ph_Tru_etcone40;
   std::vector<double>  *Ph_Tru_etcone50;
   std::vector<double>  *Ph_Tru_etcone60;
   std::vector<double>  *Ph_Tru_etcone70;
   std::vector<double>  *Ph_Tru_status;
   std::vector<double>  *Mu_Tru_DeltaR;
   std::vector<double>  *Mu_Tru_E;
   std::vector<double>  *Mu_Tru_p_T;
   std::vector<double>  *Mu_Tru_phi;
   std::vector<double>  *Mu_Tru_eta;
   std::vector<double>  *Mu_Tru_px;
   std::vector<double>  *Mu_Tru_py;
   std::vector<double>  *Mu_Tru_pz;
   std::vector<double>  *Mu_Tru_m;
   std::vector<double>  *Mu_Tru_charge;
   std::vector<double>  *Mu_Tru_pdgId;
   std::vector<double>  *Mu_Tru_etcone;
   std::vector<double>  *Mu_Tru_etcone10;
   std::vector<double>  *Mu_Tru_etcone20;
   std::vector<double>  *Mu_Tru_etcone30;
   std::vector<double>  *Mu_Tru_etcone40;
   std::vector<double>  *Mu_Tru_etcone50;
   std::vector<double>  *Mu_Tru_etcone60;
   std::vector<double>  *Mu_Tru_etcone70;
   std::vector<double>  *Mu_Tru_status;
   std::vector<double>  *Tau_Tru_DeltaR;
   std::vector<double>  *Tau_Tru_E;
   std::vector<double>  *Tau_Tru_p_T;
   std::vector<double>  *Tau_Tru_phi;
   std::vector<double>  *Tau_Tru_eta;
   std::vector<double>  *Tau_Tru_px;
   std::vector<double>  *Tau_Tru_py;
   std::vector<double>  *Tau_Tru_pz;
   std::vector<double>  *Tau_Tru_m;
   std::vector<double>  *Tau_Tru_charge;
   std::vector<double>  *Tau_Tru_pdgId;
   std::vector<double>  *Tau_Tru_tauEta;
   std::vector<double>  *Tau_Tru_tauPhi;
   std::vector<double>  *Tau_Tru_tauE;
   std::vector<double>  *Tau_Tru_tauET;
   std::vector<double>  *Tau_Tru_tauPT;
   std::vector<double>  *Tau_Tru_tauPx;
   std::vector<double>  *Tau_Tru_tauPy;
   std::vector<double>  *Tau_Tru_tauPz;
   std::vector<double>  *Tau_Tru_stablesEAll;
   std::vector<double>  *Tau_Tru_stablesEAllCharged;
   std::vector<double>  *Tau_Tru_stablesEAllPion;
   std::vector<double>  *Tau_Tru_stablesEAllKaon;
   std::vector<double>  *Tau_Tru_stablesEAllOtherCh;
   std::vector<double>  *Tau_Tru_stablesEAllNeut;
   std::vector<double>  *Tau_Tru_stablesEAllPiZ;
   std::vector<double>  *Tau_Tru_stablesEAllGamma;
   std::vector<double>  *Tau_Tru_stablesEAllOtherNeut;
   std::vector<double>  *Tau_Tru_stablesETAll;
   std::vector<double>  *Tau_Tru_stablesETAllCharged;
   std::vector<double>  *Tau_Tru_stablesETAllPion;
   std::vector<double>  *Tau_Tru_stablesETAllKaon;
   std::vector<double>  *Tau_Tru_stablesETAllOtherCh;
   std::vector<double>  *Tau_Tru_stablesETAllNeut;
   std::vector<double>  *Tau_Tru_stablesETAllPiZ;
   std::vector<double>  *Tau_Tru_stablesETAllGamma;
   std::vector<double>  *Tau_Tru_stablesETAllOtherNeut;
   std::vector<double>  *Tau_Tru_stablesPTAll;
   std::vector<double>  *Tau_Tru_stablesPTAllCharged;
   std::vector<double>  *Tau_Tru_stablesPTAllPion;
   std::vector<double>  *Tau_Tru_stablesPTAllKaon;
   std::vector<double>  *Tau_Tru_stablesPTAllOtherCh;
   std::vector<double>  *Tau_Tru_stablesPTAllNeut;
   std::vector<double>  *Tau_Tru_stablesPTAllPiZ;
   std::vector<double>  *Tau_Tru_stablesPTAllGamma;
   std::vector<double>  *Tau_Tru_stablesPTAllOtherNeut;
   std::vector<int>     *PJet_FlatBTagged;
   std::vector<int>     *PJet_FlatBTaggedPlus;
   std::vector<int>     *PJet_FlatBTaggedMinus;
   std::vector<int>     *PJet_NonLooseJet;
   std::vector<int>     *PJet_GoodJet;
   std::vector<int>     *El_Loose;
   std::vector<int>     *El_GoodLepton;
   std::vector<int>     *Mu_Loose;
   std::vector<int>     *Mu_GoodLepton;
   std::vector<int>     *Mu_MuonInJet;
   std::vector<int>     *PJet_TRFIP2DWp3_1incl;
   std::vector<int>     *PJet_TRFIP2DWp3_1excl;
   std::vector<int>     *PJet_TRFIP2DWp3_2incl;
   std::vector<int>     *PJet_TRFIP2DWp3_2excl;
   std::vector<int>     *PJet_TRFSV1Wp3_1incl;
   std::vector<int>     *PJet_TRFSV1Wp3_1excl;
   std::vector<int>     *PJet_TRFSV1Wp3_2incl;
   std::vector<int>     *PJet_TRFSV1Wp3_2excl;
   std::vector<int>     *PJet_TRFIP2DWp3_Plus_1incl;
   std::vector<int>     *PJet_TRFIP2DWp3_Plus_1excl;
   std::vector<int>     *PJet_TRFIP2DWp3_Plus_2incl;
   std::vector<int>     *PJet_TRFIP2DWp3_Plus_2excl;
   std::vector<int>     *PJet_TRFSV1Wp3_Plus_1incl;
   std::vector<int>     *PJet_TRFSV1Wp3_Plus_1excl;
   std::vector<int>     *PJet_TRFSV1Wp3_Plus_2incl;
   std::vector<int>     *PJet_TRFSV1Wp3_Plus_2excl;
   std::vector<int>     *PJet_TRFIP2DWp3_Minus_1incl;
   std::vector<int>     *PJet_TRFIP2DWp3_Minus_1excl;
   std::vector<int>     *PJet_TRFIP2DWp3_Minus_2incl;
   std::vector<int>     *PJet_TRFIP2DWp3_Minus_2excl;
   std::vector<int>     *PJet_TRFSV1Wp3_Minus_1incl;
   std::vector<int>     *PJet_TRFSV1Wp3_Minus_1excl;
   std::vector<int>     *PJet_TRFSV1Wp3_Minus_2incl;
   std::vector<int>     *PJet_TRFSV1Wp3_Minus_2excl;
   std::vector<int>     *PJet_Tru_Matched;
   std::vector<int>     *PJet_Tru_barcode;
   std::vector<int>     *PJet_TruJet_Matched;
   std::vector<int>     *El_Tru_Matched;
   std::vector<int>     *El_Tru_barcode;
   std::vector<int>     *El_TruJet_Matched;
   std::vector<int>     *El_TruTau_Matched;
   std::vector<int>     *Ph_Tru_Matched;
   std::vector<int>     *Ph_Tru_barcode;
   std::vector<int>     *Mu_Tru_Matched;
   std::vector<int>     *Mu_Tru_barcode;
   std::vector<int>     *Tau_Tru_Matched;
   std::vector<int>     *Tau_Tru_compPartIndex;
   std::vector<int>     *Tau_Tru_test_FailHasCharge;
   std::vector<int>     *Tau_Tru_stablesN;
   std::vector<int>     *Tau_Tru_nDaughters;
   std::vector<int>     *Tau_Tru_nStable;
   std::vector<int>     *Tau_Tru_nStableCharged;
   std::vector<int>     *Tau_Tru_nStablePion;
   std::vector<int>     *Tau_Tru_nStableKaon;
   std::vector<int>     *Tau_Tru_nStableOtherCh;
   std::vector<int>     *Tau_Tru_nStableNeut;
   std::vector<int>     *Tau_Tru_nStablePiZ;
   std::vector<int>     *Tau_Tru_nStableGamma;
   std::vector<int>     *Tau_Tru_nStableOtherNeut;
   std::vector<int>     *Tau_Tru_nStableCounted;
   std::vector<int>     *Tau_Tru_nStableNoneMissing;
 //TStd::Vector3        *Event_ThirdAxis;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Double_t        fX;
   Double_t        fY;
   Double_t        fZ;
   Int_t           FastSimInstance;
   Int_t           FastSimCounter;
   Int_t           FastSimEventAndInstance;
   Int_t           FastSimRunNumber;
   Int_t           FastSimNInstance;

   // List of branches
   TBranch        *b_AtlfastMissingEt_ex;   //!
   TBranch        *b_AtlfastMissingEt_ey;   //!
   TBranch        *b_AtlfastMissingEt_et;   //!
   TBranch        *b_AtlfastMissingEt_sumet;   //!
   TBranch        *b_MissingEx;   //!
   TBranch        *b_MissingEy;   //!
   TBranch        *b_MissingEt;   //!
   TBranch        *b_SumEt;   //!
   TBranch        *b_eventWeight;   //!
   TBranch        *b_eventWeightMCatNLO;   //!
   TBranch        *b_eventWeightIP2DWp3_1incl;   //!
   TBranch        *b_eventWeightIP2DWp3_1excl;   //!
   TBranch        *b_eventWeightIP2DWp3_2incl;   //!
   TBranch        *b_eventWeightIP2DWp3_2excl;   //!
   TBranch        *b_eventWeightSV1Wp3_1incl;   //!
   TBranch        *b_eventWeightSV1Wp3_1excl;   //!
   TBranch        *b_eventWeightSV1Wp3_2incl;   //!
   TBranch        *b_eventWeightSV1Wp3_2excl;   //!
   TBranch        *b_eventWeightIP2DWp3_1incl_Plus;   //!
   TBranch        *b_eventWeightIP2DWp3_1excl_Plus;   //!
   TBranch        *b_eventWeightIP2DWp3_2incl_Plus;   //!
   TBranch        *b_eventWeightIP2DWp3_2excl_Plus;   //!
   TBranch        *b_eventWeightSV1Wp3_1incl_Plus;   //!
   TBranch        *b_eventWeightSV1Wp3_1excl_Plus;   //!
   TBranch        *b_eventWeightSV1Wp3_2incl_Plus;   //!
   TBranch        *b_eventWeightSV1Wp3_2excl_Plus;   //!
   TBranch        *b_eventWeightIP2DWp3_1incl_Minus;   //!
   TBranch        *b_eventWeightIP2DWp3_1excl_Minus;   //!
   TBranch        *b_eventWeightIP2DWp3_2incl_Minus;   //!
   TBranch        *b_eventWeightIP2DWp3_2excl_Minus;   //!
   TBranch        *b_eventWeightSV1Wp3_1incl_Minus;   //!
   TBranch        *b_eventWeightSV1Wp3_1excl_Minus;   //!
   TBranch        *b_eventWeightSV1Wp3_2incl_Minus;   //!
   TBranch        *b_eventWeightSV1Wp3_2excl_Minus;   //!
   TBranch        *b_MTotal;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_Event_Sphericity;   //!
   TBranch        *b_Event_Aplanarity;   //!
   TBranch        *b_Event_Planarity;   //!
   TBranch        *b_N;   //!
   TBranch        *b_PJet_N;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   //TBranch        *b_N;   //!
   TBranch        *b_GoodLepton_N;   //!
   TBranch        *b_GoodJet_N;   //!
   TBranch        *b_NonLooseJet_N;   //!
   TBranch        *b_El_N;   //!
   TBranch        *b_Mu_N;   //!
   TBranch        *b_Ph_N;   //!
   TBranch        *b_Tau_N;   //!
   TBranch        *b_VectSumAll_N;   //!
   TBranch        *b_ForwardJet_N;   //!
   TBranch        *b_PJet_E;   //!
   TBranch        *b_PJet_p_T;   //!
   TBranch        *b_PJet_phi;   //!
   TBranch        *b_PJet_eta;   //!
   TBranch        *b_PJet_px;   //!
   TBranch        *b_PJet_py;   //!
   TBranch        *b_PJet_pz;   //!
   TBranch        *b_PJet_m;   //!
   TBranch        *b_PJet_charge;   //!
   TBranch        *b_PJet_pdgId;   //!
   TBranch        *b_PJet_MCCalibFn;   //!
   TBranch        *b_PJet_MCCalib;   //!
   TBranch        *b_PJet_MCCalibError;   //!
   TBranch        *b_PJet_BHypoCalFactor;   //!
   TBranch        *b_PJet_LightHypoCalFactor;   //!
   TBranch        *b_PJet_TauHypoCalFactor;   //!
   TBranch        *b_PJet_Tau1P3PHypoCalFactor;   //!
   TBranch        *b_PJet_DRTruB;   //!
   TBranch        *b_PJet_DRTruC;   //!
   TBranch        *b_PJet_DRTruT;   //!
   TBranch        *b_El_E;   //!
   TBranch        *b_El_p_T;   //!
   TBranch        *b_El_phi;   //!
   TBranch        *b_El_eta;   //!
   TBranch        *b_El_px;   //!
   TBranch        *b_El_py;   //!
   TBranch        *b_El_pz;   //!
   TBranch        *b_El_m;   //!
   TBranch        *b_El_charge;   //!
   TBranch        *b_El_pdgId;   //!
   TBranch        *b_Mu_E;   //!
   TBranch        *b_Mu_p_T;   //!
   TBranch        *b_Mu_phi;   //!
   TBranch        *b_Mu_eta;   //!
   TBranch        *b_Mu_px;   //!
   TBranch        *b_Mu_py;   //!
   TBranch        *b_Mu_pz;   //!
   TBranch        *b_Mu_m;   //!
   TBranch        *b_Mu_charge;   //!
   TBranch        *b_Mu_pdgId;   //!
   TBranch        *b_Ph_E;   //!
   TBranch        *b_Ph_p_T;   //!
   TBranch        *b_Ph_phi;   //!
   TBranch        *b_Ph_eta;   //!
   TBranch        *b_Ph_px;   //!
   TBranch        *b_Ph_py;   //!
   TBranch        *b_Ph_pz;   //!
   TBranch        *b_Ph_m;   //!
   TBranch        *b_Ph_charge;   //!
   TBranch        *b_Ph_pdgId;   //!
   TBranch        *b_Tau_E;   //!
   TBranch        *b_Tau_p_T;   //!
   TBranch        *b_Tau_phi;   //!
   TBranch        *b_Tau_eta;   //!
   TBranch        *b_Tau_px;   //!
   TBranch        *b_Tau_py;   //!
   TBranch        *b_Tau_pz;   //!
   TBranch        *b_Tau_m;   //!
   TBranch        *b_Tau_charge;   //!
   TBranch        *b_Tau_pdgId;   //!
   TBranch        *b_PJet_TRFIP2DWp3;   //!
   TBranch        *b_PJet_TRFSV1Wp3;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Plus;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Plus;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Minus;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Minus;   //!
   TBranch        *b_El_TransM;   //!
   TBranch        *b_Mu_TransM;   //!
   TBranch        *b_VectSumAll_E;   //!
   TBranch        *b_VectSumAll_p_T;   //!
   TBranch        *b_VectSumAll_phi;   //!
   TBranch        *b_VectSumAll_eta;   //!
   TBranch        *b_VectSumAll_px;   //!
   TBranch        *b_VectSumAll_py;   //!
   TBranch        *b_VectSumAll_pz;   //!
   TBranch        *b_VectSumAll_m;   //!
   TBranch        *b_VectSumAll_charge;   //!
   TBranch        *b_VectSumAll_pdgId;   //!
   TBranch        *b_PJet_cosThetaStar;   //!
   TBranch        *b_PJet_Tru_DeltaR;   //!
   TBranch        *b_PJet_Tru_E;   //!
   TBranch        *b_PJet_Tru_p_T;   //!
   TBranch        *b_PJet_Tru_phi;   //!
   TBranch        *b_PJet_Tru_eta;   //!
   TBranch        *b_PJet_Tru_px;   //!
   TBranch        *b_PJet_Tru_py;   //!
   TBranch        *b_PJet_Tru_pz;   //!
   TBranch        *b_PJet_Tru_m;   //!
   TBranch        *b_PJet_Tru_charge;   //!
   TBranch        *b_PJet_Tru_pdgId;   //!
   TBranch        *b_PJet_Tru_etcone;   //!
   TBranch        *b_PJet_Tru_etcone10;   //!
   TBranch        *b_PJet_Tru_etcone20;   //!
   TBranch        *b_PJet_Tru_etcone30;   //!
   TBranch        *b_PJet_Tru_etcone40;   //!
   TBranch        *b_PJet_Tru_etcone50;   //!
   TBranch        *b_PJet_Tru_etcone60;   //!
   TBranch        *b_PJet_Tru_etcone70;   //!
   TBranch        *b_PJet_Tru_status;   //!
   TBranch        *b_PJet_TruJet_DeltaR;   //!
   TBranch        *b_PJet_TruJet_E;   //!
   TBranch        *b_PJet_TruJet_p_T;   //!
   TBranch        *b_PJet_TruJet_phi;   //!
   TBranch        *b_PJet_TruJet_eta;   //!
   TBranch        *b_PJet_TruJet_px;   //!
   TBranch        *b_PJet_TruJet_py;   //!
   TBranch        *b_PJet_TruJet_pz;   //!
   TBranch        *b_PJet_TruJet_m;   //!
   TBranch        *b_PJet_TruJet_charge;   //!
   TBranch        *b_PJet_TruJet_pdgId;   //!
   TBranch        *b_El_Tru_DeltaR;   //!
   TBranch        *b_El_Tru_E;   //!
   TBranch        *b_El_Tru_p_T;   //!
   TBranch        *b_El_Tru_phi;   //!
   TBranch        *b_El_Tru_eta;   //!
   TBranch        *b_El_Tru_px;   //!
   TBranch        *b_El_Tru_py;   //!
   TBranch        *b_El_Tru_pz;   //!
   TBranch        *b_El_Tru_m;   //!
   TBranch        *b_El_Tru_charge;   //!
   TBranch        *b_El_Tru_pdgId;   //!
   TBranch        *b_El_Tru_etcone;   //!
   TBranch        *b_El_Tru_etcone10;   //!
   TBranch        *b_El_Tru_etcone20;   //!
   TBranch        *b_El_Tru_etcone30;   //!
   TBranch        *b_El_Tru_etcone40;   //!
   TBranch        *b_El_Tru_etcone50;   //!
   TBranch        *b_El_Tru_etcone60;   //!
   TBranch        *b_El_Tru_etcone70;   //!
   TBranch        *b_El_Tru_status;   //!
   TBranch        *b_El_TruJet_DeltaR;   //!
   TBranch        *b_El_TruTau_DeltaR;   //!
   TBranch        *b_Ph_Tru_DeltaR;   //!
   TBranch        *b_Ph_Tru_E;   //!
   TBranch        *b_Ph_Tru_p_T;   //!
   TBranch        *b_Ph_Tru_phi;   //!
   TBranch        *b_Ph_Tru_eta;   //!
   TBranch        *b_Ph_Tru_px;   //!
   TBranch        *b_Ph_Tru_py;   //!
   TBranch        *b_Ph_Tru_pz;   //!
   TBranch        *b_Ph_Tru_m;   //!
   TBranch        *b_Ph_Tru_charge;   //!
   TBranch        *b_Ph_Tru_pdgId;   //!
   TBranch        *b_Ph_Tru_etcone;   //!
   TBranch        *b_Ph_Tru_etcone10;   //!
   TBranch        *b_Ph_Tru_etcone20;   //!
   TBranch        *b_Ph_Tru_etcone30;   //!
   TBranch        *b_Ph_Tru_etcone40;   //!
   TBranch        *b_Ph_Tru_etcone50;   //!
   TBranch        *b_Ph_Tru_etcone60;   //!
   TBranch        *b_Ph_Tru_etcone70;   //!
   TBranch        *b_Ph_Tru_status;   //!
   TBranch        *b_Mu_Tru_DeltaR;   //!
   TBranch        *b_Mu_Tru_E;   //!
   TBranch        *b_Mu_Tru_p_T;   //!
   TBranch        *b_Mu_Tru_phi;   //!
   TBranch        *b_Mu_Tru_eta;   //!
   TBranch        *b_Mu_Tru_px;   //!
   TBranch        *b_Mu_Tru_py;   //!
   TBranch        *b_Mu_Tru_pz;   //!
   TBranch        *b_Mu_Tru_m;   //!
   TBranch        *b_Mu_Tru_charge;   //!
   TBranch        *b_Mu_Tru_pdgId;   //!
   TBranch        *b_Mu_Tru_etcone;   //!
   TBranch        *b_Mu_Tru_etcone10;   //!
   TBranch        *b_Mu_Tru_etcone20;   //!
   TBranch        *b_Mu_Tru_etcone30;   //!
   TBranch        *b_Mu_Tru_etcone40;   //!
   TBranch        *b_Mu_Tru_etcone50;   //!
   TBranch        *b_Mu_Tru_etcone60;   //!
   TBranch        *b_Mu_Tru_etcone70;   //!
   TBranch        *b_Mu_Tru_status;   //!
   TBranch        *b_Tau_Tru_DeltaR;   //!
   TBranch        *b_Tau_Tru_E;   //!
   TBranch        *b_Tau_Tru_p_T;   //!
   TBranch        *b_Tau_Tru_phi;   //!
   TBranch        *b_Tau_Tru_eta;   //!
   TBranch        *b_Tau_Tru_px;   //!
   TBranch        *b_Tau_Tru_py;   //!
   TBranch        *b_Tau_Tru_pz;   //!
   TBranch        *b_Tau_Tru_m;   //!
   TBranch        *b_Tau_Tru_charge;   //!
   TBranch        *b_Tau_Tru_pdgId;   //!
   TBranch        *b_Tau_Tru_tauEta;   //!
   TBranch        *b_Tau_Tru_tauPhi;   //!
   TBranch        *b_Tau_Tru_tauE;   //!
   TBranch        *b_Tau_Tru_tauET;   //!
   TBranch        *b_Tau_Tru_tauPT;   //!
   TBranch        *b_Tau_Tru_tauPx;   //!
   TBranch        *b_Tau_Tru_tauPy;   //!
   TBranch        *b_Tau_Tru_tauPz;   //!
   TBranch        *b_Tau_Tru_stablesEAll;   //!
   TBranch        *b_Tau_Tru_stablesEAllCharged;   //!
   TBranch        *b_Tau_Tru_stablesEAllPion;   //!
   TBranch        *b_Tau_Tru_stablesEAllKaon;   //!
   TBranch        *b_Tau_Tru_stablesEAllOtherCh;   //!
   TBranch        *b_Tau_Tru_stablesEAllNeut;   //!
   TBranch        *b_Tau_Tru_stablesEAllPiZ;   //!
   TBranch        *b_Tau_Tru_stablesEAllGamma;   //!
   TBranch        *b_Tau_Tru_stablesEAllOtherNeut;   //!
   TBranch        *b_Tau_Tru_stablesETAll;   //!
   TBranch        *b_Tau_Tru_stablesETAllCharged;   //!
   TBranch        *b_Tau_Tru_stablesETAllPion;   //!
   TBranch        *b_Tau_Tru_stablesETAllKaon;   //!
   TBranch        *b_Tau_Tru_stablesETAllOtherCh;   //!
   TBranch        *b_Tau_Tru_stablesETAllNeut;   //!
   TBranch        *b_Tau_Tru_stablesETAllPiZ;   //!
   TBranch        *b_Tau_Tru_stablesETAllGamma;   //!
   TBranch        *b_Tau_Tru_stablesETAllOtherNeut;   //!
   TBranch        *b_Tau_Tru_stablesPTAll;   //!
   TBranch        *b_Tau_Tru_stablesPTAllCharged;   //!
   TBranch        *b_Tau_Tru_stablesPTAllPion;   //!
   TBranch        *b_Tau_Tru_stablesPTAllKaon;   //!
   TBranch        *b_Tau_Tru_stablesPTAllOtherCh;   //!
   TBranch        *b_Tau_Tru_stablesPTAllNeut;   //!
   TBranch        *b_Tau_Tru_stablesPTAllPiZ;   //!
   TBranch        *b_Tau_Tru_stablesPTAllGamma;   //!
   TBranch        *b_Tau_Tru_stablesPTAllOtherNeut;   //!
   TBranch        *b_PJet_FlatBTagged;   //!
   TBranch        *b_PJet_FlatBTaggedPlus;   //!
   TBranch        *b_PJet_FlatBTaggedMinus;   //!
   TBranch        *b_PJet_NonLooseJet;   //!
   TBranch        *b_PJet_GoodJet;   //!
   TBranch        *b_El_Loose;   //!
   TBranch        *b_El_GoodLepton;   //!
   TBranch        *b_Mu_Loose;   //!
   TBranch        *b_Mu_GoodLepton;   //!
   TBranch        *b_Mu_MuonInJet;   //!
   TBranch        *b_PJet_TRFIP2DWp3_1incl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_1excl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_2incl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_2excl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_1incl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_1excl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_2incl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_2excl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Plus_1incl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Plus_1excl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Plus_2incl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Plus_2excl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Plus_1incl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Plus_1excl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Plus_2incl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Plus_2excl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Minus_1incl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Minus_1excl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Minus_2incl;   //!
   TBranch        *b_PJet_TRFIP2DWp3_Minus_2excl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Minus_1incl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Minus_1excl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Minus_2incl;   //!
   TBranch        *b_PJet_TRFSV1Wp3_Minus_2excl;   //!
   TBranch        *b_PJet_Tru_Matched;   //!
   TBranch        *b_PJet_Tru_barcode;   //!
   TBranch        *b_PJet_TruJet_Matched;   //!
   TBranch        *b_El_Tru_Matched;   //!
   TBranch        *b_El_Tru_barcode;   //!
   TBranch        *b_El_TruJet_Matched;   //!
   TBranch        *b_El_TruTau_Matched;   //!
   TBranch        *b_Ph_Tru_Matched;   //!
   TBranch        *b_Ph_Tru_barcode;   //!
   TBranch        *b_Mu_Tru_Matched;   //!
   TBranch        *b_Mu_Tru_barcode;   //!
   TBranch        *b_Tau_Tru_Matched;   //!
   TBranch        *b_Tau_Tru_compPartIndex;   //!
   TBranch        *b_Tau_Tru_test_FailHasCharge;   //!
   TBranch        *b_Tau_Tru_stablesN;   //!
   TBranch        *b_Tau_Tru_nDaughters;   //!
   TBranch        *b_Tau_Tru_nStable;   //!
   TBranch        *b_Tau_Tru_nStableCharged;   //!
   TBranch        *b_Tau_Tru_nStablePion;   //!
   TBranch        *b_Tau_Tru_nStableKaon;   //!
   TBranch        *b_Tau_Tru_nStableOtherCh;   //!
   TBranch        *b_Tau_Tru_nStableNeut;   //!
   TBranch        *b_Tau_Tru_nStablePiZ;   //!
   TBranch        *b_Tau_Tru_nStableGamma;   //!
   TBranch        *b_Tau_Tru_nStableOtherNeut;   //!
   TBranch        *b_Tau_Tru_nStableCounted;   //!
   TBranch        *b_Tau_Tru_nStableNoneMissing;   //!
   TBranch        *b_Event_SphericityAxis_fUniqueID;   //!
   TBranch        *b_Event_SphericityAxis_fBits;   //!
   TBranch        *b_Event_SphericityAxis_fX;   //!
   TBranch        *b_Event_SphericityAxis_fY;   //!
   TBranch        *b_Event_SphericityAxis_fZ;   //!
   TBranch        *b_Event_SecondAxis_fUniqueID;   //!
   TBranch        *b_Event_SecondAxis_fBits;   //!
   TBranch        *b_Event_SecondAxis_fX;   //!
   TBranch        *b_Event_SecondAxis_fY;   //!
   TBranch        *b_Event_SecondAxis_fZ;   //!
   TBranch        *b_Event_ThirdAxis_fUniqueID;   //!
   TBranch        *b_Event_ThirdAxis_fBits;   //!
   TBranch        *b_Event_ThirdAxis_fX;   //!
   TBranch        *b_Event_ThirdAxis_fY;   //!
   TBranch        *b_Event_ThirdAxis_fZ;   //!
   TBranch        *b_FastSimInstance;   //!
   TBranch        *b_FastSimCounter;   //!
   TBranch        *b_FastSimEventAndInstance;   //!
   TBranch        *b_FastSimRunNumber;   //!
   TBranch        *b_FastSimNInstance;   //!

   FastSimTree(TTree *tree=0);
   virtual ~FastSimTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef FastSimTree_cxx
FastSimTree::FastSimTree(TTree *tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
     std::cout<<"Error in FastSim tree constructor: no valid input tree!"<<std::endl;
     return;
   }
   Init(tree);
}

FastSimTree::~FastSimTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t FastSimTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t FastSimTree::LoadTree(Long64_t entry)
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

void FastSimTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normaly not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PJet_E = 0;
   PJet_p_T = 0;
   PJet_phi = 0;
   PJet_eta = 0;
   PJet_px = 0;
   PJet_py = 0;
   PJet_pz = 0;
   PJet_m = 0;
   PJet_charge = 0;
   PJet_pdgId = 0;
   PJet_MCCalibFn = 0;
   PJet_MCCalib = 0;
   PJet_MCCalibError = 0;
   PJet_BHypoCalFactor = 0;
   PJet_LightHypoCalFactor = 0;
   PJet_TauHypoCalFactor = 0;
   PJet_Tau1P3PHypoCalFactor = 0;
   PJet_DRTruB = 0;
   PJet_DRTruC = 0;
   PJet_DRTruT = 0;
   El_E = 0;
   El_p_T = 0;
   El_phi = 0;
   El_eta = 0;
   El_px = 0;
   El_py = 0;
   El_pz = 0;
   El_m = 0;
   El_charge = 0;
   El_pdgId = 0;
   Mu_E = 0;
   Mu_p_T = 0;
   Mu_phi = 0;
   Mu_eta = 0;
   Mu_px = 0;
   Mu_py = 0;
   Mu_pz = 0;
   Mu_m = 0;
   Mu_charge = 0;
   Mu_pdgId = 0;
   Ph_E = 0;
   Ph_p_T = 0;
   Ph_phi = 0;
   Ph_eta = 0;
   Ph_px = 0;
   Ph_py = 0;
   Ph_pz = 0;
   Ph_m = 0;
   Ph_charge = 0;
   Ph_pdgId = 0;
   Tau_E = 0;
   Tau_p_T = 0;
   Tau_phi = 0;
   Tau_eta = 0;
   Tau_px = 0;
   Tau_py = 0;
   Tau_pz = 0;
   Tau_m = 0;
   Tau_charge = 0;
   Tau_pdgId = 0;
   PJet_TRFIP2DWp3 = 0;
   PJet_TRFSV1Wp3 = 0;
   PJet_TRFIP2DWp3_Plus = 0;
   PJet_TRFSV1Wp3_Plus = 0;
   PJet_TRFIP2DWp3_Minus = 0;
   PJet_TRFSV1Wp3_Minus = 0;
   El_TransM = 0;
   Mu_TransM = 0;
   VectSumAll_E = 0;
   VectSumAll_p_T = 0;
   VectSumAll_phi = 0;
   VectSumAll_eta = 0;
   VectSumAll_px = 0;
   VectSumAll_py = 0;
   VectSumAll_pz = 0;
   VectSumAll_m = 0;
   VectSumAll_charge = 0;
   VectSumAll_pdgId = 0;
   PJet_cosThetaStar = 0;
   PJet_Tru_DeltaR = 0;
   PJet_Tru_E = 0;
   PJet_Tru_p_T = 0;
   PJet_Tru_phi = 0;
   PJet_Tru_eta = 0;
   PJet_Tru_px = 0;
   PJet_Tru_py = 0;
   PJet_Tru_pz = 0;
   PJet_Tru_m = 0;
   PJet_Tru_charge = 0;
   PJet_Tru_pdgId = 0;
   PJet_Tru_etcone = 0;
   PJet_Tru_etcone10 = 0;
   PJet_Tru_etcone20 = 0;
   PJet_Tru_etcone30 = 0;
   PJet_Tru_etcone40 = 0;
   PJet_Tru_etcone50 = 0;
   PJet_Tru_etcone60 = 0;
   PJet_Tru_etcone70 = 0;
   PJet_Tru_status = 0;
   PJet_TruJet_DeltaR = 0;
   PJet_TruJet_E = 0;
   PJet_TruJet_p_T = 0;
   PJet_TruJet_phi = 0;
   PJet_TruJet_eta = 0;
   PJet_TruJet_px = 0;
   PJet_TruJet_py = 0;
   PJet_TruJet_pz = 0;
   PJet_TruJet_m = 0;
   PJet_TruJet_charge = 0;
   PJet_TruJet_pdgId = 0;
   El_Tru_DeltaR = 0;
   El_Tru_E = 0;
   El_Tru_p_T = 0;
   El_Tru_phi = 0;
   El_Tru_eta = 0;
   El_Tru_px = 0;
   El_Tru_py = 0;
   El_Tru_pz = 0;
   El_Tru_m = 0;
   El_Tru_charge = 0;
   El_Tru_pdgId = 0;
   El_Tru_etcone = 0;
   El_Tru_etcone10 = 0;
   El_Tru_etcone20 = 0;
   El_Tru_etcone30 = 0;
   El_Tru_etcone40 = 0;
   El_Tru_etcone50 = 0;
   El_Tru_etcone60 = 0;
   El_Tru_etcone70 = 0;
   El_Tru_status = 0;
   El_TruJet_DeltaR = 0;
   El_TruTau_DeltaR = 0;
   Ph_Tru_DeltaR = 0;
   Ph_Tru_E = 0;
   Ph_Tru_p_T = 0;
   Ph_Tru_phi = 0;
   Ph_Tru_eta = 0;
   Ph_Tru_px = 0;
   Ph_Tru_py = 0;
   Ph_Tru_pz = 0;
   Ph_Tru_m = 0;
   Ph_Tru_charge = 0;
   Ph_Tru_pdgId = 0;
   Ph_Tru_etcone = 0;
   Ph_Tru_etcone10 = 0;
   Ph_Tru_etcone20 = 0;
   Ph_Tru_etcone30 = 0;
   Ph_Tru_etcone40 = 0;
   Ph_Tru_etcone50 = 0;
   Ph_Tru_etcone60 = 0;
   Ph_Tru_etcone70 = 0;
   Ph_Tru_status = 0;
   Mu_Tru_DeltaR = 0;
   Mu_Tru_E = 0;
   Mu_Tru_p_T = 0;
   Mu_Tru_phi = 0;
   Mu_Tru_eta = 0;
   Mu_Tru_px = 0;
   Mu_Tru_py = 0;
   Mu_Tru_pz = 0;
   Mu_Tru_m = 0;
   Mu_Tru_charge = 0;
   Mu_Tru_pdgId = 0;
   Mu_Tru_etcone = 0;
   Mu_Tru_etcone10 = 0;
   Mu_Tru_etcone20 = 0;
   Mu_Tru_etcone30 = 0;
   Mu_Tru_etcone40 = 0;
   Mu_Tru_etcone50 = 0;
   Mu_Tru_etcone60 = 0;
   Mu_Tru_etcone70 = 0;
   Mu_Tru_status = 0;
   Tau_Tru_DeltaR = 0;
   Tau_Tru_E = 0;
   Tau_Tru_p_T = 0;
   Tau_Tru_phi = 0;
   Tau_Tru_eta = 0;
   Tau_Tru_px = 0;
   Tau_Tru_py = 0;
   Tau_Tru_pz = 0;
   Tau_Tru_m = 0;
   Tau_Tru_charge = 0;
   Tau_Tru_pdgId = 0;
   Tau_Tru_tauEta = 0;
   Tau_Tru_tauPhi = 0;
   Tau_Tru_tauE = 0;
   Tau_Tru_tauET = 0;
   Tau_Tru_tauPT = 0;
   Tau_Tru_tauPx = 0;
   Tau_Tru_tauPy = 0;
   Tau_Tru_tauPz = 0;
   Tau_Tru_stablesEAll = 0;
   Tau_Tru_stablesEAllCharged = 0;
   Tau_Tru_stablesEAllPion = 0;
   Tau_Tru_stablesEAllKaon = 0;
   Tau_Tru_stablesEAllOtherCh = 0;
   Tau_Tru_stablesEAllNeut = 0;
   Tau_Tru_stablesEAllPiZ = 0;
   Tau_Tru_stablesEAllGamma = 0;
   Tau_Tru_stablesEAllOtherNeut = 0;
   Tau_Tru_stablesETAll = 0;
   Tau_Tru_stablesETAllCharged = 0;
   Tau_Tru_stablesETAllPion = 0;
   Tau_Tru_stablesETAllKaon = 0;
   Tau_Tru_stablesETAllOtherCh = 0;
   Tau_Tru_stablesETAllNeut = 0;
   Tau_Tru_stablesETAllPiZ = 0;
   Tau_Tru_stablesETAllGamma = 0;
   Tau_Tru_stablesETAllOtherNeut = 0;
   Tau_Tru_stablesPTAll = 0;
   Tau_Tru_stablesPTAllCharged = 0;
   Tau_Tru_stablesPTAllPion = 0;
   Tau_Tru_stablesPTAllKaon = 0;
   Tau_Tru_stablesPTAllOtherCh = 0;
   Tau_Tru_stablesPTAllNeut = 0;
   Tau_Tru_stablesPTAllPiZ = 0;
   Tau_Tru_stablesPTAllGamma = 0;
   Tau_Tru_stablesPTAllOtherNeut = 0;
   PJet_FlatBTagged = 0;
   PJet_FlatBTaggedPlus = 0;
   PJet_FlatBTaggedMinus = 0;
   PJet_NonLooseJet = 0;
   PJet_GoodJet = 0;
   El_Loose = 0;
   El_GoodLepton = 0;
   Mu_Loose = 0;
   Mu_GoodLepton = 0;
   Mu_MuonInJet = 0;
   PJet_TRFIP2DWp3_1incl = 0;
   PJet_TRFIP2DWp3_1excl = 0;
   PJet_TRFIP2DWp3_2incl = 0;
   PJet_TRFIP2DWp3_2excl = 0;
   PJet_TRFSV1Wp3_1incl = 0;
   PJet_TRFSV1Wp3_1excl = 0;
   PJet_TRFSV1Wp3_2incl = 0;
   PJet_TRFSV1Wp3_2excl = 0;
   PJet_TRFIP2DWp3_Plus_1incl = 0;
   PJet_TRFIP2DWp3_Plus_1excl = 0;
   PJet_TRFIP2DWp3_Plus_2incl = 0;
   PJet_TRFIP2DWp3_Plus_2excl = 0;
   PJet_TRFSV1Wp3_Plus_1incl = 0;
   PJet_TRFSV1Wp3_Plus_1excl = 0;
   PJet_TRFSV1Wp3_Plus_2incl = 0;
   PJet_TRFSV1Wp3_Plus_2excl = 0;
   PJet_TRFIP2DWp3_Minus_1incl = 0;
   PJet_TRFIP2DWp3_Minus_1excl = 0;
   PJet_TRFIP2DWp3_Minus_2incl = 0;
   PJet_TRFIP2DWp3_Minus_2excl = 0;
   PJet_TRFSV1Wp3_Minus_1incl = 0;
   PJet_TRFSV1Wp3_Minus_1excl = 0;
   PJet_TRFSV1Wp3_Minus_2incl = 0;
   PJet_TRFSV1Wp3_Minus_2excl = 0;
   PJet_Tru_Matched = 0;
   PJet_Tru_barcode = 0;
   PJet_TruJet_Matched = 0;
   El_Tru_Matched = 0;
   El_Tru_barcode = 0;
   El_TruJet_Matched = 0;
   El_TruTau_Matched = 0;
   Ph_Tru_Matched = 0;
   Ph_Tru_barcode = 0;
   Mu_Tru_Matched = 0;
   Mu_Tru_barcode = 0;
   Tau_Tru_Matched = 0;
   Tau_Tru_compPartIndex = 0;
   Tau_Tru_test_FailHasCharge = 0;
   Tau_Tru_stablesN = 0;
   Tau_Tru_nDaughters = 0;
   Tau_Tru_nStable = 0;
   Tau_Tru_nStableCharged = 0;
   Tau_Tru_nStablePion = 0;
   Tau_Tru_nStableKaon = 0;
   Tau_Tru_nStableOtherCh = 0;
   Tau_Tru_nStableNeut = 0;
   Tau_Tru_nStablePiZ = 0;
   Tau_Tru_nStableGamma = 0;
   Tau_Tru_nStableOtherNeut = 0;
   Tau_Tru_nStableCounted = 0;
   Tau_Tru_nStableNoneMissing = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("AtlfastMissingEt_ex", &AtlfastMissingEt_ex, &b_AtlfastMissingEt_ex);
   fChain->SetBranchAddress("AtlfastMissingEt_ey", &AtlfastMissingEt_ey, &b_AtlfastMissingEt_ey);
   fChain->SetBranchAddress("AtlfastMissingEt_et", &AtlfastMissingEt_et, &b_AtlfastMissingEt_et);
   fChain->SetBranchAddress("AtlfastMissingEt_sumet", &AtlfastMissingEt_sumet, &b_AtlfastMissingEt_sumet);
   fChain->SetBranchAddress("MissingEx", &MissingEx, &b_MissingEx);
   fChain->SetBranchAddress("MissingEy", &MissingEy, &b_MissingEy);
   fChain->SetBranchAddress("MissingEt", &MissingEt, &b_MissingEt);
   fChain->SetBranchAddress("SumEt", &SumEt, &b_SumEt);
   fChain->SetBranchAddress("eventWeight", &eventWeight, &b_eventWeight);
   fChain->SetBranchAddress("eventWeightMCatNLO", &eventWeightMCatNLO, &b_eventWeightMCatNLO);
   fChain->SetBranchAddress("eventWeightIP2DWp3_1incl", &eventWeightIP2DWp3_1incl, &b_eventWeightIP2DWp3_1incl);
   fChain->SetBranchAddress("eventWeightIP2DWp3_1excl", &eventWeightIP2DWp3_1excl, &b_eventWeightIP2DWp3_1excl);
   fChain->SetBranchAddress("eventWeightIP2DWp3_2incl", &eventWeightIP2DWp3_2incl, &b_eventWeightIP2DWp3_2incl);
   fChain->SetBranchAddress("eventWeightIP2DWp3_2excl", &eventWeightIP2DWp3_2excl, &b_eventWeightIP2DWp3_2excl);
   fChain->SetBranchAddress("eventWeightSV1Wp3_1incl", &eventWeightSV1Wp3_1incl, &b_eventWeightSV1Wp3_1incl);
   fChain->SetBranchAddress("eventWeightSV1Wp3_1excl", &eventWeightSV1Wp3_1excl, &b_eventWeightSV1Wp3_1excl);
   fChain->SetBranchAddress("eventWeightSV1Wp3_2incl", &eventWeightSV1Wp3_2incl, &b_eventWeightSV1Wp3_2incl);
   fChain->SetBranchAddress("eventWeightSV1Wp3_2excl", &eventWeightSV1Wp3_2excl, &b_eventWeightSV1Wp3_2excl);
   fChain->SetBranchAddress("eventWeightIP2DWp3_1incl_Plus", &eventWeightIP2DWp3_1incl_Plus, &b_eventWeightIP2DWp3_1incl_Plus);
   fChain->SetBranchAddress("eventWeightIP2DWp3_1excl_Plus", &eventWeightIP2DWp3_1excl_Plus, &b_eventWeightIP2DWp3_1excl_Plus);
   fChain->SetBranchAddress("eventWeightIP2DWp3_2incl_Plus", &eventWeightIP2DWp3_2incl_Plus, &b_eventWeightIP2DWp3_2incl_Plus);
   fChain->SetBranchAddress("eventWeightIP2DWp3_2excl_Plus", &eventWeightIP2DWp3_2excl_Plus, &b_eventWeightIP2DWp3_2excl_Plus);
   fChain->SetBranchAddress("eventWeightSV1Wp3_1incl_Plus", &eventWeightSV1Wp3_1incl_Plus, &b_eventWeightSV1Wp3_1incl_Plus);
   fChain->SetBranchAddress("eventWeightSV1Wp3_1excl_Plus", &eventWeightSV1Wp3_1excl_Plus, &b_eventWeightSV1Wp3_1excl_Plus);
   fChain->SetBranchAddress("eventWeightSV1Wp3_2incl_Plus", &eventWeightSV1Wp3_2incl_Plus, &b_eventWeightSV1Wp3_2incl_Plus);
   fChain->SetBranchAddress("eventWeightSV1Wp3_2excl_Plus", &eventWeightSV1Wp3_2excl_Plus, &b_eventWeightSV1Wp3_2excl_Plus);
   fChain->SetBranchAddress("eventWeightIP2DWp3_1incl_Minus", &eventWeightIP2DWp3_1incl_Minus, &b_eventWeightIP2DWp3_1incl_Minus);
   fChain->SetBranchAddress("eventWeightIP2DWp3_1excl_Minus", &eventWeightIP2DWp3_1excl_Minus, &b_eventWeightIP2DWp3_1excl_Minus);
   fChain->SetBranchAddress("eventWeightIP2DWp3_2incl_Minus", &eventWeightIP2DWp3_2incl_Minus, &b_eventWeightIP2DWp3_2incl_Minus);
   fChain->SetBranchAddress("eventWeightIP2DWp3_2excl_Minus", &eventWeightIP2DWp3_2excl_Minus, &b_eventWeightIP2DWp3_2excl_Minus);
   fChain->SetBranchAddress("eventWeightSV1Wp3_1incl_Minus", &eventWeightSV1Wp3_1incl_Minus, &b_eventWeightSV1Wp3_1incl_Minus);
   fChain->SetBranchAddress("eventWeightSV1Wp3_1excl_Minus", &eventWeightSV1Wp3_1excl_Minus, &b_eventWeightSV1Wp3_1excl_Minus);
   fChain->SetBranchAddress("eventWeightSV1Wp3_2incl_Minus", &eventWeightSV1Wp3_2incl_Minus, &b_eventWeightSV1Wp3_2incl_Minus);
   fChain->SetBranchAddress("eventWeightSV1Wp3_2excl_Minus", &eventWeightSV1Wp3_2excl_Minus, &b_eventWeightSV1Wp3_2excl_Minus);
   fChain->SetBranchAddress("MTotal", &MTotal, &b_MTotal);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("Event_Sphericity", &Event_Sphericity, &b_Event_Sphericity);
   fChain->SetBranchAddress("Event_Aplanarity", &Event_Aplanarity, &b_Event_Aplanarity);
   fChain->SetBranchAddress("Event_Planarity", &Event_Planarity, &b_Event_Planarity);
   fChain->SetBranchAddress("N", &N, &b_N);
   fChain->SetBranchAddress("PJet_N", &PJet_N, &b_PJet_N);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("N", &N, &b_N);
   fChain->SetBranchAddress("GoodLepton_N", &GoodLepton_N, &b_GoodLepton_N);
   fChain->SetBranchAddress("GoodJet_N", &GoodJet_N, &b_GoodJet_N);
   fChain->SetBranchAddress("NonLooseJet_N", &NonLooseJet_N, &b_NonLooseJet_N);
   fChain->SetBranchAddress("El_N", &El_N, &b_El_N);
   fChain->SetBranchAddress("Mu_N", &Mu_N, &b_Mu_N);
   fChain->SetBranchAddress("Ph_N", &Ph_N, &b_Ph_N);
   fChain->SetBranchAddress("Tau_N", &Tau_N, &b_Tau_N);
   fChain->SetBranchAddress("VectSumAll_N", &VectSumAll_N, &b_VectSumAll_N);
   fChain->SetBranchAddress("ForwardJet_N", &ForwardJet_N, &b_ForwardJet_N);
   fChain->SetBranchAddress("PJet_E", &PJet_E, &b_PJet_E);
   fChain->SetBranchAddress("PJet_p_T", &PJet_p_T, &b_PJet_p_T);
   fChain->SetBranchAddress("PJet_phi", &PJet_phi, &b_PJet_phi);
   fChain->SetBranchAddress("PJet_eta", &PJet_eta, &b_PJet_eta);
   fChain->SetBranchAddress("PJet_px", &PJet_px, &b_PJet_px);
   fChain->SetBranchAddress("PJet_py", &PJet_py, &b_PJet_py);
   fChain->SetBranchAddress("PJet_pz", &PJet_pz, &b_PJet_pz);
   fChain->SetBranchAddress("PJet_m", &PJet_m, &b_PJet_m);
   fChain->SetBranchAddress("PJet_charge", &PJet_charge, &b_PJet_charge);
   fChain->SetBranchAddress("PJet_pdgId", &PJet_pdgId, &b_PJet_pdgId);
   fChain->SetBranchAddress("PJet_MCCalibFn", &PJet_MCCalibFn, &b_PJet_MCCalibFn);
   fChain->SetBranchAddress("PJet_MCCalib", &PJet_MCCalib, &b_PJet_MCCalib);
   fChain->SetBranchAddress("PJet_MCCalibError", &PJet_MCCalibError, &b_PJet_MCCalibError);
   fChain->SetBranchAddress("PJet_BHypoCalFactor", &PJet_BHypoCalFactor, &b_PJet_BHypoCalFactor);
   fChain->SetBranchAddress("PJet_LightHypoCalFactor", &PJet_LightHypoCalFactor, &b_PJet_LightHypoCalFactor);
   fChain->SetBranchAddress("PJet_TauHypoCalFactor", &PJet_TauHypoCalFactor, &b_PJet_TauHypoCalFactor);
   fChain->SetBranchAddress("PJet_Tau1P3PHypoCalFactor", &PJet_Tau1P3PHypoCalFactor, &b_PJet_Tau1P3PHypoCalFactor);
   fChain->SetBranchAddress("PJet_DRTruB", &PJet_DRTruB, &b_PJet_DRTruB);
   fChain->SetBranchAddress("PJet_DRTruC", &PJet_DRTruC, &b_PJet_DRTruC);
   fChain->SetBranchAddress("PJet_DRTruT", &PJet_DRTruT, &b_PJet_DRTruT);
   fChain->SetBranchAddress("El_E", &El_E, &b_El_E);
   fChain->SetBranchAddress("El_p_T", &El_p_T, &b_El_p_T);
   fChain->SetBranchAddress("El_phi", &El_phi, &b_El_phi);
   fChain->SetBranchAddress("El_eta", &El_eta, &b_El_eta);
   fChain->SetBranchAddress("El_px", &El_px, &b_El_px);
   fChain->SetBranchAddress("El_py", &El_py, &b_El_py);
   fChain->SetBranchAddress("El_pz", &El_pz, &b_El_pz);
   fChain->SetBranchAddress("El_m", &El_m, &b_El_m);
   fChain->SetBranchAddress("El_charge", &El_charge, &b_El_charge);
   fChain->SetBranchAddress("El_pdgId", &El_pdgId, &b_El_pdgId);
   fChain->SetBranchAddress("Mu_E", &Mu_E, &b_Mu_E);
   fChain->SetBranchAddress("Mu_p_T", &Mu_p_T, &b_Mu_p_T);
   fChain->SetBranchAddress("Mu_phi", &Mu_phi, &b_Mu_phi);
   fChain->SetBranchAddress("Mu_eta", &Mu_eta, &b_Mu_eta);
   fChain->SetBranchAddress("Mu_px", &Mu_px, &b_Mu_px);
   fChain->SetBranchAddress("Mu_py", &Mu_py, &b_Mu_py);
   fChain->SetBranchAddress("Mu_pz", &Mu_pz, &b_Mu_pz);
   fChain->SetBranchAddress("Mu_m", &Mu_m, &b_Mu_m);
   fChain->SetBranchAddress("Mu_charge", &Mu_charge, &b_Mu_charge);
   fChain->SetBranchAddress("Mu_pdgId", &Mu_pdgId, &b_Mu_pdgId);
   fChain->SetBranchAddress("Ph_E", &Ph_E, &b_Ph_E);
   fChain->SetBranchAddress("Ph_p_T", &Ph_p_T, &b_Ph_p_T);
   fChain->SetBranchAddress("Ph_phi", &Ph_phi, &b_Ph_phi);
   fChain->SetBranchAddress("Ph_eta", &Ph_eta, &b_Ph_eta);
   fChain->SetBranchAddress("Ph_px", &Ph_px, &b_Ph_px);
   fChain->SetBranchAddress("Ph_py", &Ph_py, &b_Ph_py);
   fChain->SetBranchAddress("Ph_pz", &Ph_pz, &b_Ph_pz);
   fChain->SetBranchAddress("Ph_m", &Ph_m, &b_Ph_m);
   fChain->SetBranchAddress("Ph_charge", &Ph_charge, &b_Ph_charge);
   fChain->SetBranchAddress("Ph_pdgId", &Ph_pdgId, &b_Ph_pdgId);
   fChain->SetBranchAddress("Tau_E", &Tau_E, &b_Tau_E);
   fChain->SetBranchAddress("Tau_p_T", &Tau_p_T, &b_Tau_p_T);
   fChain->SetBranchAddress("Tau_phi", &Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_eta", &Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_px", &Tau_px, &b_Tau_px);
   fChain->SetBranchAddress("Tau_py", &Tau_py, &b_Tau_py);
   fChain->SetBranchAddress("Tau_pz", &Tau_pz, &b_Tau_pz);
   fChain->SetBranchAddress("Tau_m", &Tau_m, &b_Tau_m);
   fChain->SetBranchAddress("Tau_charge", &Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_pdgId", &Tau_pdgId, &b_Tau_pdgId);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3", &PJet_TRFIP2DWp3, &b_PJet_TRFIP2DWp3);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3", &PJet_TRFSV1Wp3, &b_PJet_TRFSV1Wp3);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Plus", &PJet_TRFIP2DWp3_Plus, &b_PJet_TRFIP2DWp3_Plus);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Plus", &PJet_TRFSV1Wp3_Plus, &b_PJet_TRFSV1Wp3_Plus);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Minus", &PJet_TRFIP2DWp3_Minus, &b_PJet_TRFIP2DWp3_Minus);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Minus", &PJet_TRFSV1Wp3_Minus, &b_PJet_TRFSV1Wp3_Minus);
   fChain->SetBranchAddress("El_TransM", &El_TransM, &b_El_TransM);
   fChain->SetBranchAddress("Mu_TransM", &Mu_TransM, &b_Mu_TransM);
   fChain->SetBranchAddress("VectSumAll_E", &VectSumAll_E, &b_VectSumAll_E);
   fChain->SetBranchAddress("VectSumAll_p_T", &VectSumAll_p_T, &b_VectSumAll_p_T);
   fChain->SetBranchAddress("VectSumAll_phi", &VectSumAll_phi, &b_VectSumAll_phi);
   fChain->SetBranchAddress("VectSumAll_eta", &VectSumAll_eta, &b_VectSumAll_eta);
   fChain->SetBranchAddress("VectSumAll_px", &VectSumAll_px, &b_VectSumAll_px);
   fChain->SetBranchAddress("VectSumAll_py", &VectSumAll_py, &b_VectSumAll_py);
   fChain->SetBranchAddress("VectSumAll_pz", &VectSumAll_pz, &b_VectSumAll_pz);
   fChain->SetBranchAddress("VectSumAll_m", &VectSumAll_m, &b_VectSumAll_m);
   fChain->SetBranchAddress("VectSumAll_charge", &VectSumAll_charge, &b_VectSumAll_charge);
   fChain->SetBranchAddress("VectSumAll_pdgId", &VectSumAll_pdgId, &b_VectSumAll_pdgId);
   fChain->SetBranchAddress("PJet_cosThetaStar", &PJet_cosThetaStar, &b_PJet_cosThetaStar);
   fChain->SetBranchAddress("PJet_Tru_DeltaR", &PJet_Tru_DeltaR, &b_PJet_Tru_DeltaR);
   fChain->SetBranchAddress("PJet_Tru_E", &PJet_Tru_E, &b_PJet_Tru_E);
   fChain->SetBranchAddress("PJet_Tru_p_T", &PJet_Tru_p_T, &b_PJet_Tru_p_T);
   fChain->SetBranchAddress("PJet_Tru_phi", &PJet_Tru_phi, &b_PJet_Tru_phi);
   fChain->SetBranchAddress("PJet_Tru_eta", &PJet_Tru_eta, &b_PJet_Tru_eta);
   fChain->SetBranchAddress("PJet_Tru_px", &PJet_Tru_px, &b_PJet_Tru_px);
   fChain->SetBranchAddress("PJet_Tru_py", &PJet_Tru_py, &b_PJet_Tru_py);
   fChain->SetBranchAddress("PJet_Tru_pz", &PJet_Tru_pz, &b_PJet_Tru_pz);
   fChain->SetBranchAddress("PJet_Tru_m", &PJet_Tru_m, &b_PJet_Tru_m);
   fChain->SetBranchAddress("PJet_Tru_charge", &PJet_Tru_charge, &b_PJet_Tru_charge);
   fChain->SetBranchAddress("PJet_Tru_pdgId", &PJet_Tru_pdgId, &b_PJet_Tru_pdgId);
   fChain->SetBranchAddress("PJet_Tru_etcone", &PJet_Tru_etcone, &b_PJet_Tru_etcone);
   fChain->SetBranchAddress("PJet_Tru_etcone10", &PJet_Tru_etcone10, &b_PJet_Tru_etcone10);
   fChain->SetBranchAddress("PJet_Tru_etcone20", &PJet_Tru_etcone20, &b_PJet_Tru_etcone20);
   fChain->SetBranchAddress("PJet_Tru_etcone30", &PJet_Tru_etcone30, &b_PJet_Tru_etcone30);
   fChain->SetBranchAddress("PJet_Tru_etcone40", &PJet_Tru_etcone40, &b_PJet_Tru_etcone40);
   fChain->SetBranchAddress("PJet_Tru_etcone50", &PJet_Tru_etcone50, &b_PJet_Tru_etcone50);
   fChain->SetBranchAddress("PJet_Tru_etcone60", &PJet_Tru_etcone60, &b_PJet_Tru_etcone60);
   fChain->SetBranchAddress("PJet_Tru_etcone70", &PJet_Tru_etcone70, &b_PJet_Tru_etcone70);
   fChain->SetBranchAddress("PJet_Tru_status", &PJet_Tru_status, &b_PJet_Tru_status);
   fChain->SetBranchAddress("PJet_TruJet_DeltaR", &PJet_TruJet_DeltaR, &b_PJet_TruJet_DeltaR);
   fChain->SetBranchAddress("PJet_TruJet_E", &PJet_TruJet_E, &b_PJet_TruJet_E);
   fChain->SetBranchAddress("PJet_TruJet_p_T", &PJet_TruJet_p_T, &b_PJet_TruJet_p_T);
   fChain->SetBranchAddress("PJet_TruJet_phi", &PJet_TruJet_phi, &b_PJet_TruJet_phi);
   fChain->SetBranchAddress("PJet_TruJet_eta", &PJet_TruJet_eta, &b_PJet_TruJet_eta);
   fChain->SetBranchAddress("PJet_TruJet_px", &PJet_TruJet_px, &b_PJet_TruJet_px);
   fChain->SetBranchAddress("PJet_TruJet_py", &PJet_TruJet_py, &b_PJet_TruJet_py);
   fChain->SetBranchAddress("PJet_TruJet_pz", &PJet_TruJet_pz, &b_PJet_TruJet_pz);
   fChain->SetBranchAddress("PJet_TruJet_m", &PJet_TruJet_m, &b_PJet_TruJet_m);
   fChain->SetBranchAddress("PJet_TruJet_charge", &PJet_TruJet_charge, &b_PJet_TruJet_charge);
   fChain->SetBranchAddress("PJet_TruJet_pdgId", &PJet_TruJet_pdgId, &b_PJet_TruJet_pdgId);
   fChain->SetBranchAddress("El_Tru_DeltaR", &El_Tru_DeltaR, &b_El_Tru_DeltaR);
   fChain->SetBranchAddress("El_Tru_E", &El_Tru_E, &b_El_Tru_E);
   fChain->SetBranchAddress("El_Tru_p_T", &El_Tru_p_T, &b_El_Tru_p_T);
   fChain->SetBranchAddress("El_Tru_phi", &El_Tru_phi, &b_El_Tru_phi);
   fChain->SetBranchAddress("El_Tru_eta", &El_Tru_eta, &b_El_Tru_eta);
   fChain->SetBranchAddress("El_Tru_px", &El_Tru_px, &b_El_Tru_px);
   fChain->SetBranchAddress("El_Tru_py", &El_Tru_py, &b_El_Tru_py);
   fChain->SetBranchAddress("El_Tru_pz", &El_Tru_pz, &b_El_Tru_pz);
   fChain->SetBranchAddress("El_Tru_m", &El_Tru_m, &b_El_Tru_m);
   fChain->SetBranchAddress("El_Tru_charge", &El_Tru_charge, &b_El_Tru_charge);
   fChain->SetBranchAddress("El_Tru_pdgId", &El_Tru_pdgId, &b_El_Tru_pdgId);
   fChain->SetBranchAddress("El_Tru_etcone", &El_Tru_etcone, &b_El_Tru_etcone);
   fChain->SetBranchAddress("El_Tru_etcone10", &El_Tru_etcone10, &b_El_Tru_etcone10);
   fChain->SetBranchAddress("El_Tru_etcone20", &El_Tru_etcone20, &b_El_Tru_etcone20);
   fChain->SetBranchAddress("El_Tru_etcone30", &El_Tru_etcone30, &b_El_Tru_etcone30);
   fChain->SetBranchAddress("El_Tru_etcone40", &El_Tru_etcone40, &b_El_Tru_etcone40);
   fChain->SetBranchAddress("El_Tru_etcone50", &El_Tru_etcone50, &b_El_Tru_etcone50);
   fChain->SetBranchAddress("El_Tru_etcone60", &El_Tru_etcone60, &b_El_Tru_etcone60);
   fChain->SetBranchAddress("El_Tru_etcone70", &El_Tru_etcone70, &b_El_Tru_etcone70);
   fChain->SetBranchAddress("El_Tru_status", &El_Tru_status, &b_El_Tru_status);
   fChain->SetBranchAddress("El_TruJet_DeltaR", &El_TruJet_DeltaR, &b_El_TruJet_DeltaR);
   fChain->SetBranchAddress("El_TruTau_DeltaR", &El_TruTau_DeltaR, &b_El_TruTau_DeltaR);
   fChain->SetBranchAddress("Ph_Tru_DeltaR", &Ph_Tru_DeltaR, &b_Ph_Tru_DeltaR);
   fChain->SetBranchAddress("Ph_Tru_E", &Ph_Tru_E, &b_Ph_Tru_E);
   fChain->SetBranchAddress("Ph_Tru_p_T", &Ph_Tru_p_T, &b_Ph_Tru_p_T);
   fChain->SetBranchAddress("Ph_Tru_phi", &Ph_Tru_phi, &b_Ph_Tru_phi);
   fChain->SetBranchAddress("Ph_Tru_eta", &Ph_Tru_eta, &b_Ph_Tru_eta);
   fChain->SetBranchAddress("Ph_Tru_px", &Ph_Tru_px, &b_Ph_Tru_px);
   fChain->SetBranchAddress("Ph_Tru_py", &Ph_Tru_py, &b_Ph_Tru_py);
   fChain->SetBranchAddress("Ph_Tru_pz", &Ph_Tru_pz, &b_Ph_Tru_pz);
   fChain->SetBranchAddress("Ph_Tru_m", &Ph_Tru_m, &b_Ph_Tru_m);
   fChain->SetBranchAddress("Ph_Tru_charge", &Ph_Tru_charge, &b_Ph_Tru_charge);
   fChain->SetBranchAddress("Ph_Tru_pdgId", &Ph_Tru_pdgId, &b_Ph_Tru_pdgId);
   fChain->SetBranchAddress("Ph_Tru_etcone", &Ph_Tru_etcone, &b_Ph_Tru_etcone);
   fChain->SetBranchAddress("Ph_Tru_etcone10", &Ph_Tru_etcone10, &b_Ph_Tru_etcone10);
   fChain->SetBranchAddress("Ph_Tru_etcone20", &Ph_Tru_etcone20, &b_Ph_Tru_etcone20);
   fChain->SetBranchAddress("Ph_Tru_etcone30", &Ph_Tru_etcone30, &b_Ph_Tru_etcone30);
   fChain->SetBranchAddress("Ph_Tru_etcone40", &Ph_Tru_etcone40, &b_Ph_Tru_etcone40);
   fChain->SetBranchAddress("Ph_Tru_etcone50", &Ph_Tru_etcone50, &b_Ph_Tru_etcone50);
   fChain->SetBranchAddress("Ph_Tru_etcone60", &Ph_Tru_etcone60, &b_Ph_Tru_etcone60);
   fChain->SetBranchAddress("Ph_Tru_etcone70", &Ph_Tru_etcone70, &b_Ph_Tru_etcone70);
   fChain->SetBranchAddress("Ph_Tru_status", &Ph_Tru_status, &b_Ph_Tru_status);
   fChain->SetBranchAddress("Mu_Tru_DeltaR", &Mu_Tru_DeltaR, &b_Mu_Tru_DeltaR);
   fChain->SetBranchAddress("Mu_Tru_E", &Mu_Tru_E, &b_Mu_Tru_E);
   fChain->SetBranchAddress("Mu_Tru_p_T", &Mu_Tru_p_T, &b_Mu_Tru_p_T);
   fChain->SetBranchAddress("Mu_Tru_phi", &Mu_Tru_phi, &b_Mu_Tru_phi);
   fChain->SetBranchAddress("Mu_Tru_eta", &Mu_Tru_eta, &b_Mu_Tru_eta);
   fChain->SetBranchAddress("Mu_Tru_px", &Mu_Tru_px, &b_Mu_Tru_px);
   fChain->SetBranchAddress("Mu_Tru_py", &Mu_Tru_py, &b_Mu_Tru_py);
   fChain->SetBranchAddress("Mu_Tru_pz", &Mu_Tru_pz, &b_Mu_Tru_pz);
   fChain->SetBranchAddress("Mu_Tru_m", &Mu_Tru_m, &b_Mu_Tru_m);
   fChain->SetBranchAddress("Mu_Tru_charge", &Mu_Tru_charge, &b_Mu_Tru_charge);
   fChain->SetBranchAddress("Mu_Tru_pdgId", &Mu_Tru_pdgId, &b_Mu_Tru_pdgId);
   fChain->SetBranchAddress("Mu_Tru_etcone", &Mu_Tru_etcone, &b_Mu_Tru_etcone);
   fChain->SetBranchAddress("Mu_Tru_etcone10", &Mu_Tru_etcone10, &b_Mu_Tru_etcone10);
   fChain->SetBranchAddress("Mu_Tru_etcone20", &Mu_Tru_etcone20, &b_Mu_Tru_etcone20);
   fChain->SetBranchAddress("Mu_Tru_etcone30", &Mu_Tru_etcone30, &b_Mu_Tru_etcone30);
   fChain->SetBranchAddress("Mu_Tru_etcone40", &Mu_Tru_etcone40, &b_Mu_Tru_etcone40);
   fChain->SetBranchAddress("Mu_Tru_etcone50", &Mu_Tru_etcone50, &b_Mu_Tru_etcone50);
   fChain->SetBranchAddress("Mu_Tru_etcone60", &Mu_Tru_etcone60, &b_Mu_Tru_etcone60);
   fChain->SetBranchAddress("Mu_Tru_etcone70", &Mu_Tru_etcone70, &b_Mu_Tru_etcone70);
   fChain->SetBranchAddress("Mu_Tru_status", &Mu_Tru_status, &b_Mu_Tru_status);
   fChain->SetBranchAddress("Tau_Tru_DeltaR", &Tau_Tru_DeltaR, &b_Tau_Tru_DeltaR);
   fChain->SetBranchAddress("Tau_Tru_E", &Tau_Tru_E, &b_Tau_Tru_E);
   fChain->SetBranchAddress("Tau_Tru_p_T", &Tau_Tru_p_T, &b_Tau_Tru_p_T);
   fChain->SetBranchAddress("Tau_Tru_phi", &Tau_Tru_phi, &b_Tau_Tru_phi);
   fChain->SetBranchAddress("Tau_Tru_eta", &Tau_Tru_eta, &b_Tau_Tru_eta);
   fChain->SetBranchAddress("Tau_Tru_px", &Tau_Tru_px, &b_Tau_Tru_px);
   fChain->SetBranchAddress("Tau_Tru_py", &Tau_Tru_py, &b_Tau_Tru_py);
   fChain->SetBranchAddress("Tau_Tru_pz", &Tau_Tru_pz, &b_Tau_Tru_pz);
   fChain->SetBranchAddress("Tau_Tru_m", &Tau_Tru_m, &b_Tau_Tru_m);
   fChain->SetBranchAddress("Tau_Tru_charge", &Tau_Tru_charge, &b_Tau_Tru_charge);
   fChain->SetBranchAddress("Tau_Tru_pdgId", &Tau_Tru_pdgId, &b_Tau_Tru_pdgId);
   fChain->SetBranchAddress("Tau_Tru_tauEta", &Tau_Tru_tauEta, &b_Tau_Tru_tauEta);
   fChain->SetBranchAddress("Tau_Tru_tauPhi", &Tau_Tru_tauPhi, &b_Tau_Tru_tauPhi);
   fChain->SetBranchAddress("Tau_Tru_tauE", &Tau_Tru_tauE, &b_Tau_Tru_tauE);
   fChain->SetBranchAddress("Tau_Tru_tauET", &Tau_Tru_tauET, &b_Tau_Tru_tauET);
   fChain->SetBranchAddress("Tau_Tru_tauPT", &Tau_Tru_tauPT, &b_Tau_Tru_tauPT);
   fChain->SetBranchAddress("Tau_Tru_tauPx", &Tau_Tru_tauPx, &b_Tau_Tru_tauPx);
   fChain->SetBranchAddress("Tau_Tru_tauPy", &Tau_Tru_tauPy, &b_Tau_Tru_tauPy);
   fChain->SetBranchAddress("Tau_Tru_tauPz", &Tau_Tru_tauPz, &b_Tau_Tru_tauPz);
   fChain->SetBranchAddress("Tau_Tru_stablesEAll", &Tau_Tru_stablesEAll, &b_Tau_Tru_stablesEAll);
   fChain->SetBranchAddress("Tau_Tru_stablesEAllCharged", &Tau_Tru_stablesEAllCharged, &b_Tau_Tru_stablesEAllCharged);
   fChain->SetBranchAddress("Tau_Tru_stablesEAllPion", &Tau_Tru_stablesEAllPion, &b_Tau_Tru_stablesEAllPion);
   fChain->SetBranchAddress("Tau_Tru_stablesEAllKaon", &Tau_Tru_stablesEAllKaon, &b_Tau_Tru_stablesEAllKaon);
   fChain->SetBranchAddress("Tau_Tru_stablesEAllOtherCh", &Tau_Tru_stablesEAllOtherCh, &b_Tau_Tru_stablesEAllOtherCh);
   fChain->SetBranchAddress("Tau_Tru_stablesEAllNeut", &Tau_Tru_stablesEAllNeut, &b_Tau_Tru_stablesEAllNeut);
   fChain->SetBranchAddress("Tau_Tru_stablesEAllPiZ", &Tau_Tru_stablesEAllPiZ, &b_Tau_Tru_stablesEAllPiZ);
   fChain->SetBranchAddress("Tau_Tru_stablesEAllGamma", &Tau_Tru_stablesEAllGamma, &b_Tau_Tru_stablesEAllGamma);
   fChain->SetBranchAddress("Tau_Tru_stablesEAllOtherNeut", &Tau_Tru_stablesEAllOtherNeut, &b_Tau_Tru_stablesEAllOtherNeut);
   fChain->SetBranchAddress("Tau_Tru_stablesETAll", &Tau_Tru_stablesETAll, &b_Tau_Tru_stablesETAll);
   fChain->SetBranchAddress("Tau_Tru_stablesETAllCharged", &Tau_Tru_stablesETAllCharged, &b_Tau_Tru_stablesETAllCharged);
   fChain->SetBranchAddress("Tau_Tru_stablesETAllPion", &Tau_Tru_stablesETAllPion, &b_Tau_Tru_stablesETAllPion);
   fChain->SetBranchAddress("Tau_Tru_stablesETAllKaon", &Tau_Tru_stablesETAllKaon, &b_Tau_Tru_stablesETAllKaon);
   fChain->SetBranchAddress("Tau_Tru_stablesETAllOtherCh", &Tau_Tru_stablesETAllOtherCh, &b_Tau_Tru_stablesETAllOtherCh);
   fChain->SetBranchAddress("Tau_Tru_stablesETAllNeut", &Tau_Tru_stablesETAllNeut, &b_Tau_Tru_stablesETAllNeut);
   fChain->SetBranchAddress("Tau_Tru_stablesETAllPiZ", &Tau_Tru_stablesETAllPiZ, &b_Tau_Tru_stablesETAllPiZ);
   fChain->SetBranchAddress("Tau_Tru_stablesETAllGamma", &Tau_Tru_stablesETAllGamma, &b_Tau_Tru_stablesETAllGamma);
   fChain->SetBranchAddress("Tau_Tru_stablesETAllOtherNeut", &Tau_Tru_stablesETAllOtherNeut, &b_Tau_Tru_stablesETAllOtherNeut);
   fChain->SetBranchAddress("Tau_Tru_stablesPTAll", &Tau_Tru_stablesPTAll, &b_Tau_Tru_stablesPTAll);
   fChain->SetBranchAddress("Tau_Tru_stablesPTAllCharged", &Tau_Tru_stablesPTAllCharged, &b_Tau_Tru_stablesPTAllCharged);
   fChain->SetBranchAddress("Tau_Tru_stablesPTAllPion", &Tau_Tru_stablesPTAllPion, &b_Tau_Tru_stablesPTAllPion);
   fChain->SetBranchAddress("Tau_Tru_stablesPTAllKaon", &Tau_Tru_stablesPTAllKaon, &b_Tau_Tru_stablesPTAllKaon);
   fChain->SetBranchAddress("Tau_Tru_stablesPTAllOtherCh", &Tau_Tru_stablesPTAllOtherCh, &b_Tau_Tru_stablesPTAllOtherCh);
   fChain->SetBranchAddress("Tau_Tru_stablesPTAllNeut", &Tau_Tru_stablesPTAllNeut, &b_Tau_Tru_stablesPTAllNeut);
   fChain->SetBranchAddress("Tau_Tru_stablesPTAllPiZ", &Tau_Tru_stablesPTAllPiZ, &b_Tau_Tru_stablesPTAllPiZ);
   fChain->SetBranchAddress("Tau_Tru_stablesPTAllGamma", &Tau_Tru_stablesPTAllGamma, &b_Tau_Tru_stablesPTAllGamma);
   fChain->SetBranchAddress("Tau_Tru_stablesPTAllOtherNeut", &Tau_Tru_stablesPTAllOtherNeut, &b_Tau_Tru_stablesPTAllOtherNeut);
   fChain->SetBranchAddress("PJet_FlatBTagged", &PJet_FlatBTagged, &b_PJet_FlatBTagged);
   fChain->SetBranchAddress("PJet_FlatBTaggedPlus", &PJet_FlatBTaggedPlus, &b_PJet_FlatBTaggedPlus);
   fChain->SetBranchAddress("PJet_FlatBTaggedMinus", &PJet_FlatBTaggedMinus, &b_PJet_FlatBTaggedMinus);
   fChain->SetBranchAddress("PJet_NonLooseJet", &PJet_NonLooseJet, &b_PJet_NonLooseJet);
   fChain->SetBranchAddress("PJet_GoodJet", &PJet_GoodJet, &b_PJet_GoodJet);
   fChain->SetBranchAddress("El_Loose", &El_Loose, &b_El_Loose);
   fChain->SetBranchAddress("El_GoodLepton", &El_GoodLepton, &b_El_GoodLepton);
   fChain->SetBranchAddress("Mu_Loose", &Mu_Loose, &b_Mu_Loose);
   fChain->SetBranchAddress("Mu_GoodLepton", &Mu_GoodLepton, &b_Mu_GoodLepton);
   fChain->SetBranchAddress("Mu_MuonInJet", &Mu_MuonInJet, &b_Mu_MuonInJet);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_1incl", &PJet_TRFIP2DWp3_1incl, &b_PJet_TRFIP2DWp3_1incl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_1excl", &PJet_TRFIP2DWp3_1excl, &b_PJet_TRFIP2DWp3_1excl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_2incl", &PJet_TRFIP2DWp3_2incl, &b_PJet_TRFIP2DWp3_2incl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_2excl", &PJet_TRFIP2DWp3_2excl, &b_PJet_TRFIP2DWp3_2excl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_1incl", &PJet_TRFSV1Wp3_1incl, &b_PJet_TRFSV1Wp3_1incl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_1excl", &PJet_TRFSV1Wp3_1excl, &b_PJet_TRFSV1Wp3_1excl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_2incl", &PJet_TRFSV1Wp3_2incl, &b_PJet_TRFSV1Wp3_2incl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_2excl", &PJet_TRFSV1Wp3_2excl, &b_PJet_TRFSV1Wp3_2excl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Plus_1incl", &PJet_TRFIP2DWp3_Plus_1incl, &b_PJet_TRFIP2DWp3_Plus_1incl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Plus_1excl", &PJet_TRFIP2DWp3_Plus_1excl, &b_PJet_TRFIP2DWp3_Plus_1excl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Plus_2incl", &PJet_TRFIP2DWp3_Plus_2incl, &b_PJet_TRFIP2DWp3_Plus_2incl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Plus_2excl", &PJet_TRFIP2DWp3_Plus_2excl, &b_PJet_TRFIP2DWp3_Plus_2excl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Plus_1incl", &PJet_TRFSV1Wp3_Plus_1incl, &b_PJet_TRFSV1Wp3_Plus_1incl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Plus_1excl", &PJet_TRFSV1Wp3_Plus_1excl, &b_PJet_TRFSV1Wp3_Plus_1excl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Plus_2incl", &PJet_TRFSV1Wp3_Plus_2incl, &b_PJet_TRFSV1Wp3_Plus_2incl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Plus_2excl", &PJet_TRFSV1Wp3_Plus_2excl, &b_PJet_TRFSV1Wp3_Plus_2excl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Minus_1incl", &PJet_TRFIP2DWp3_Minus_1incl, &b_PJet_TRFIP2DWp3_Minus_1incl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Minus_1excl", &PJet_TRFIP2DWp3_Minus_1excl, &b_PJet_TRFIP2DWp3_Minus_1excl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Minus_2incl", &PJet_TRFIP2DWp3_Minus_2incl, &b_PJet_TRFIP2DWp3_Minus_2incl);
   fChain->SetBranchAddress("PJet_TRFIP2DWp3_Minus_2excl", &PJet_TRFIP2DWp3_Minus_2excl, &b_PJet_TRFIP2DWp3_Minus_2excl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Minus_1incl", &PJet_TRFSV1Wp3_Minus_1incl, &b_PJet_TRFSV1Wp3_Minus_1incl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Minus_1excl", &PJet_TRFSV1Wp3_Minus_1excl, &b_PJet_TRFSV1Wp3_Minus_1excl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Minus_2incl", &PJet_TRFSV1Wp3_Minus_2incl, &b_PJet_TRFSV1Wp3_Minus_2incl);
   fChain->SetBranchAddress("PJet_TRFSV1Wp3_Minus_2excl", &PJet_TRFSV1Wp3_Minus_2excl, &b_PJet_TRFSV1Wp3_Minus_2excl);
   fChain->SetBranchAddress("PJet_Tru_Matched", &PJet_Tru_Matched, &b_PJet_Tru_Matched);
   fChain->SetBranchAddress("PJet_Tru_barcode", &PJet_Tru_barcode, &b_PJet_Tru_barcode);
   fChain->SetBranchAddress("PJet_TruJet_Matched", &PJet_TruJet_Matched, &b_PJet_TruJet_Matched);
   fChain->SetBranchAddress("El_Tru_Matched", &El_Tru_Matched, &b_El_Tru_Matched);
   fChain->SetBranchAddress("El_Tru_barcode", &El_Tru_barcode, &b_El_Tru_barcode);
   fChain->SetBranchAddress("El_TruJet_Matched", &El_TruJet_Matched, &b_El_TruJet_Matched);
   fChain->SetBranchAddress("El_TruTau_Matched", &El_TruTau_Matched, &b_El_TruTau_Matched);
   fChain->SetBranchAddress("Ph_Tru_Matched", &Ph_Tru_Matched, &b_Ph_Tru_Matched);
   fChain->SetBranchAddress("Ph_Tru_barcode", &Ph_Tru_barcode, &b_Ph_Tru_barcode);
   fChain->SetBranchAddress("Mu_Tru_Matched", &Mu_Tru_Matched, &b_Mu_Tru_Matched);
   fChain->SetBranchAddress("Mu_Tru_barcode", &Mu_Tru_barcode, &b_Mu_Tru_barcode);
   fChain->SetBranchAddress("Tau_Tru_Matched", &Tau_Tru_Matched, &b_Tau_Tru_Matched);
   fChain->SetBranchAddress("Tau_Tru_compPartIndex", &Tau_Tru_compPartIndex, &b_Tau_Tru_compPartIndex);
   fChain->SetBranchAddress("Tau_Tru_test_FailHasCharge", &Tau_Tru_test_FailHasCharge, &b_Tau_Tru_test_FailHasCharge);
   fChain->SetBranchAddress("Tau_Tru_stablesN", &Tau_Tru_stablesN, &b_Tau_Tru_stablesN);
   fChain->SetBranchAddress("Tau_Tru_nDaughters", &Tau_Tru_nDaughters, &b_Tau_Tru_nDaughters);
   fChain->SetBranchAddress("Tau_Tru_nStable", &Tau_Tru_nStable, &b_Tau_Tru_nStable);
   fChain->SetBranchAddress("Tau_Tru_nStableCharged", &Tau_Tru_nStableCharged, &b_Tau_Tru_nStableCharged);
   fChain->SetBranchAddress("Tau_Tru_nStablePion", &Tau_Tru_nStablePion, &b_Tau_Tru_nStablePion);
   fChain->SetBranchAddress("Tau_Tru_nStableKaon", &Tau_Tru_nStableKaon, &b_Tau_Tru_nStableKaon);
   fChain->SetBranchAddress("Tau_Tru_nStableOtherCh", &Tau_Tru_nStableOtherCh, &b_Tau_Tru_nStableOtherCh);
   fChain->SetBranchAddress("Tau_Tru_nStableNeut", &Tau_Tru_nStableNeut, &b_Tau_Tru_nStableNeut);
   fChain->SetBranchAddress("Tau_Tru_nStablePiZ", &Tau_Tru_nStablePiZ, &b_Tau_Tru_nStablePiZ);
   fChain->SetBranchAddress("Tau_Tru_nStableGamma", &Tau_Tru_nStableGamma, &b_Tau_Tru_nStableGamma);
   fChain->SetBranchAddress("Tau_Tru_nStableOtherNeut", &Tau_Tru_nStableOtherNeut, &b_Tau_Tru_nStableOtherNeut);
   fChain->SetBranchAddress("Tau_Tru_nStableCounted", &Tau_Tru_nStableCounted, &b_Tau_Tru_nStableCounted);
   fChain->SetBranchAddress("Tau_Tru_nStableNoneMissing", &Tau_Tru_nStableNoneMissing, &b_Tau_Tru_nStableNoneMissing);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Event_SphericityAxis_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_Event_SphericityAxis_fBits);
   fChain->SetBranchAddress("fX", &fX, &b_Event_SphericityAxis_fX);
   fChain->SetBranchAddress("fY", &fY, &b_Event_SphericityAxis_fY);
   fChain->SetBranchAddress("fZ", &fZ, &b_Event_SphericityAxis_fZ);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Event_SecondAxis_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_Event_SecondAxis_fBits);
   fChain->SetBranchAddress("fX", &fX, &b_Event_SecondAxis_fX);
   fChain->SetBranchAddress("fY", &fY, &b_Event_SecondAxis_fY);
   fChain->SetBranchAddress("fZ", &fZ, &b_Event_SecondAxis_fZ);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Event_ThirdAxis_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_Event_ThirdAxis_fBits);
   fChain->SetBranchAddress("fX", &fX, &b_Event_ThirdAxis_fX);
   fChain->SetBranchAddress("fY", &fY, &b_Event_ThirdAxis_fY);
   fChain->SetBranchAddress("fZ", &fZ, &b_Event_ThirdAxis_fZ);
   fChain->SetBranchAddress("FastSimInstance", &FastSimInstance, &b_FastSimInstance);
   fChain->SetBranchAddress("FastSimCounter", &FastSimCounter, &b_FastSimCounter);
   fChain->SetBranchAddress("FastSimEventAndInstance", &FastSimEventAndInstance, &b_FastSimEventAndInstance);
   fChain->SetBranchAddress("FastSimRunNumber", &FastSimRunNumber, &b_FastSimRunNumber);
   fChain->SetBranchAddress("FastSimNInstance", &FastSimNInstance, &b_FastSimNInstance);
   Notify();
}

Bool_t FastSimTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normaly not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void FastSimTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t FastSimTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef FastSimTree_cxx
