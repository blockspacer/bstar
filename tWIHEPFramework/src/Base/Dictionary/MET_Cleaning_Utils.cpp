#include "SingleTopRootAnalysis/Base/Dictionary/MET_Cleaning_Utils.hpp"
#include <iostream>
#include <cmath>

// Constructor
MET_Cleaning_Utils::MET_Cleaning_Utils():
  GeV(1000)
  ,jet_THRESHOLD(20*GeV)
  ,jet_THRESHOLD_1sigma(4*GeV)
  ,jet_iso(0.1)
{
  eta_min.push_back(0.0);eta_max.push_back(1.45);phi_min.push_back(-0.78847);phi_max.push_back(-0.59213);
}

//destructor
MET_Cleaning_Utils::~MET_Cleaning_Utils(){

}

// Common funtions
float MET_Cleaning_Utils::DeltaPhi(float phi1, float phi2){
  float dphi=fabs(phi1-phi2);
  if (dphi>TMath::Pi())   dphi=2*TMath::Pi()-dphi;
  double sign=1;
  if(phi1-phi2<0) sign=-1;
  return dphi*sign;
}
float MET_Cleaning_Utils::DeltaEta(float eta1, float eta2){
  float deta=fabs(eta1-eta2);
  return deta;
}
float MET_Cleaning_Utils::dr(float phi1, float eta1,float phi2, float eta2){
  return sqrt(pow(DeltaEta(eta1,eta2),2.0)+DeltaPhi(phi1,phi2));
}

float MET_Cleaning_Utils::dphi_MET_MPT(float &met_phi,std::vector<float> *trk_pt,std::vector<float> *trk_eta,
				       std::vector<float> *trk_d0_wrtPV,std::vector<float> *trk_z0_wrtPV,
				       std::vector<float> *trk_theta_wrtPV,std::vector<float> *trk_phi_wrtPV,
				       std::vector<int> *trk_nPixHits,std::vector<int> *trk_nSCTHits){
  float mpt_phi=MPT_Loose(trk_pt,trk_eta,trk_d0_wrtPV,trk_z0_wrtPV,trk_theta_wrtPV,trk_phi_wrtPV,trk_nPixHits,trk_nSCTHits).Phi();
  return DeltaPhi(met_phi,mpt_phi);
}

float MET_Cleaning_Utils::JetTrackJetphi(TLorentzVector &jet,std::vector<TLorentzVector> &trackjet,
					 std::vector<TLorentzVector> &jettruth,float smear, bool isMC){
  if(jet.Pt()<20*GeV) return -1000;
  unsigned int tjidx= Jetmatching(jet,trackjet,0.4);
  if(tjidx<0 || tjidx>=trackjet.size())return -999;
    float jet_phi=jet.Phi();
  if(smear>0 && isMC){
    unsigned int mcidx=Jetmatching(jet,jettruth,0.4);
    if(mcidx<0 || mcidx>=jettruth.size())return -998;
    float dphi=DeltaPhi(jet_phi,jettruth.at(mcidx).Phi());
    jet_phi+=dphi*(1+smear);
    float jetphi=jet_phi;
    if(jet_phi>TMath::Pi())jet_phi=2*TMath::Pi()-jetphi;
    if(jet_phi<-TMath::Pi())jet_phi=-2*TMath::Pi()+jetphi;
  }
  return DeltaPhi(jet_phi,trackjet.at(tjidx).Phi());
}

float MET_Cleaning_Utils::dphi_MET_MPT(float met_phi,float met_truthphi, float smear, bool isMC,
				       std::vector<float> *trk_pt,std::vector<float> *trk_eta,
				       std::vector<float> *trk_d0_wrtPV,std::vector<float> *trk_z0_wrtPV,
				       std::vector<float> *trk_theta_wrtPV,std::vector<float> *trk_phi_wrtPV,
				       std::vector<int> *trk_nPixHits,std::vector<int> *trk_nSCTHits){
  if(smear>0 && isMC){
    float dphi= DeltaPhi(met_phi,met_truthphi);
    met_phi+=dphi*(1+smear);
    float metphi=met_phi;
    if(met_phi>TMath::Pi())met_phi=2*TMath::Pi()-metphi;
    if(met_phi<-TMath::Pi())met_phi=-2*TMath::Pi()+metphi;
  }
  return dphi_MET_MPT(met_phi,trk_pt,trk_eta,trk_d0_wrtPV,trk_z0_wrtPV,trk_theta_wrtPV,trk_phi_wrtPV,trk_nPixHits,trk_nSCTHits);
}


TLorentzVector MET_Cleaning_Utils::MPT_Loose(std::vector<float> *trk_pt,std::vector<float> *trk_eta,std::vector<float> *trk_d0_wrtPV,
					     std::vector<float> *trk_z0_wrtPV,std::vector<float> *trk_theta_wrtPV,
					     std::vector<float> *trk_phi_wrtPV,std::vector<int> *trk_nPixHits,
					     std::vector<int> *trk_nSCTHits ){
  TLorentzVector P(0,0,0,0);
  for(unsigned int i=0; i<trk_pt->size();i++){
    if(trk_pt->at(i)>500 &&  /*fabs(trk_eta->at(i))<2.5 &&*/ fabs(trk_d0_wrtPV->at(i))<2.0 &&
       fabs(trk_z0_wrtPV->at(i)*sin(trk_theta_wrtPV->at(i)))<10 && trk_nPixHits->at(i)+ trk_nSCTHits->at(i)>=7 ){
      TLorentzVector LT(0,0,0,0);
      LT.SetPtEtaPhiM(trk_pt->at(i),0,trk_phi_wrtPV->at(i),0);
      P-=LT;
    }
  }
  return P;
}



bool MET_Cleaning_Utils::iswithinDeadFEB(float eta, float phi, float objCone, float objIso,bool dosquare){
  for(unsigned int i=0;i<eta_min.size();i++){
    //std::cout << "iswithinDeadFEB" << std::endl;
    // compute if the jet cone overlaps with a rectangular dead feb region
    //step 1 compute rectangular cut
    float etamid=(eta_max.at(i)+eta_min.at(i))/2;
    float phimid=(phi_max.at(i)+phi_min.at(i))/2;
    float etadiff=fabs(eta_max.at(i)-eta_min.at(i))/2+objIso;
    float phidiff=fabs(DeltaPhi(phi_max.at(i),phi_min.at(i)))/2+objIso;
    float dphi=fabs(DeltaPhi(phi,phimid));
    float deta=fabs(eta-etamid);
    if(fabs(etadiff+objCone)<fabs(deta) || fabs(phidiff+objCone)<fabs(dphi))break;
    if(dosquare) return true;
    // step to now look at corners to take into accont that you have a square and a circle
    if(fabs(etadiff+objCone)>fabs(deta) && fabs(phidiff)>fabs(dphi)) return true;
    if(fabs(etadiff)>deta && fabs(phidiff+objCone)>fabs(dphi)) return true;
    if(sqrt(pow(deta-etadiff,2.0)+pow(dphi-phidiff,2.0)) < objCone)return true;
  }
  return false;
}


bool MET_Cleaning_Utils::isJetAffected(std::vector<TLorentzVector> &jet, float iso, float JCone, float minpt ){
  for(unsigned int i=0;i<jet.size();i++){
    if(isJetAffected(jet.at(i),iso,JCone,minpt))return true;
  }
  return false;
}

bool MET_Cleaning_Utils::isJetAffected(TLorentzVector &jet, float iso, float JCone, float minpt ){
  if(jet.Pt()>minpt){
    if(iswithinDeadFEB(jet.Eta(),jet.Phi(),JCone,iso)){
      return true;
    }
  }
  return false;
}


bool MET_Cleaning_Utils::hasUnmatchTrackJet(std::vector<TLorentzVector> &jet,float JCone,std::vector<TLorentzVector> &trackjet,
					    float TJCone,float minpt){
  for(unsigned int i=0;i<trackjet.size();i++){
    float dr_min=10;
    if(isJetAffected(trackjet.at(i),0.0,(JCone+TJCone)/4,minpt)){
      for(unsigned int j=0; j<jet.size();j++){
        float currentdr=dr(trackjet.at(i).Phi(),trackjet.at(i).Eta(),
                           jet.at(j).Phi(),jet.at(j).Eta());
        if(currentdr<dr_min && currentdr<JCone+TJCone){
          dr_min=currentdr;
        }
      }
      if((JCone+TJCone)/4<dr_min) return true;
    }
  }
  return false;
}

int MET_Cleaning_Utils::Jetmatching(TLorentzVector &jet,std::vector<TLorentzVector> &jetColl2,
					     float Cone){
  unsigned int matchidx=-1;
  for(unsigned int i=0;i<jetColl2.size();i++){
    float dr_min=10;
    float currentdr=dr(jetColl2.at(i).Phi(),jetColl2.at(i).Eta(),jet.Phi(),jet.Eta());
    if(currentdr<dr_min && currentdr<Cone){
      dr_min=currentdr;
      matchidx=i;
    }
  }
  return matchidx;
}



float MET_Cleaning_Utils::JetTrackJetRatio(std::vector<TLorentzVector> &jet,TLorentzVector &trackjet){
  float JCone=0.4;
  float TJCone=0.4;
  float dr_min=999;
  int jidx=-1;
  for(unsigned int j=0; j<jet.size();j++){
    float currentdr=dr(trackjet.Phi(),trackjet.Eta(),
		       jet.at(j).Phi(),jet.at(j).Eta());
    if(currentdr<dr_min && currentdr<(JCone+TJCone)/4.0){
      dr_min=currentdr;
      jidx=(int)j;
    }
  }
  if(jidx!=-1) return jet.at(jidx).Pt()/trackjet.Pt();
  return -1.0;
}

bool MET_Cleaning_Utils::isJetSafe(std::vector<TLorentzVector> &jet, std::vector<float> *BCH_CORR_JET,
                                   std::vector<float> *BCH_CORR_CELL, bool isdata, float shift,int RunNumber){
  if(HasLArFailure(RunNumber) || !isdata)return isJetSafe(jet,BCH_CORR_JET,BCH_CORR_CELL,isdata,shift);
  return true;
}

bool MET_Cleaning_Utils::HasLArFailure(int RunNumber){
  return (180614<=RunNumber && RunNumber<185353);
}


bool MET_Cleaning_Utils::isJetSafe(std::vector<TLorentzVector> &jet, std::vector<float> *BCH_CORR_JET, 
				   std::vector<float> *BCH_CORR_CELL, bool isdata, float shift){
  if(jet.size()!=BCH_CORR_CELL->size() || jet.size()!=BCH_CORR_JET->size()){
    std::cout << "FATAL: jets BCH_CORR_CELL and BCH_CORR_JET do not match! FIX YOU INPUTS  " << std::endl;
  }
  for(unsigned int i=0;i<jet.size();i++){
    float THRESHOLD=jet_THRESHOLD+jet_THRESHOLD_1sigma*shift;
    if(isdata)THRESHOLD*=(1-BCH_CORR_JET->at(i))/(1-BCH_CORR_CELL->at(i));
    if(isJetAffected(jet.at(i),0.0,jet_iso,THRESHOLD))return false;
  }
  return true;
}


float MET_Cleaning_Utils::dphi_Jet_MET(std::vector<TLorentzVector> &jet, float met_phi){
  float dphi_min=10;
  for(unsigned int i=0;i<jet.size();i++){
    float Dphi= fabs(DeltaPhi(jet.at(i).Phi(),met_phi));
    if(dphi_min>Dphi){
      dphi_min=Dphi;
    }
  }
  return dphi_min;
}
