#include "SingleTopRootAnalysis/Base/Dictionary/Top_MET.hpp"

#include <iostream>
using std::cout;
using std::endl;

double Top_MET::SoftJetEM_FractionUncert[]={0.105};
double Top_MET::SoftJetEM_SumETBins[]={1000};
unsigned int Top_MET::SoftJetEM_NBins=1;
double Top_MET::CellOutEM_FractionUncert[]={0.132};
double Top_MET::CellOutEM_SumETBins[]={1000};
unsigned int Top_MET::CellOutEM_NBins=1;

//Setup uncerianties for the CellOut
double Top_MET::SoftJetLC_FractionUncert[]={0, 0.00871549, 0.0897552, 0.0892624, 0.0841947, 0.0872207, 0.0894101, 0.090524, 0.0902687, 0.0909246, 0.0926338, 0.0917935, 0.0932977, 0.0892082, 0.093984, 0.0974575, 0.100168, 0.100101, 0.0979492, 0.100863, 0.102229, 0.102693, 0.108248, 0.10524, 0.100301, 0.104939, 0.110576, 0.103886, 0.0982754, 0.102644, 0.105031, 0.0989276, 0.101526, 0.0975329, 0.102695, 0.112649, 0.108279, 0.0992751, 0.111916, 0.107718, 0.115603, 0.115638, 0.107733, 0.111464, 0.103725, 0.102329, 0.117397, 0.112757, 0.103333, 0.0925, 0.085404, 0.114804, 0.150765, 0.117442, 0.104841, 0.112848, 0.109667, 0.116017, 0.113594, 0.107885, 0.111078, 0.112407}; 
double Top_MET::SoftJetLC_SumETBins[]={-5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0, 165.0, 170.0, 175.0, 180.0, 185.0, 190.0, 195.0, 200.0, 205.0, 210.0, 215.0, 220.0, 225.0, 230.0, 235.0, 240.0, 245.0, 250.0, 255.0, 260.0, 265.0, 270.0, 275.0, 280.0, 285.0, 290.0, 295.0, 300.0}; 
unsigned int Top_MET::SoftJetLC_NBins=62;
double Top_MET::CellOutEFlow_FractionUncert[]={0.136409, 0.153857, 0.153689, 0.149694, 0.138453, 0.139799, 0.133924, 0.139531, 0.138844, 0.134789, 0.138834, 0.132662, 0.142114, 0.13908, 0.13661, 0.128625, 0.133791, 0.125408, 0.122293, 0.133035, 0.132411, 0.131334, 0.128113, 0.13843, 0.12271, 0.129, 0.125528, 0.125005, 0.130233, 0.125374, 0.120071, 0.125658, 0.121132, 0.121722};
  double Top_MET::CellOutEFlow_SumETBins[]={-5.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 125.0, 130.0, 135.0, 140.0, 145.0, 150.0, 155.0, 160.0}; 
unsigned int Top_MET::CellOutEFlow_NBins=34;


//constructor
Top_MET::Top_MET():
  JetPtThreshold(20*GeV)
  ,GeV(1000)
  ,METComOk(false)
  ,METWeightsOK(false)
  ,pileupsigma(0.0)
{
  Reset();
}

//destructor
Top_MET::~Top_MET(){}

void Top_MET::Reset(){
  e_scaled.clear();
  e.clear();
  ph_scaled.clear();
  ph.clear();
  jet_scaled.clear();
  jet.clear();
  jet_em.clear();
  mucb_scaled.clear();
  mucb.clear();
  mums_scaled.clear();
  mums.clear();
  muid_scaled.clear();
  muid.clear();
  Jets_JESup.clear();
  SoftJet_sigma=0.0;
  CellOut_sigma=0.0;
  pileupsigma=0.0;
}

void Top_MET::Set_Electrons(std::vector<TLorentzVector> e_scaled_, std::vector<TLorentzVector> e_){
  e_scaled=e_scaled_;
  e=e_;
}

void Top_MET::Set_Photons(std::vector<TLorentzVector> ph_scaled_, std::vector<TLorentzVector> ph_){
  ph_scaled=ph_scaled_;
  ph=ph_;
}

void Top_MET::Set_Jets(std::vector<TLorentzVector> jet_scaled_, std::vector<TLorentzVector> jet_, std::vector<TLorentzVector> jet_em_){
  jet_scaled=jet_scaled_;
  jet=jet_;
  jet_em=jet_em_;
}

void Top_MET::Set_Muons(std::vector<TLorentzVector> mucb_scaled_, std::vector<TLorentzVector> mucb_,
			std::vector<TLorentzVector> mums_scaled_, std::vector<TLorentzVector> mums_,
			std::vector<TLorentzVector> muid_scaled_, std::vector<TLorentzVector> muid_){
  mucb_scaled=mucb_scaled_;
  mucb=mucb_;
  mums_scaled=mums_scaled_;
  mums=mums_;
  muid_scaled=muid_scaled_;
  muid=muid_;
}

void Top_MET::ApplySoftJetUncertainty(float sigma, std::vector<TLorentzVector> Jets_JESup_, int softjet_sysType){
  SoftJet_sigma=sigma;
  Jets_JESup=Jets_JESup_;
  if(softjet_sysType==EM || softjet_sysType==LC) SoftJet_SysType=softjet_sysType;
  else{
    std::cout << "Top_MET::ApplyCellOutUncertainty Invalid value for SoftJet_sysType. Using EM" << std::endl;
    SoftJet_SysType=EM;
  }
}

void Top_MET::ApplyCellOutUncertainty(float sigma,int cellout_sysType){
  CellOut_sigma=sigma;
  if(cellout_sysType==EM || cellout_sysType==EFlow) CellOut_SysType=cellout_sysType;
  else{
    std::cout << "Top_MET::ApplyCellOutUncertainty Invalid value for CellOut_sysType. Using EM" << std::endl;
    CellOut_SysType=EM;
  }
}

void Top_MET::ApplyPileupUncertainty(float sigma){
  pileupsigma=sigma;
}


float Top_MET::MET_EtMiss(bool correct,bool METComp,int type){
  float sumet(0),ex(0),ey(0);
  Get_MET(sumet,ex,ey,correct,METComp,type);
  return sqrt(ex*ex+ey*ey);
}

float Top_MET::MET_EyMiss(bool correct,bool METComp,int type){
  float sumet(0),ex(0),ey(0);
  Get_MET(sumet,ex,ey,correct,METComp,type);
  return ey;
}

float Top_MET::MET_ExMiss(bool correct,bool METComp,int type){
  float sumet(0),ex(0),ey(0);
  Get_MET(sumet,ex,ey,correct,METComp,type);
  return ex;
}

float Top_MET::MET_SumEt(bool correct,bool METComp,int type){
  float sumet(0),ex(0),ey(0);
  Get_MET(sumet,ex,ey,correct,METComp,type);
  return sumet;
}

float Top_MET::MET_MetPhi(bool correct,bool METComp,int type){
  float sumet(0),ex(0),ey(0);
  Get_MET(sumet,ex,ey,correct,METComp,type);
  return TMath::ATan2(ey,ex);
}

void Top_MET::Get_MET(float &met_et, float &met_sumet, float &met_etx, float &met_ety, float &met_phi, 
		      bool correct, bool METComp, int type){
  Get_MET(met_sumet,met_etx,met_ety,correct,METComp,type);
  met_phi=TMath::ATan2(met_ety,met_etx);
  met_et=sqrt(met_etx*met_etx+met_ety*met_ety);
}

void Top_MET::Get_MET(float &met_sumet, float &met_etx, float &met_ety, 
		      bool correct, bool METComp, int type){
  //set Default for MET composition or Error.
  met_sumet=0;
  met_etx=0;
  met_ety=0;

  if(!METComOk){
    std::cout << "Error Met composition has not been set (Top_MET::Get_MET)" << std::endl;
    return;
  }
  if(!METWeightsOK){
    std::cout << "Error Met Weights have not been set (Top_MET::Get_MET)" << std::endl;
    return;
  }
  if(METComp){
    met_sumet=0;
    met_etx=0;
    met_ety=0;
  }
  else if(type==MET_All){
    met_sumet=(*MET_sumet);
    met_etx=(*MET_etx);
    met_ety=(*MET_ety);
  }
  else if(type==MET_e){
    met_sumet=(*MET_RefEle_sumet);
    met_etx=(*MET_RefEle_etx);
    met_ety=(*MET_RefEle_ety);
  }
  else if(type==MET_ph){
    met_sumet=(*MET_RefGamma_sumet);
    met_etx=(*MET_RefGamma_etx);
    met_ety=(*MET_RefGamma_ety);
  }
  else if(type==MET_jet){
    met_sumet=(*MET_RefJet_sumet);
    met_etx=(*MET_RefJet_etx);
    met_ety=(*MET_RefJet_ety);
  }
  else if(type==MET_softjet){
    met_sumet=(*MET_SoftJets_sumet);
    met_etx=(*MET_SoftJets_etx);
    met_ety=(*MET_SoftJets_ety);
  }
  else if(type==MET_mu){
    met_sumet=(*MET_Muon_Total_Muid_sumet);
    met_etx=(*MET_Muon_Total_Muid_etx);
    met_ety=(*MET_Muon_Total_Muid_ety);
  }
  else if(type==MET_cellout){
    met_sumet=(*MET_CellOut_sumet);
    met_etx=(*MET_CellOut_etx);
    met_ety=(*MET_CellOut_ety);
  }
  else{
    std::cout << "Error incorrect MET Type selected: " << type << " --- using default settings: results maybe invalid" << std::endl;
    Get_MET(met_sumet,met_etx,met_ety);
    return;
  }
  if(!correct && !METComp){
    return;
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // Start Composition
  //
  // Modify MET for scaled Electrons
  if(type==MET_All || type==MET_e){
    for(unsigned int i=0; i<e.size();i++){
      float obj_pt=e.at(i).Pt();
      float obj_px=e.at(i).Px();
      float obj_py=e.at(i).Py();
      float obj_wx=el_MET_wpx->at(i).at(0);
      float obj_wy=el_MET_wpy->at(i).at(0);
      float obj_wt=el_MET_wet->at(i).at(0);
      float obj_scale=1;
      if(e.at(i).Pt()!=0 ) obj_scale=e_scaled.at(i).Pt()/e.at(i).Pt();
      CorrectMet(met_sumet,met_etx,met_ety,obj_pt,obj_px,obj_py,obj_wt,obj_wx,obj_wy,obj_scale,!METComp);
    }
  }
  //Modify MET for scaled Photons
  if(type==MET_All || type==MET_ph){
    for(unsigned int i=0; i<ph.size();i++){
      float obj_pt=ph.at(i).Pt();
      float obj_px=ph.at(i).Px();
      float obj_py=ph.at(i).Py();
      float obj_wx=ph_MET_wpx->at(i).at(0);
      float obj_wy=ph_MET_wpy->at(i).at(0);
      float obj_wt=ph_MET_wet->at(i).at(0);
      float obj_scale=1;
      if(ph_scaled.at(i).Pt()!=0 && ph.at(i).Pt()!=0 ) obj_scale=ph_scaled.at(i).Pt()/ph.at(i).Pt();
      CorrectMet(met_sumet,met_etx,met_ety,obj_pt,obj_px,obj_py,obj_wt,obj_wx,obj_wy,obj_scale,!METComp);
    }
  }
  //Modify MET for scaled Jets
  if(type==MET_All || type==MET_jet || type==MET_softjet){
    for(unsigned int i=0; i<jet.size();i++){
      float obj_pt=jet.at(i).Pt();
      float obj_px=jet.at(i).Px();
      float obj_py=jet.at(i).Py();
      float obj_wx=jet_wpx->at(i).at(0);
      float obj_wy=jet_wpy->at(i).at(0);
      float obj_wt=jet_wet->at(i).at(0);
      float obj_scale=1;
      if(jet_scaled.at(i).Pt()!=0 && jet.at(i).Pt()!=0){
	obj_scale=jet_scaled.at(i).Pt()/jet.at(i).Pt();
	if(jet_scaled.at(i).Pt()>20*GeV && jet.at(i).Pt()<20*GeV)obj_scale*=jet.at(i).Pt()/jet_em.at(i).Pt();
	if(jet_scaled.at(i).Pt()<20*GeV && jet.at(i).Pt()>20*GeV)obj_scale*=jet_em.at(i).Pt()/jet.at(i).Pt();
      }
      if(jet_scaled.at(i).Pt()<20*GeV && jet.at(i).Pt()<20*GeV)obj_scale=1; //do not scale soft jets
      // remove outof range jets when only looking at the JetRef or the SoftJet term 
      if(type==MET_softjet && ((jet_scaled.at(i).Pt()>=20*GeV && correct) ||  (!correct &&  jet.at(i).Pt()>=20*GeV )))  obj_scale=0;
      if(type==MET_jet     && ((jet_scaled.at(i).Pt()<20*GeV  && correct) ||  (!correct && jet.at(i).Pt()<20*GeV)))  obj_scale=0;
      // allow for efficiency systematic uncertainty and put it in the SoftJet term
      if(jet_scaled.at(i).Pt()==0.0 && jet.at(i).Pt()!=0 && correct){
	obj_scale=1; // Soft jet term is at EM scale
	if(jet.at(i).Pt()>20*GeV) obj_scale=jet_em.at(i).Pt()/jet.at(i).Pt();
	if(type==MET_jet) obj_scale=0;
      } 
      CorrectMet(met_sumet,met_etx,met_ety,obj_pt,obj_px,obj_py,obj_wt,obj_wx,obj_wy,obj_scale,!METComp);
    }
  }
  //Modify MET for scaled muons
  if(type==MET_All || type==MET_mu){
    for(unsigned int i=0; i<mucb.size();i++){
      float obj_pt=mucb.at(i).Pt();
      float obj_px=mucb.at(i).Px();
      float obj_py=mucb.at(i).Py();
      float obj_scale=1;
      if(mucb.at(i).Pt()!=0) obj_scale=mucb_scaled.at(i).Pt()/mucb.at(i).Pt();
      if ((mu_muid_MET_statusWord->at(i).at(0)) &DEFAULT) {
	//already setup as default
      } 
      else if (mu_muid_MET_statusWord->at(i).at(0) &SPECTRO){
	// use the muon spectrometer track px and py
        obj_pt=mums.at(i).Pt();
	obj_px=mums.at(i).Px();
	obj_py=mums.at(i).Py();
	if(mums.at(i).Pt()!=0)obj_scale=mums_scaled.at(i).Pt()/mums.at(i).Pt();
      } 
      else if (mu_muid_MET_statusWord->at(i).at(0) &TRACK){
	// use the muon track px and py
	obj_pt=muid.at(i).Pt();
	obj_px=muid.at(i).Px();
	obj_py=muid.at(i).Py();
	if(muid.at(i).Pt()!=0) obj_scale=muid_scaled.at(i).Pt()/muid.at(i).Pt();
      }
      float obj_wx=mu_muid_MET_wpx->at(i).at(0);
      float obj_wy=mu_muid_MET_wpy->at(i).at(0);
      float obj_wt=mu_muid_MET_wet->at(i).at(0); // muons do not contribute to the sumET
      if(type!=MET_mu) obj_wt=0;
      CorrectMet(met_sumet,met_etx,met_ety,obj_pt,obj_px,obj_py,obj_wt,obj_wx,obj_wy,obj_scale,!METComp);
    }
  }
  // Add the cell Out for METComp
  if((type==MET_All || type==MET_cellout) && METComp){
    met_sumet+=(*MET_CellOut_sumet);
    met_etx+=(*MET_CellOut_etx);
    met_ety+=(*MET_CellOut_ety);
  }
  //SoftJet Uncertainty
  if(SoftJet_sigma!=0.0){
    float obj_pt=(*MET_SoftJets_sumet);
    float obj_px=(*MET_SoftJets_etx);
    float obj_py=(*MET_SoftJets_ety);
    float obj_wx=1;
    float obj_wy=1;
    float obj_wt=1;
    float obj_scale=1;
    float SoftJetSys=0;
    if(SoftJet_SysType==EM){
      SoftJetSys=SoftJetEM_FractionUncert[0];
      for(unsigned int j=0;j+1<SoftJetEM_NBins;j++){
	SoftJetSys=SoftJetEM_FractionUncert[j];
	if(obj_pt/GeV<SoftJetEM_SumETBins[j+1]) break;
      }
    }
    else if(SoftJet_SysType==LC){
      SoftJetSys=SoftJetLC_FractionUncert[0];
      for(unsigned int j=0;j+1<SoftJetLC_NBins;j++){
	SoftJetSys=SoftJetLC_FractionUncert[j];
	if(obj_pt/GeV<SoftJetLC_SumETBins[j+1]) break;
      }
    }
    obj_scale=1+SoftJet_sigma*SoftJetSys;
    CorrectMet(met_sumet,met_etx,met_ety,obj_pt,obj_px,obj_py,obj_wt,obj_wx,obj_wy,obj_scale,true);
  }
  // CellOut uncertainty
  if(CellOut_sigma!=0.0){
    float obj_pt=(*MET_CellOut_sumet);
    float obj_px=(*MET_CellOut_etx);
    float obj_py=(*MET_CellOut_ety);
    float obj_wx=1;
    float obj_wy=1;
    float obj_wt=1;
    float obj_scale=1;
    float CellOutSys=1;
    if(CellOut_SysType==EM){
      CellOutSys=CellOutEM_FractionUncert[0];
      for(unsigned int j=0;j+1<CellOutEM_NBins;j++){
        CellOutSys=CellOutEM_FractionUncert[j];
        if(obj_pt/GeV<CellOutEM_SumETBins[j+1]) break;
      }
      obj_scale=1+CellOut_sigma*CellOutSys;
      CorrectMet(met_sumet,met_etx,met_ety,obj_pt,obj_px,obj_py,obj_wt,obj_wx,obj_wy,obj_scale,true);
    }
    else if(CellOut_SysType==EFlow){
      CellOutSys=CellOutEFlow_FractionUncert[0];
      for(unsigned int j=0;j+1<CellOutEFlow_NBins;j++){
        CellOutSys=CellOutEFlow_FractionUncert[j];
        if(obj_pt/GeV<CellOutEFlow_SumETBins[j+1]) break;
      }
      obj_scale=1+CellOut_sigma*CellOutSys;
      CorrectMet(met_sumet,met_etx,met_ety,obj_pt,obj_px,obj_py,obj_wt,obj_wx,obj_wy,obj_scale,true);
    }
  }
  if(pileupsigma!=0.0){
    float obj_pt=(*MET_SoftJets_sumet);
    float obj_px=(*MET_SoftJets_etx);
    float obj_py=(*MET_SoftJets_ety);
    float obj_wx=1;
    float obj_wy=1;
    float obj_wt=1;
    float obj_scale=1;
    obj_scale=1+pileupsigma*0.1;
    CorrectMet(met_sumet,met_etx,met_ety,obj_pt,obj_px,obj_py,obj_wt,obj_wx,obj_wy,obj_scale,true);

    //cellout
    obj_pt=(*MET_CellOut_sumet);
    obj_px=(*MET_CellOut_etx);
    obj_py=(*MET_CellOut_ety);
    obj_scale=1+pileupsigma*0.1;
    CorrectMet(met_sumet,met_etx,met_ety,obj_pt,obj_px,obj_py,obj_wt,obj_wx,obj_wy,obj_scale,true);

  }
}

void Top_MET::CheckConsistency(){
  for(unsigned int i=MET_All; i<=MET_cellout;i++){
    if((MET_EtMiss(false,false,i)-MET_EtMiss(false,true,i))/MET_EtMiss(false,false,i)>1E-4){
      std::cout << "Top_MET::CheckConsistency Error MET Et for enum MET_Types=" << i << " not consistent at :"  
		<<(MET_EtMiss(false,false,i)-MET_EtMiss(false,true,i)) << " D3PD value: " 
		<< MET_EtMiss(false,false,i) << " METComposition value " 
		<< MET_EtMiss(false,true,i)  << " Please check the variables corresponding to the MET Type " << i << std::endl; 
    }
    if((MET_ExMiss(false,false,i)-MET_ExMiss(false,true,i))/MET_EtMiss(false,false,i)>1E-4){
      std::cout << "Top_MET::CheckConsistency Error MET Ex  for enum MET_Types=" << i << " not consistent "
                <<(MET_ExMiss(false,false,i)-MET_ExMiss(false,true,i)) << " D3PD value: "
                << MET_ExMiss(false,false,i) << " METComposition value "
                << MET_ExMiss(false,true,i)  << " Please check the variables corresponding to the MET Type " << i << std::endl;
    }
    if((MET_EyMiss(false,false,i)-MET_EyMiss(false,true,i))/MET_EtMiss(false,false,i)>1E-4){
      std::cout << "Top_MET::CheckConsistency Error MET Ey  for enum MET_Types=" << i << " not consistent "
                <<(MET_EyMiss(false,false,i)-MET_EyMiss(false,true,i)) << " D3PD value: "
                << MET_EyMiss(false,false,i) << " METComposition value "
                << MET_EyMiss(false,true,i)  << " Please check the variables corresponding to the MET Type " << i << std::endl;
    }
    if((MET_SumEt(false,false,i)-MET_SumEt(false,true,i))/MET_SumEt(false,false,i)>1E-4){
      std::cout << "Top_MET::CheckConsistency Error MET SumEt  for enum MET_Types=" << i << " not consistent "
                <<(MET_SumEt(false,false,i)-MET_SumEt(false,true,i)) << " D3PD value: "
                << MET_SumEt(false,false,i) << " METComposition value "
                << MET_SumEt(false,true,i)  << " Please check the variables corresponding to the MET Type " << i << std::endl;
    }
  }
}



void Top_MET::CorrectMet(float &met_sumet, float &met_etx, float &met_ety, float obj_pt, float obj_px, float obj_py,
			 float obj_wt,  float obj_wx,   float obj_wy,   float obj_scale,bool remove){
  if(remove){
    met_etx+= obj_px*obj_wx;
    met_ety+= obj_py*obj_wy;
    met_sumet-= obj_pt*obj_wt;
  }
  met_etx-= obj_px*obj_wx*obj_scale;
  met_ety-= obj_py*obj_wy*obj_scale;
  met_sumet += obj_pt*obj_wt*obj_scale;
}

void Top_MET::Set_METComposition(  float &MET_etx_,
				   float &MET_ety_,
				   float &MET_sumet_,
				   float &MET_RefEle_sumet_,
				   float &MET_RefEle_etx_,
				   float &MET_RefEle_ety_,
				   float &MET_RefGamma_sumet_,
				   float &MET_RefGamma_etx_,
				   float &MET_RefGamma_ety_,
				   float &MET_RefJet_sumet_,
				   float &MET_RefJet_etx_,
				   float &MET_RefJet_ety_,
				   float &MET_SoftJets_sumet_,
				   float &MET_SoftJets_etx_,
				   float &MET_SoftJets_ety_,
				   float &MET_CellOut_sumet_,
				   float &MET_CellOut_etx_,
				   float &MET_CellOut_ety_,
				   float &MET_Muon_Total_Muid_sumet_,
				   float &MET_Muon_Total_Muid_etx_,
				   float &MET_Muon_Total_Muid_ety_
				   ){
  MET_etx=&MET_etx_;
  MET_ety=&MET_ety_;
  MET_sumet=&MET_sumet_;
  MET_RefEle_sumet=&MET_RefEle_sumet_;
  MET_RefEle_etx=&MET_RefEle_etx_;
  MET_RefEle_ety=&MET_RefEle_ety_;
  MET_RefGamma_sumet=&MET_RefGamma_sumet_;
  MET_RefGamma_etx=&MET_RefGamma_etx_;
  MET_RefGamma_ety=&MET_RefGamma_ety_;
  MET_RefJet_sumet=&MET_RefJet_sumet_;
  MET_RefJet_etx=&MET_RefJet_etx_;
  MET_RefJet_ety=&MET_RefJet_ety_;
  MET_SoftJets_sumet=&MET_SoftJets_sumet_;
  MET_SoftJets_etx=&MET_SoftJets_etx_;
  MET_SoftJets_ety=&MET_SoftJets_ety_;
  MET_CellOut_sumet=&MET_CellOut_sumet_;
  MET_CellOut_etx=&MET_CellOut_etx_;
  MET_CellOut_ety=&MET_CellOut_ety_;
  MET_Muon_Total_Muid_sumet=&MET_Muon_Total_Muid_sumet_;
  MET_Muon_Total_Muid_etx=&MET_Muon_Total_Muid_etx_;
  MET_Muon_Total_Muid_ety=&MET_Muon_Total_Muid_ety_;  
  METComOk=true;
}



void Top_MET::Set_METWeights( std::vector<std::vector<float> >*  el_MET_wpx_,
			      std::vector<std::vector<float> >*  el_MET_wpy_,
			      std::vector<std::vector<float> >*  el_MET_wet_,
			      std::vector<std::vector<float> >*  ph_MET_wpx_,
			      std::vector<std::vector<float> >*  ph_MET_wpy_,
			      std::vector<std::vector<float> >*  ph_MET_wet_,
			      std::vector<std::vector<float> >*  jet_wpx_,
			      std::vector<std::vector<float> >*  jet_wpy_,
			      std::vector<std::vector<float> >*  jet_wet_,
			      std::vector<std::vector<unsigned int> >*  mu_muid_MET_statusWord_,
			      std::vector<std::vector<float> >*  mu_muid_MET_wpx_,
			      std::vector<std::vector<float> >*  mu_muid_MET_wpy_,
			      std::vector<std::vector<float> >*  mu_muid_MET_wet_
			      ){
  el_MET_wpx=el_MET_wpx_;
  el_MET_wpy=el_MET_wpy_;
  el_MET_wet=el_MET_wet_;
  ph_MET_wpx=ph_MET_wpx_;
  ph_MET_wpy=ph_MET_wpy_;
  ph_MET_wet=ph_MET_wet_;
  jet_wpx=jet_wpx_;
  jet_wpy=jet_wpy_;
  jet_wet=jet_wet_;
  mu_muid_MET_statusWord=mu_muid_MET_statusWord_;
  mu_muid_MET_wpx=mu_muid_MET_wpx_;
  mu_muid_MET_wpy=mu_muid_MET_wpy_;
  mu_muid_MET_wet=mu_muid_MET_wet_;
  METWeightsOK=true;
}


void Top_MET::Get_ElectronWeights(unsigned int i,float  &wet, float  &wpx, float  &wpy){
  wet=0;
  wpx=0;
  wpy=0;
  if(el_MET_wpx->size()>i && i>=0){
    wet=el_MET_wet->at(i).at(0);
    wpx=el_MET_wpx->at(i).at(0);
    wpy=el_MET_wpy->at(i).at(0);
  }
}

void Top_MET::Get_PhotonWeights(unsigned int i,float  &wet, float  &wpx, float  &wpy){
  wet=0;
  wpx=0;
  wpy=0;
  if(ph_MET_wpx->size()>i && i>=0){
    wet=ph_MET_wet->at(i).at(0);
    wpx=ph_MET_wpx->at(i).at(0);
    wpy=ph_MET_wpy->at(i).at(0);
  }
}

void Top_MET::Get_JetWeights(unsigned int i,float  &wet, float  &wpx, float  &wpy){
  wet=0;
  wpx=0;
  wpy=0;
  if(jet_wpx->size()>i && i>=0){
    wet=jet_wet->at(i).at(0);
    wpx=jet_wpx->at(i).at(0);
    wpy=jet_wpy->at(i).at(0);
  }
}
void Top_MET::Get_MuonWeights(unsigned int i,float  &wet, float  &wpx, float  &wpy,unsigned int &mu_statusWord){
  wet=0;
  wpx=0;
  wpy=0;
  mu_statusWord=0;
  if(mu_muid_MET_wpx->size()>i && i>=0){
    wet=mu_muid_MET_wet->at(i).at(0);
    wpx=mu_muid_MET_wpx->at(i).at(0);
    wpy=mu_muid_MET_wpy->at(i).at(0);
    mu_statusWord=mu_muid_MET_statusWord->at(i).at(0);
  }
}

