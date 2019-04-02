/******************************************************************************
 * testVar.cpp                                                                *
 *
 * Testing out some variables being added to the skim tree.
 *
 * History
 *      10th Nov 2016 - Created by D. Leggat. Just a test, so probably not that important
 ******************************************************************************/

#include "SingleTopRootAnalysis/Vars/bstarVars.hpp"
#include <iostream>

//Test out a couple of variables, one int and one float I guess
bstarVars::bstarVars(bool makeHistos){
  SetName("bstarVars");
  //some Vars will use in bstar analyize
  _floatVars["M_Mass_bstar"] = 500.;
  _floatVars["M_Mass_top"] = 500.;
  _floatVars["M_Mass_W"]       = 300.;
  _floatVars["M_Mass_Lepton"]       = 300.;
  _floatVars["M_Mass_Jet"]       = 300.;
  _floatVars["M_Mass_BJet"]       = 300.;
  _floatVars["M_Mass_BoostedJet"]       = 300.;
  _floatVars["M_prunedmass_BoostedJet"]       = 0.;
  _floatVars["M_tau21_BoostedJet"]       = 1.0;
  _floatVars["M_tau32_BoostedJet"]       = 1.0;
  _floatVars["M_E_W"]          = 500.;
  _floatVars["M_Pt_W"]         = 500.;
  _floatVars["M_Eta_W"]        = -6.;
  _floatVars["M_Phi_W"]        = -3.2;
  _floatVars["M_E_Top"]          =500.; 
  _floatVars["M_Pt_Top"]         = 500.;
  _floatVars["M_Eta_Top"]        = -6.;
  _floatVars["M_Phi_Top"]        = -3.2;
  _floatVars["M_E_bstar"]          =500.;
  _floatVars["M_Pt_bstar"]         = 500.;
  _floatVars["M_Eta_bstar"]        = -6.;
  _floatVars["M_Phi_bstar"]        = -3.2;
  _floatVars["M_E_BoostedJet"]          =500.;
  _floatVars["M_Pt_BoostedJet"]         = 500.;
  _floatVars["M_Eta_BoostedJet"]        = -6.;
  _floatVars["M_Phi_BoostedJet"]        = -3.2;
  _floatVars["M_E_Jet"]   =500.;
  _floatVars["M_Pt_Jet"]  = 500.;
  _floatVars["M_Eta_Jet"] = -6.;
  _floatVars["M_Phi_Jet"] = -3.2; 
  _floatVars["M_E_BJet"]   =500.;
  _floatVars["M_Pt_BJet"]  = 500.;
  _floatVars["M_Eta_BJet"] = -6.;
  _floatVars["M_Phi_BJet"] = -3.2; 
   _floatVars["M_E_Lepton"]   =500.;
  _floatVars["M_Pt_Lepton"]  = 500.;
  _floatVars["M_Eta_Lepton"] = -6.;
  _floatVars["M_Phi_Lepton"] = -3.2; 


 //others
  
  SetDoHists(makeHistos);

}

void bstarVars::FillBranches(EventContainer * evtObj){

  //Initialise vectors;
//  BoostedJet.clear();
  Jet.clear();
  BJet.clear();

  sumJet.SetPtEtaPhiE(0,0,0,0);
  sumBJet.SetPtEtaPhiE(0,0,0,0);
  totalJets.SetPtEtaPhiE(0,0,0,0);
  
  Int_t _bstar = 0;
  TEnv *config = evtObj-> GetConfig();
  _bstar = config -> GetValue("ifbstar",0);
  EventTree* eventTree = evtObj->GetEventTree();
  if(_bstar){
	  //set up lepton and met
	  TLorentzVector Lepton(00,0,0,0);		
	  TLorentzVector lep_p4(0,0,0,0);		
	  for(int i=0;i<evtObj->muonsToUsePtr->size();i++){
		  Lepton.SetPtEtaPhiE(evtObj->muonsToUsePtr->at(i).Pt(),evtObj->muonsToUsePtr->at(i).Eta(),evtObj->muonsToUsePtr->at(i).Phi(),evtObj->muonsToUsePtr->at(i).E());
		  lep_p4 = TLorentzVector(evtObj->muonsToUsePtr->at(i).Px(),evtObj->muonsToUsePtr->at(i).Py(),evtObj->muonsToUsePtr->at(i).Pz(),evtObj->muonsToUsePtr->at(i).E()); 

		  _floatVars["M_Mass_Lepton"]       = (Lepton).M();
		  _floatVars["M_E_Lepton"]          = (Lepton).E();
		  _floatVars["M_Pt_Lepton"]         = (Lepton).Pt();
		  _floatVars["M_Eta_Lepton"]        = (Lepton).Eta();
		  _floatVars["M_Phi_Lepton"]        = (Lepton).Phi();
	  }

	  TLorentzVector met(00,0,0,0);
	  TLorentzVector W(00,0,0,0);
	  TLorentzVector lep_neutrino(0,0,0,0);

	  met.SetPtEtaPhiE(evtObj->missingEt,0,evtObj->missingPhi,0.);
	  //Get W
	  double M_W  = 80.4;
	  double M_mu =  0.10566;
	  double M_e = 0.511e-3;
	  double M_lepton = M_mu;
	  double emu = Lepton.E();
	  double pxmu = Lepton.Px();
	  double pymu = Lepton.Py();
	  double pzmu = Lepton.Pz();
	  double pxnu = met.Px();
	  double pynu = met.Py();
	  double pznu = 0.;
	  double a = M_W*M_W - M_lepton*M_lepton + 2.0*(pxmu*pxnu + pymu*pynu);
	  double A = 4.0*(emu*emu - pzmu*pzmu);
	  double B = -4.0*a*pzmu;
	  double C = 4.0*emu*emu*(pxnu*pxnu + pynu*pynu) - a*a;

	  double tmproot = B*B - 4.0*A*C;
	  bool isComplex_ = false;
	  if (tmproot<0) {
		  isComplex_= true;
		  pznu = - B/(2*A); // take real part of complex roots
	  }
	  else {
		  isComplex_ = false;
		  double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0*A);
		  double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0*A);
		  if (TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2;}
		  else pznu = tmpsol1;
		  // if pznu is > 300 pick the most central root
		  if ( pznu > 300. ) {
			  if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) pznu = tmpsol1;
			  else pznu = tmpsol2;
		  }
	  } 

	  float mu_nupz=pznu ;
	  float  mu_nuE  = sqrt( met.Px()*met.Px() + met.Py()*met.Py()+mu_nupz*mu_nupz); 
	  lep_neutrino  = TLorentzVector(met.Px(), met.Py(), mu_nupz, mu_nuE );
	  //W  = lep_neutrino+Lepton;
	  W  = lep_neutrino+lep_p4;
	  _floatVars["M_Mass_W"]       = (W).M();
	  _floatVars["M_E_W"]          = (W).E();
	  _floatVars["M_Pt_W"]         = (W).Pt();
	  _floatVars["M_Eta_W"]        = (W).Eta();
	  _floatVars["M_Phi_W"]        = (W).Phi();

	  //Get Top
	  TLorentzVector Jet_p4(0,0,0,0);
	  TLorentzVector Top_p4(0,0,0,0);
	  TLorentzVector Top(0,0,0,0);

	  for(int i=0;i<evtObj->taggedjetsToUsePtr->size();i++){
		  Jet_p4 = TLorentzVector(evtObj->taggedjetsToUsePtr->at(i).Px(),evtObj->taggedjetsToUsePtr->at(i).Py(),evtObj->taggedjetsToUsePtr->at(i).Pz(),evtObj->taggedjetsToUsePtr->at(i).E());
		  //		 Jet = TLorentzVector(evtObj->jetsToUsePtr->at(i).Pt(),evtObj->jetsToUsePtr->at(i).Eta(),evtObj->jetsToUsePtr->at(i).Phi(),evtObj->jetsToUsePtr->at(i).E());
		  Top = Jet_p4+W;
	  	  Top_p4.SetPtEtaPhiM(Top.Pt(),Top.Eta(),Top.Phi(),Top.M());


		  _floatVars["M_Mass_top"]       = (Top).M();
		  _floatVars["M_E_Top"]          = (Top).E();
		  _floatVars["M_Pt_Top"]         = (Top).Pt();
		  _floatVars["M_Eta_Top"]        = (Top).Eta();
		  _floatVars["M_Phi_Top"]        = (Top).Phi();	
		
		
		  _floatVars["M_Mass_BJet"]     = (Jet_p4).M();
		  _floatVars["M_E_BJet"]     = (Jet_p4).E();
        	  _floatVars["M_Pt_BJet"]     = (Jet_p4).Pt();
     	          _floatVars["M_Eta_BJet"]     = (Jet_p4).Eta();
         	  _floatVars["M_Phi_BJet"]     = (Jet_p4).Phi();


	  }

	  //Get Jet and BJet pt eta phi 
	  TLorentzVector tempjet(0,0,0,0);
	  for (auto jet : evtObj->taggedJets){
		  tempjet.SetPtEtaPhiE(jet.Pt(),jet.Eta(),jet.Phi(),jet.E());
		  BJet.push_back(tempjet);
	  }
	  for (auto jet : evtObj->jets){
		  tempjet.SetPtEtaPhiE(jet.Pt(),jet.Eta(),jet.Phi(),jet.E());
		  Jet.push_back(tempjet);
	  }

	for(int i=0;i<BJet.size();i++)sumBJet+=BJet.at(i);
	for(int i=0;i<Jet.size();i++)sumJet+=Jet.at(i);
/*
	  _floatVars["M_E_BJet"]     = (sumBJet).E();
	  _floatVars["M_Pt_BJet"]     = (sumBJet).Pt();
	  _floatVars["M_Eta_BJet"]     = (sumBJet).Eta();
	  _floatVars["M_Phi_BJet"]     = (sumBJet).Phi();
*/

	  _floatVars["M_Mass_Jet"]          = (sumJet).M();
	  _floatVars["M_E_Jet"]          = (sumJet).E();
	  _floatVars["M_Pt_Jet"]         = (sumJet).Pt();
	  _floatVars["M_Eta_Jet"]        = (sumJet).Eta();
	  _floatVars["M_Phi_Jet"]        = (sumJet).Phi();



	  //Get bstar
	  TLorentzVector bstar(0,0,0,0);
	 // TLorentzVector Top_p4(0,0,0,0);
	  TLorentzVector BoostedJet(0,0,0,0);
	  //Top_p4.SetPtEtaPhiM(Top.Pt(),Top.Eta(),Top.Phi(),Top.M());
	  for(int i=0;i<evtObj->boostedjetsToUsePtr->size();i++){	
		  BoostedJet.SetPtEtaPhiM(evtObj->boostedjetsToUsePtr->at(i).Pt(),evtObj->boostedjetsToUsePtr->at(i).Eta(),evtObj->boostedjetsToUsePtr->at(i).Phi(),evtObj->boostedjetsToUsePtr->at(i).M());
		//  BoostedJet.Setprunedmass(evtObj->boostedjetsToUsePtr->at(i).prunedmass());
		  bstar = Top_p4+BoostedJet;

		  _floatVars["M_Mass_BoostedJet"]          = (BoostedJet).M();
		  _floatVars["M_prunedmass_BoostedJet"]          = evtObj->boostedjetsToUsePtr->at(i).prunedmass();
		  _floatVars["M_tau21_BoostedJet"]          = evtObj->boostedjetsToUsePtr->at(i).tau2()/evtObj->boostedjetsToUsePtr->at(i).tau1();
		  _floatVars["M_tau32_BoostedJet"]          = evtObj->boostedjetsToUsePtr->at(i).tau3()/evtObj->boostedjetsToUsePtr->at(i).tau2();
		  _floatVars["M_E_BoostedJet"]          = (BoostedJet).E();
		  _floatVars["M_Pt_BoostedJet"]         = (BoostedJet).Pt();
		  _floatVars["M_Eta_BoostedJet"]        = (BoostedJet).Eta();
		  _floatVars["M_Phi_BoostedJet"]        = (BoostedJet).Phi();

		  _floatVars["M_Mass_bstar"]       = (bstar).M();
		  _floatVars["M_E_bstar"]          = (bstar).E();
		  _floatVars["M_Pt_bstar"]         = (bstar).Pt();
		  _floatVars["M_Eta_bstar"]        = (bstar).Eta();
		  _floatVars["M_Phi_bstar"]        = (bstar).Phi();
	  }


	}
	}
