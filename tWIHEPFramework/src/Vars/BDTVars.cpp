/******************************************************************************
 * testVar.cpp                                                                *
 *
 * Testing out some variables being added to the skim tree.
 *
 * History
 *      10th Nov 2016 - Created by D. Leggat. Just a test, so probably not that important
 ******************************************************************************/

#include "SingleTopRootAnalysis/Vars/BDTVars.hpp"

#include <iostream>

//Test out a couple of variables, one int and one float I guess
BDTVars::BDTVars(bool makeHistos){

  SetName("BDTVars");
  
  _intVars["M_nJet2030"] = 10;
  _intVars["M_nJet3040"] = 10;
  _intVars["M_nJet4000"] = 10;
  _intVars["M_nJet2040"] = 10;
  _intVars["M_nBJet2030"] = 10;
  _intVars["M_nBJet3040"] = 10;
  _intVars["M_nBJet4000"] = 10;
  _intVars["M_nBJet2040"] = 10;
  _intVars["M_nJet2030e24"] = 10;
  _intVars["M_nJet3040e24"] = 10;
  _intVars["M_nJet4000e24"] = 10;
  _intVars["M_nJet2040e24"] = 10;
  _intVars["M_nJet20"] = 10;
  _intVars["M_nJet25"] = 10;
  _intVars["M_nJet30"] = 10;
  _intVars["M_nJet35"] = 10;
  _intVars["M_nJet40"] = 10;
  _intVars["M_nBJet20"] = 10;
  _intVars["M_nBJet25"] = 10;
  _intVars["M_nBJet30"] = 10;
  _intVars["M_nBJet35"] = 10;
  _intVars["M_nBJet40"] = 10;
  _intVars["M_nBJet50"] = 10;
  _intVars["M_nBJet70"] = 10;
  _intVars["M_nBLJet20"] = 10;
  _intVars["M_nBLJet25"] = 10;
  _intVars["M_nBLJet30"] = 10;
  _intVars["M_nBLJet35"] = 10;
  _intVars["M_nBLJet40"] = 10;
  _intVars["M_nJet3500e24"] = 10;

  _floatVars["M_likelihood_BJet4000"] = 1.;
  _floatVars["M_likelihood_BJet2030"] = 1.;
  _floatVars["M_likelihood_BJet3040"] = 1.;
  _floatVars["M_likelihood_BJet2040"] = 1.;

  _floatVars["M_Mass_Wlv"] = 300.;
  _floatVars["M_E_Wlv"] = 500.;
  _floatVars["M_Pt_Wlv"] = 500.;
  _floatVars["M_Eta_Wlv"] = -6.;
  _floatVars["M_Phi_Wlv"] = -3.2;
  _floatVars["M_Mass_AllJets"] = 1500.;
  _floatVars["M_E_AllJets"] = 500.;
  _floatVars["M_Pt_AllJets"] = 500.;
  _floatVars["M_Eta_AllJets"] = -6.;
  _floatVars["M_Phi_AllJets"] = -3.2;

  _floatVars["M_Mass_AllJets2030"] = 2000;
  _floatVars["M_E_AllJets2030"] = 100;
  _floatVars["M_Pt_AllJets2030"] = 100;
  _floatVars["M_Eta_AllJets2030"] = -6.;
  _floatVars["M_Phi_AllJets2030"] = -3.2;
  _floatVars["M_Mass_AllJets2040"] = 3000.;
  _floatVars["M_E_AllJets2040"] = 100.;
  _floatVars["M_Pt_AllJets2040"] = 100.;
  _floatVars["M_Eta_AllJets2040"] = -6.;
  _floatVars["M_Phi_AllJets2040"] = -3.2;
  _floatVars["M_Mass_AllJets3040"] = 2000.;
  _floatVars["M_E_AllJets3040"] = 100.;
  _floatVars["M_Pt_AllJets3040"] = 100.;
  _floatVars["M_Eta_AllJets3040"] = -6.;
  _floatVars["M_Phi_AllJets3040"] = -3.2;
  _floatVars["M_Mass_AllJets4000"] = 5000.;
  _floatVars["M_E_AllJets4000"] = 300.;
  _floatVars["M_Pt_AllJets4000"] = 300.;
  _floatVars["M_Eta_AllJets4000"] = -6.;
  _floatVars["M_Phi_AllJets4000"] = -3.2;

  _floatVars["M_Mass_LeadJet"] = 100.;
  _floatVars["M_E_LeadJet"] = 300.;
  _floatVars["M_Pt_LeadJet"] = 300.;
  _floatVars["M_Eta_LeadJet"] = -6.;
  _floatVars["M_Phi_LeadJet"] = -3.2;

  _floatVars["M_Pt_Jet1_2030"] = 30.;
  _floatVars["M_Eta_Jet1_2030"] = -6.;
  _floatVars["M_Phi_Jet1_2030"] = -3.2;
  _floatVars["M_Pt_Jet1_2040"] = 40.;
  _floatVars["M_Eta_Jet1_2040"] = -6.;
  _floatVars["M_Phi_Jet1_2040"] = -3.2;
  _floatVars["M_Pt_Jet1_3040"] = 40.;
  _floatVars["M_Eta_Jet1_3040"] = -6.;
  _floatVars["M_Phi_Jet1_3040"] = -3.2;
  _floatVars["M_Pt_Jet1_4000"] = 300;
  _floatVars["M_Eta_Jet1_4000"] = -6.;
  _floatVars["M_Phi_Jet1_4000"] = -3.2;
  
  _floatVars["M_Pt_BJet1_2030"] = 30.;
  _floatVars["M_Eta_BJet1_2030"] = -6.;
  _floatVars["M_Phi_BJet1_2030"] = -3.2;
  _floatVars["M_Pt_BJet1_2040"] = 40.;
  _floatVars["M_Eta_BJet1_2040"] = -6.;
  _floatVars["M_Phi_BJet1_2040"] = -3.2;
  _floatVars["M_Pt_BJet1_3040"] = 40.;
  _floatVars["M_Eta_BJet1_3040"] = -6.;
  _floatVars["M_Phi_BJet1_3040"] = -3.2;
  _floatVars["M_Pt_BJet1_4000"] = 300.;
  _floatVars["M_Eta_BJet1_4000"] = -6.;
  _floatVars["M_Phi_BJet1_4000"] = -3.2;
  
  _floatVars["M_Pt_Jet1_2030e24"] = 30;
  _floatVars["M_Eta_Jet1_2030e24"] = -6.;
  _floatVars["M_Phi_Jet1_2030e24"] = -3.2;
  _floatVars["M_Pt_Jet1_2040e24"] = 40;
  _floatVars["M_Eta_Jet1_2040e24"] = -6.;
  _floatVars["M_Phi_Jet1_2040e24"] = -3.2;
  _floatVars["M_Pt_Jet1_3040e24"] = 40;
  _floatVars["M_Eta_Jet1_3040e24"] = -6.;
  _floatVars["M_Phi_Jet1_3040e24"] = -3.2;
  _floatVars["M_Pt_Jet1_4000e24"] = 300;
  _floatVars["M_Eta_Jet1_4000e24"] = -6.;
  _floatVars["M_Phi_Jet1_4000e24"] = -3.2;
  
  _floatVars["M_Mass_Lepton"] = 100.;
  _floatVars["M_E_Lepton"] = 200.;
  _floatVars["M_Pt_Lepton"] = 300.;
  _floatVars["M_Eta_Lepton"] = -6.;
  _floatVars["M_Phi_Lepton"] = -3.2;
  _floatVars["M_Pt_sys"] = 300;
  _floatVars["M_Pt_sys_sig"] = 70;
  _floatVars["M_Pt_sys_LeptonJet1"] = 250;
  _floatVars["M_eta_sys_LeptonJet1"] = -6.;
  
  _floatVars["M_Mass_AllJetsLepton"] = 1200;
  _floatVars["M_E_AllJetsLepton"] = 300;
  _floatVars["M_Pt_AllJetsLepton"] = 300;
  _floatVars["M_Eta_AllJetsLepton"] = -6.;
  _floatVars["M_Phi_AllJetsLepton"] = -3.2;
  _floatVars["M_Mass_AllJetsMET"] = 1800;
  _floatVars["M_E_AllJetsMET"] = 300;
  _floatVars["M_Pt_AllJetsMET"] = 300;
  _floatVars["M_Eta_AllJetsMET"] = -6.;
  _floatVars["M_Phi_AllJetsMET"] = -3.2;
  _floatVars["M_Mass_AllJetsLeptonMET"] = 2000;
  _floatVars["M_E_AllJetsLeptonMET"] = 300;
  _floatVars["M_Pt_AllJetsLeptonMET"] = 300;
  _floatVars["M_Eta_AllJetsLeptonMET"] = -6.;
  _floatVars["M_Phi_AllJetsLeptonMET"] = -3.2;
  
  _floatVars["M_AllJetsLepton_Centrality"] = 1.;
  _floatVars["M_Shat"] = 1000.;
  _floatVars["M_DeltaEtaWlvLepton"] = 6.;
  _floatVars["M_DeltaPhiWlvLepton"] = 3.2;
  _floatVars["M_DeltaRWlvLepton"] = 6.;
  
  _floatVars["M_DeltaEtaWlvMET"] = 6.;
  _floatVars["M_DeltaPhiWlvMET"] = 3.2;
  _floatVars["M_DeltaRWlvMET"] = 6.;
  
  _floatVars["M_DeltaEtaLeptonMET"] = 6.;
  _floatVars["M_DeltaPhiLeptonMET"] = 3.2;
  _floatVars["M_DeltaRLeptonMET"] = 6.;
  
  _floatVars["M_DeltaEtaJet1MET"] = 6.;
  _floatVars["M_DeltaPhiJet1MET"] = 3.2;
  _floatVars["M_DeltaRJet1MET"] = 6.;
  
  _floatVars["M_Mass_LeptonJet1"] = 1000.;
  _floatVars["M_E_LeptonJet1"] = 300;
  _floatVars["M_Pt_LeptonJet1"] = 300;
  _floatVars["M_Eta_LeptonJet1"] = -6.;
  _floatVars["M_Phi_LeptonJet1"] = -3.2;
  
  _floatVars["M_Mass_Jet1MET"] = 1200.;
  _floatVars["M_E_Jet1MET"] = 300;
  _floatVars["M_Pt_Jet1MET"] = 300;
  _floatVars["M_Eta_Jet1MET"] = -6.;
  _floatVars["M_Phi_Jet1MET"] = -3.2;
  
  _floatVars["M_Mass_Jet1LeptonMET"] = 1300.;
  _floatVars["M_E_Jet1LeptonMET"] = 300;
  _floatVars["M_Pt_Jet1LeptonMET"] = 300;
  _floatVars["M_Eta_Jet1LeptonMET"] = -6.;
  _floatVars["M_Phi_Jet1LeptonMET"] = -3.2;
  
  _floatVars["M_Pt_RHT"] = 15.;
  
  _floatVars["M_DeltaEtaLeptonJet1"] = 6.;
  _floatVars["M_DeltaPhiLeptonJet1"] = 3.4;
  _floatVars["M_DeltaRLeptonJet1"] = 6.;
  
  _floatVars["M_DeltaEtaWlvJet1"] = 6.;
  _floatVars["M_DeltaPhiWlvJet1"] = 3.4;
  _floatVars["M_DeltaRWlvJet1"] = 6.;
  
  _floatVars["M_Mass_Jet1Jet2"] = 1000.;
  _floatVars["M_E_Jet1Jet2"] = 300;
  _floatVars["M_Pt_Jet1Jet2"] = 300;
  _floatVars["M_Eta_Jet1Jet2"] = -6.;
  _floatVars["M_Phi_Jet1Jet2"] = -3.2;
  
  _floatVars["M_Mass_Jet1Jet3"] = 1000.;
  _floatVars["M_E_Jet1Jet3"] = 300;
  _floatVars["M_Pt_Jet1Jet3"] = 300;
  _floatVars["M_Eta_Jet1Jet3"] = -6.;
  _floatVars["M_Phi_Jet1Jet3"] = -3.2;
  
  _floatVars["M_Mass_Jet2Jet3"] = 1000.;
  _floatVars["M_E_Jet2Jet3"] = 300;
  _floatVars["M_Pt_Jet2Jet3"] = 300;
  _floatVars["M_Eta_Jet2Jet3"] = -6.;
  _floatVars["M_Phi_Jet2Jet3"] = -3.2;
  
  _floatVars["M_Mass_Jet1Jet2Jet3"] = 2000.;
  _floatVars["M_E_Jet1Jet2Jet3"] = 300;
  _floatVars["M_Pt_Jet1Jet2Jet3"] = 300;
  _floatVars["M_Eta_Jet1Jet2Jet3"] = -6.;
  _floatVars["M_Phi_Jet1Jet2Jet3"] = -3.2;
  
  _floatVars["M_Mass_LeptonJet2"] = 1000.;
  _floatVars["M_E_LeptonJet2"] = 300;
  _floatVars["M_Pt_LeptonJet2"] = 300;
  _floatVars["M_Eta_LeptonJet2"] = -6.;
  _floatVars["M_Phi_LeptonJet2"] = -3.2;
  
  _floatVars["M_Mass_LeptonJet3"] = 500.;
  _floatVars["M_E_LeptonJet3"] = 300;
  _floatVars["M_Pt_LeptonJet3"] = 300;
  _floatVars["M_Eta_LeptonJet3"] = -6.;
  _floatVars["M_Phi_LeptonJet3"] = -3.2;
  
  _floatVars["M_Mass_Jet2MET"] = 800.;
  _floatVars["M_E_Jet2MET"] = 300;
  _floatVars["M_Pt_Jet2MET"] = 300;
  _floatVars["M_Eta_Jet2MET"] = -6.;
  _floatVars["M_Phi_Jet2MET"] = -3.2;
  
  _floatVars["M_Mass_Jet3MET"] = 500.;
  _floatVars["M_E_Jet3MET"] = 300;
  _floatVars["M_Pt_Jet3MET"] = 300;
  _floatVars["M_Eta_Jet3MET"] = -6.;
  _floatVars["M_Phi_Jet3MET"] = -3.2;
  
  _floatVars["M_Mass_Jet1Jet2Lepton"] = 1800.;
  _floatVars["M_E_Jet1Jet2Lepton"] = 300;
  _floatVars["M_Pt_Jet1Jet2Lepton"] = 300;
  _floatVars["M_Eta_Jet1Jet2Lepton"] = -6.;
  _floatVars["M_Phi_Jet1Jet2Lepton"] = -3.2;
  
  _floatVars["M_Mass_Jet1Jet3Lepton"] = 1200.;
  _floatVars["M_E_Jet1Jet3Lepton"] = 300;
  _floatVars["M_Pt_Jet1Jet3Lepton"] = 300;
  _floatVars["M_Eta_Jet1Jet3Lepton"] = -6.;
  _floatVars["M_Phi_Jet1Jet3Lepton"] = -3.2;
  
  _floatVars["M_Mass_Jet2Jet3Lepton"] = 1200.;
  _floatVars["M_E_Jet2Jet3Lepton"] = 300;
  _floatVars["M_Pt_Jet2Jet3Lepton"] = 300;
  _floatVars["M_Eta_Jet2Jet3Lepton"] = -6.;
  _floatVars["M_Phi_Jet2Jet3Lepton"] = -3.2;
  
  _floatVars["M_Mass_Jet1Jet2Jet3Lepton"] = 2000.;
  _floatVars["M_E_Jet1Jet2Jet3Lepton"] = 300;
  _floatVars["M_Pt_Jet1Jet2Jet3Lepton"] = 300;
  _floatVars["M_Eta_Jet1Jet2Jet3Lepton"] = -6.;
  _floatVars["M_Phi_Jet1Jet2Jet3Lepton"] = -3.2;
  
  _floatVars["M_Mass_Jet1Jet2MET"] = 1800.;
  _floatVars["M_E_Jet1Jet2MET"] = 300;
  _floatVars["M_Pt_Jet1Jet2MET"] = 300;
  _floatVars["M_Eta_Jet1Jet2MET"] = -6.;
  _floatVars["M_Phi_Jet1Jet2MET"] = -3.2;
  
  _floatVars["M_Mass_Jet1Jet3MET"] = 1300.;
  _floatVars["M_E_Jet1Jet3MET"] = 300;
  _floatVars["M_Pt_Jet1Jet3MET"] = 300;
  _floatVars["M_Eta_Jet1Jet3MET"] = -6.;
  _floatVars["M_Phi_Jet1Jet3MET"] = -3.2;
  
  _floatVars["M_Mass_Jet2Jet3MET"] = 1000.;
  _floatVars["M_E_Jet2Jet3MET"] = 300;
  _floatVars["M_Pt_Jet2Jet3MET"] = 300;
  _floatVars["M_Eta_Jet2Jet3MET"] = -6.;
  _floatVars["M_Phi_Jet2Jet3MET"] = -3.2;
  
  _floatVars["M_Mass_Jet1Jet2Jet3MET"] = 1800.;
  _floatVars["M_E_Jet1Jet2Jet3MET"] = 300;
  _floatVars["M_Pt_Jet1Jet2Jet3MET"] = 300;
  _floatVars["M_Eta_Jet1Jet2Jet3MET"] = -6.;
  _floatVars["M_Phi_Jet1Jet2Jet3MET"] = -3.2;
  
  _floatVars["M_Mass_LeptonJet2MET"] = 1000.;
  _floatVars["M_E_LeptonJet2MET"] = 300;
  _floatVars["M_Pt_LeptonJet2MET"] = 300;
  _floatVars["M_Eta_LeptonJet2MET"] = -6.;
  _floatVars["M_Phi_LeptonJet2MET"] = -3.2;
  
  _floatVars["M_Mass_LeptonJet1Jet2MET"] = 1800.;
  _floatVars["M_E_LeptonJet1Jet2MET"] = 300;
  _floatVars["M_Pt_LeptonJet1Jet2MET"] = 300;
  _floatVars["M_Eta_LeptonJet1Jet2MET"] = -6.;
  _floatVars["M_Phi_LeptonJet1Jet2MET"] = -3.2;
  
  _floatVars["M_Mass_LeptonJet1Jet3MET"] = 1500.;
  _floatVars["M_E_LeptonJet1Jet3MET"] = 300;
  _floatVars["M_Pt_LeptonJet1Jet3MET"] = 300;
  _floatVars["M_Eta_LeptonJet1Jet3MET"] = -6.;
  _floatVars["M_Phi_LeptonJet1Jet3MET"] = -3.2;
  
  _floatVars["M_Mass_LeptonJet2Jet3MET"] = 1300.;
  _floatVars["M_E_LeptonJet2Jet3MET"] = 300;
  _floatVars["M_Pt_LeptonJet2Jet3MET"] = 300;
  _floatVars["M_Eta_LeptonJet2Jet3MET"] = -6.;
  _floatVars["M_Phi_LeptonJet2Jet3MET"] = -3.2;
  
  _floatVars["M_Mass_Jet1Jet2Jet3LeptonMET"] = 2000.;
  _floatVars["M_E_Jet1Jet2Jet3LeptonMET"] = 300;
  _floatVars["M_Pt_Jet1Jet2Jet3LeptonMET"] = 300;
  _floatVars["M_Eta_Jet1Jet2Jet3LeptonMET"] = -6.;
  _floatVars["M_Phi_Jet1Jet2Jet3LeptonMET"] = -3.2;
  
  _floatVars["M_TMass_Jet1Jet2"] = 800.;
  _floatVars["M_TMass_Jet1Jet3"] = 900.;
  _floatVars["M_TMass_Jet2Jet3"] = 500.;
  _floatVars["M_TMass_Jet1Jet2Jet3"] = 1000.;
  
  _floatVars["M_DeltaRJet1Jet2"] = 6.;
  _floatVars["M_DeltaEtaJet1Jet2"] = 6.;
  _floatVars["M_DeltaPhiJet1Jet2"] = 3.2;
  
  _floatVars["M_DeltaRJet1Jet3"] = 6.;
  _floatVars["M_DeltaEtaJet1Jet3"] = 6.;
  _floatVars["M_DeltaPhiJet1Jet3"] = 3.2;
  
  _floatVars["M_DeltaRJet2Jet3"] = 6.;
  _floatVars["M_DeltaEtaJet2Jet3"] = 6.;
  _floatVars["M_DeltaPhiJet2Jet3"] = 3.2;
  
  _floatVars["M_DeltaRJet2MET"] = 6.;
  _floatVars["M_DeltaEtaJet2MET"] = 6.;
  _floatVars["M_DeltaPhiJet2MET"] = 3.2;
  
  _floatVars["M_DeltaRJet3MET"] = 6.;
  _floatVars["M_DeltaEtaJet3MET"] = 6.;
  _floatVars["M_DeltaPhiJet3MET"] = 3.2;
  
  _floatVars["M_DeltaRLeptonJet2"] = 6.;
  _floatVars["M_DeltaEtaLeptonJet2"] = 6.;
  _floatVars["M_DeltaPhiLeptonJet2"] = 3.2;
  
  _floatVars["M_DeltaRWlvJet2"] = 6.;
  _floatVars["M_DeltaEtaWlvJet2"] = 6.;
  _floatVars["M_DeltaPhiWlvJet2"] = 3.2;
  
  _floatVars["M_DeltaRLeptonJet3"] = 6.;
  _floatVars["M_DeltaEtaLeptonJet3"] = 6.;
  _floatVars["M_DeltaPhiLeptonJet3"] = 3.2;
  
  _floatVars["M_DeltaRWlvJet3"] = 6.;
  _floatVars["M_DeltaEtaWlvJet3"] = 6.;
  _floatVars["M_DeltaPhiWlvJet3"] = 3.2;
  
  _floatVars["M_Jet1Jet2_Centrality"] = 1.;
  
  _floatVars["M_Jet1Jet3_Centrality"] = 1.;
  
  _floatVars["M_Jet2Jet3_Centrality"] = 1.;
  
  _floatVars["M_Jet1Jet2Jet3_Centrality"] = 1.;
  
  _floatVars["largestCSV"] = 1.;

  _floatVars["firstJetCSV"] = 1.;
  _floatVars["secondJetCSV"] = 1.;
  _floatVars["thirdJetCSV"] = 1.;
  
  _floatVars["taggedJetCSV"] = 1.;
  _floatVars["lightJet1CSV"] = 1.;
  _floatVars["lightJet2CSV"] = 1.;
  
  _floatVars["M_hadronicWmass"] = 500.;
  
  _floatVars["M_hadronicWPt"] = 500.;
  
  _floatVars["M_hadronicWEta"] = -6.;
  
  _floatVars["M_hadronicWPhi"] = -3.2;
  
  _floatVars["M_topMass2"] = 500.;
  _floatVars["M_topMassLep"] = 500.;
  _floatVars["M_topMass2_lep"] = 500.;
  _floatVars["M_topMassLep_had"] = 500.;
  
  _floatVars["M_cosThetaStar"] = -1.;
  _floatVars["M_cosThetaStar_lepOnly"] = -1.;
  _floatVars["M_cosThetaStar_hadOnly"] = -1.;

  _floatVars["M_DeltaRlightjets"] = 6.;
  
  _floatVars["M_DeltaRBJethadronicW"] = 6.;
  
  _floatVars["M_DeltaRBJetLepton"] = 6.;
  
  _floatVars["M_maxDeltaRBJetLightJet"] = 6.;
  
  _floatVars["M_HT"] = 800.;

  SetDoHists(makeHistos);

}

void BDTVars::FillBranches(EventContainer * evtObj){

  //Initialise vectors;
  Jet.clear();
  Jet2030.clear();
  Jet2040.clear();
  Jet3040.clear();
  Jet4000.clear();
  BJet.clear();
  BJet2030.clear();
  BJet2040.clear();
  BJet3040.clear();
  BJet4000.clear();
  Jete24.clear();
  Jet2030e24.clear();
  Jet2040e24.clear();
  Jet3040e24.clear();
  Jet4000e24.clear();
  Jet3500e24.clear();

  totalJets.SetPtEtaPhiE(0,0,0,0);
  sumJets2030.SetPtEtaPhiE(0,0,0,0);
  sumJets2040.SetPtEtaPhiE(0,0,0,0);
  sumJets3040.SetPtEtaPhiE(0,0,0,0);
  sumJets4000.SetPtEtaPhiE(0,0,0,0);

  //First thing we do is split up the jets into the vectors defined by pt
  TLorentzVector tempjet(0,0,0,0);
  EventTree* eventTree = evtObj->GetEventTree();
  for (auto jet : evtObj->alljets){
    tempjet.SetPtEtaPhiE(jet.Pt(),jet.Eta(),jet.Phi(),jet.E());
    //  for (int i = 0; i < eventTree-> Jet_pt->size(); i++){
    //tempjet.SetPtEtaPhiE(eventTree-> Jet_pt->at(i),eventTree->Jet_eta->at(i),eventTree->Jet_phi->at(i),eventTree->Jet_energy->at(i));
    if (tempjet.Pt() > 20 && tempjet.Pt() < 30) Jet2030.push_back(tempjet);
    if (tempjet.Pt() > 20 && tempjet.Pt() < 40) Jet2040.push_back(tempjet);
    if (tempjet.Pt() > 30 && tempjet.Pt() < 40) Jet3040.push_back(tempjet);
    if (tempjet.Pt() > 40) Jet4000.push_back(tempjet);
  }

  for (auto jet : evtObj->jets){
    tempjet.SetPtEtaPhiE(jet.Pt(),jet.Eta(),jet.Phi(),jet.E());
    Jet.push_back(tempjet);
    totalJets += tempjet;
    Jete24.push_back(tempjet);
    if (jet.Pt() > 20 && jet.Pt() < 30) Jet2030e24.push_back(tempjet);
    if (jet.Pt() > 20 && jet.Pt() < 40) Jet2040e24.push_back(tempjet);
    if (jet.Pt() > 30 && jet.Pt() < 40) Jet3040e24.push_back(tempjet);
    if (jet.Pt() > 40) Jet4000e24.push_back(tempjet);
  }
  
  for (auto jet : evtObj->taggedJets){
    tempjet.SetPtEtaPhiE(jet.Pt(),jet.Eta(),jet.Phi(),jet.E());
    BJet.push_back(tempjet);
    if (jet.Pt() > 20 && jet.Pt() < 30) BJet2030.push_back(tempjet);
    if (jet.Pt() > 20 && jet.Pt() < 40) BJet2040.push_back(tempjet);
    if (jet.Pt() > 30 && jet.Pt() < 40) BJet3040.push_back(tempjet);
    if (jet.Pt() > 40) BJet4000.push_back(tempjet);
  }

  //Now set up the lepton and met variables
  TLorentzVector Lepton(00,0,0,0);
  TLorentzVector Miss(00,0,0,0);
  TLorentzVector Wlv(0,0,0,0);
  
  Miss.SetPtEtaPhiE(evtObj->missingEt_xy,0,evtObj->missingPhi_xy,evtObj->missingEt_xy);
  if (evtObj->electronsToUsePtr->size() > 0){ // if this number is >0 we're in the electron channel. Otherwise use muons
    Lepton.SetPtEtaPhiE(evtObj->electronsToUsePtr->at(0).Pt(),evtObj->electronsToUsePtr->at(0).Eta(),evtObj->electronsToUsePtr->at(0).Phi(),evtObj->electronsToUsePtr->at(0).E());
  }
  else{
    Lepton.SetPtEtaPhiE(evtObj->muonsToUsePtr->at(0).Pt(),evtObj->muonsToUsePtr->at(0).Eta(),evtObj->muonsToUsePtr->at(0).Phi(),evtObj->muonsToUsePtr->at(0).E());
  }
  Wlv = Lepton+Miss;
  
  _floatVars["M_Mass_Wlv"]       = (Wlv).M();
  _floatVars["M_E_Wlv"]          = (Wlv).E();
  _floatVars["M_Pt_Wlv"]         = (Wlv).Pt();
  _floatVars["M_Eta_Wlv"]        = (Wlv).Eta();
  _floatVars["M_Phi_Wlv"]        = (Wlv).Phi();

  _floatVars["M_DeltaEtaWlvLepton"] =  fabs(Wlv.Eta()-Lepton.Eta());
  _floatVars["M_DeltaPhiWlvLepton"] =  fabs(Wlv.DeltaPhi(Lepton));
  _floatVars["M_DeltaRWlvLepton"] =  fabs(Wlv.DeltaR(Lepton));

  _floatVars["M_DeltaEtaWlvMET"] =  fabs(Miss.Eta()-Wlv.Eta());
  _floatVars["M_DeltaPhiWlvMET"] =  fabs(Miss.DeltaPhi(Wlv));
  _floatVars["M_DeltaRWlvMET"] =  fabs(Miss.DeltaR(Wlv));

  _floatVars["M_Mass_Lepton"]       = (Lepton).M();
  _floatVars["M_E_Lepton"]          = (Lepton).E();
  _floatVars["M_Pt_Lepton"]         = (Lepton).Pt();
  _floatVars["M_Eta_Lepton"]        = (Lepton).Eta();
  _floatVars["M_Phi_Lepton"]        = (Lepton).Phi();
  

  _floatVars["M_DeltaEtaLeptonMET"] =  fabs(Lepton.Eta()-Miss.Eta());
  _floatVars["M_DeltaPhiLeptonMET"] =  fabs(Miss.DeltaPhi(Lepton));
  _floatVars["M_DeltaRLeptonMET"] =  fabs(Miss.DeltaR(Lepton));


  _intVars["M_nJet2030"]   = Jet2030.size();
  for(int i=0;i<_intVars["M_nJet2030"];i++)
    {sumJets2030 += Jet2030.at(i);}


  _intVars["M_nJet2040"]   = Jet2040.size();
  for(int i=0;i<_intVars["M_nJet2040"];i++)
    {sumJets2040 += Jet2040.at(i);}


  _intVars["M_nJet3040"]   = Jet3040.size();
  for(int i=0;i<_intVars["M_nJet3040"];i++)
    {sumJets3040 += Jet3040.at(i);}

  _intVars["M_nJet4000"]   = Jet4000.size();
  for(int i=0;i<_intVars["M_nJet4000"];i++)
    {sumJets4000 += Jet4000.at(i);}

  _intVars["M_nBJet2030"]   = BJet2030.size();
  _intVars["M_nBJet2040"]   = BJet2040.size();
  _intVars["M_nBJet3040"]   = BJet3040.size();
  _intVars["M_nBJet4000"]   = BJet4000.size();

  _intVars["M_nJet2030e24"]   = Jet2030e24.size();
  _intVars["M_nJet2040e24"]   = Jet2040e24.size();
  _intVars["M_nJet3040e24"]   = Jet3040e24.size();
  _intVars["M_nJet4000e24"]   = Jet4000e24.size();

  if(_intVars["M_nJet2040"]>=1){
    _floatVars["M_Pt_Jet1_2040"]     = (Jet2040.at(0)).Pt();
    _floatVars["M_Eta_Jet1_2040"]     = (Jet2040.at(0)).Eta();
    _floatVars["M_Phi_Jet1_2040"]     = (Jet2040.at(0)).Phi();
    _floatVars["M_Mass_AllJets2040"]   = (sumJets2040).M();
    _floatVars["M_E_AllJets2040"]   = (sumJets2040).E();
    _floatVars["M_Pt_AllJets2040"]   = (sumJets2040).Pt();
    _floatVars["M_Eta_AllJets2040"]   = (sumJets2040).Eta();
    _floatVars["M_Phi_AllJets2040"]   = (sumJets2040).Phi();
  }

  if(_intVars["M_nJet2030"]>=1){
    _floatVars["M_Pt_Jet1_2030"]     = (Jet2030.at(0)).Pt();
    _floatVars["M_Eta_Jet1_2030"]     = (Jet2030.at(0)).Eta();
    _floatVars["M_Phi_Jet1_2030"]     = (Jet2030.at(0)).Phi();
    _floatVars["M_Mass_AllJets2030"]   = (sumJets2030).M();
    _floatVars["M_E_AllJets2030"]   = (sumJets2030).E();
    _floatVars["M_Pt_AllJets2030"]   = (sumJets2030).Pt();
    _floatVars["M_Eta_AllJets2030"]   = (sumJets2030).Eta();
    _floatVars["M_Phi_AllJets2030"]   = (sumJets2030).Phi();
  }



  if(_intVars["M_nJet3040"]>=1){
    _floatVars["M_Pt_Jet1_3040"]     = (Jet3040.at(0)).Pt();
    _floatVars["M_Eta_Jet1_3040"]     = (Jet3040.at(0)).Eta();
    _floatVars["M_Phi_Jet1_3040"]     = (Jet3040.at(0)).Phi();
    _floatVars["M_Mass_AllJets3040"]   = (sumJets3040).M();
    _floatVars["M_E_AllJets3040"]   = (sumJets3040).E();
    _floatVars["M_Pt_AllJets3040"]   = (sumJets3040).Pt();
    _floatVars["M_Eta_AllJets3040"]   = (sumJets3040).Eta();
    _floatVars["M_Phi_AllJets3040"]   = (sumJets3040).Phi();
  }


  if(_intVars["M_nJet4000"]>=1){
    _floatVars["M_Pt_Jet1_4000"]     = (Jet4000.at(0)).Pt();
    _floatVars["M_Eta_Jet1_4000"]     = (Jet4000.at(0)).Eta();
    _floatVars["M_Phi_Jet1_4000"]     = (Jet4000.at(0)).Phi();
    _floatVars["M_Mass_AllJets4000"]   = (sumJets4000).M();
    _floatVars["M_E_AllJets4000"]   = (sumJets4000).E();
    _floatVars["M_Pt_AllJets4000"]   = (sumJets4000).Pt();
    _floatVars["M_Eta_AllJets4000"]   = (sumJets4000).Eta();
    _floatVars["M_Phi_AllJets4000"]   = (sumJets4000).Phi();
  }


  if(_intVars["M_nBJet2040"]>=1){
    _floatVars["M_Pt_BJet1_2040"]     = (BJet2040.at(0)).Pt();
    _floatVars["M_Eta_BJet1_2040"]     = (BJet2040.at(0)).Eta();
    _floatVars["M_Phi_BJet1_2040"]     = (BJet2040.at(0)).Phi();
  }
  if(_intVars["M_nBJet2030"]>=1){
    _floatVars["M_Pt_BJet1_2030"]     = (BJet2030.at(0)).Pt();
    _floatVars["M_Eta_BJet1_2030"]     = (BJet2030.at(0)).Eta();
    _floatVars["M_Phi_BJet1_2030"]     = (BJet2030.at(0)).Phi();
  }
  if(_intVars["M_nBJet3040"]>=1){
    _floatVars["M_Pt_BJet1_3040"]     = (BJet3040.at(0)).Pt();
    _floatVars["M_Eta_BJet1_3040"]     = (BJet3040.at(0)).Eta();
    _floatVars["M_Phi_BJet1_3040"]     = (BJet3040.at(0)).Phi();
  }
  if(_intVars["M_nBJet4000"]>=1){
    _floatVars["M_Pt_BJet1_4000"]     = (BJet4000.at(0)).Pt();
    _floatVars["M_Eta_BJet1_4000"]     = (BJet4000.at(0)).Eta();
    _floatVars["M_Phi_BJet1_4000"]     = (BJet4000.at(0)).Phi();
  }

  if(_intVars["M_nJet2040e24"]>=1){
    _floatVars["M_Pt_Jet1_2040e24"]     = (Jet2040e24.at(0)).Pt();
    _floatVars["M_Eta_Jet1_2040e24"]     = (Jet2040e24.at(0)).Eta();
    _floatVars["M_Phi_Jet1_2040e24"]     = (Jet2040e24.at(0)).Phi();
  }

  if(_intVars["M_nJet2030e24"]>=1){
    _floatVars["M_Pt_Jet1_2030e24"]     = (Jet2030e24.at(0)).Pt();
    _floatVars["M_Eta_Jet1_2030e24"]     = (Jet2030e24.at(0)).Eta();
    _floatVars["M_Phi_Jet1_2030e24"]     = (Jet2030e24.at(0)).Phi();
  }

  if(_intVars["M_nJet3040e24"]>=1){
    _floatVars["M_Pt_Jet1_3040e24"]     = (Jet3040e24.at(0)).Pt();
    _floatVars["M_Eta_Jet1_3040e24"]     = (Jet3040e24.at(0)).Eta();
    _floatVars["M_Phi_Jet1_3040e24"]     = (Jet3040e24.at(0)).Phi();
  }

  if(_intVars["M_nJet4000e24"]>=1){
    _floatVars["M_Pt_Jet1_4000e24"]     = (Jet4000e24.at(0)).Pt();
    _floatVars["M_Eta_Jet1_4000e24"]     = (Jet4000e24.at(0)).Eta();
    _floatVars["M_Phi_Jet1_4000e24"]     = (Jet4000e24.at(0)).Phi();
  }

  if(Jet.size()>=1){
    _floatVars["M_Mass_AllJets"]   = (totalJets).M();
    _floatVars["M_E_AllJets"]   = (totalJets).E();
    _floatVars["M_Pt_AllJets"]   = (totalJets).Pt();
    _floatVars["M_Eta_AllJets"]   = (totalJets).Eta();
    _floatVars["M_Phi_AllJets"]   = (totalJets).Phi();
    
    _floatVars["M_Mass_AllJetsLepton"]  = (totalJets+Lepton).M();
    _floatVars["M_E_AllJetsLepton"]  = (totalJets+Lepton).E();
    _floatVars["M_Pt_AllJetsLepton"]  = (totalJets+Lepton).Pt();
    _floatVars["M_Eta_AllJetsLepton"]  = (totalJets+Lepton).Eta();
    _floatVars["M_Phi_AllJetsLepton"]  = (totalJets+Lepton).Phi();
    
    _floatVars["M_Mass_AllJetsMET"]     = (totalJets+Miss).M();
    _floatVars["M_E_AllJetsMET"]     = (totalJets+Miss).E();
    _floatVars["M_Pt_AllJetsMET"]     = (totalJets+Miss).Pt();
    _floatVars["M_Eta_AllJetsMET"]     = (totalJets+Miss).Eta();
    _floatVars["M_Phi_AllJetsMET"]     = (totalJets+Miss).Phi();
    
    _floatVars["M_Mass_AllJetsLeptonMET"]     = (totalJets + Lepton + Miss).M();
    _floatVars["M_E_AllJetsLeptonMET"]     = (totalJets + Lepton + Miss).E();
    _floatVars["M_Pt_AllJetsLeptonMET"]     = (totalJets + Lepton + Miss).Pt();
    _floatVars["M_Eta_AllJetsLeptonMET"]     = (totalJets + Lepton + Miss).Eta();
    _floatVars["M_Phi_AllJetsLeptonMET"]     = (totalJets + Lepton + Miss).Phi();
   
    _floatVars["M_Mass_LeadJet"]     = (Jet.at(0)).M();
    _floatVars["M_E_LeadJet"]     = (Jet.at(0)).E();
    _floatVars["M_Pt_LeadJet"]     = (Jet.at(0)).Pt();
    _floatVars["M_Eta_LeadJet"]     = (Jet.at(0)).Eta();
    _floatVars["M_Phi_LeadJet"]     = (Jet.at(0)).Phi();

    _floatVars["M_AllJetsLepton_Centrality"] = (totalJets.Pt() + Lepton.Pt())/ (totalJets.P() + Lepton.P());
    _floatVars["M_Shat"] = (totalJets+Lepton).M();
    double px_sys        = (Lepton + totalJets ).Px() + Miss.Px();
    double py_sys        = (Lepton + totalJets ).Py() + Miss.Px();
    _floatVars["M_Pt_sys"]              =  sqrt(pow(px_sys,2)+pow(py_sys,2));
    // the pT_sys_err is related to calobase measurement and track measurment, now only consider the Calo item since they are dominate
    double pT_sys_err = sqrt((Lepton + totalJets + Miss).Pt()+(Lepton + totalJets + Miss).Et());
    //System_Pt_significance => get missing parts in the transverse plate
    _floatVars["M_Pt_sys_sig"]            = _floatVars["M_Pt_sys"]/pT_sys_err;
    _floatVars["M_Pt_sys_LeptonJet1"]     = (Lepton + Jet.at(0)).Pt();
    _floatVars["M_eta_sys_LeptonJet1"]    = (Lepton + Jet.at(0)).Eta();
    
    



    _floatVars["M_Mass_LeptonJet1"] = (Jet.at(0)+Lepton).M();
    _floatVars["M_E_LeptonJet1"] = (Jet.at(0)+Lepton).E();
    _floatVars["M_Pt_LeptonJet1"] = (Jet.at(0)+Lepton).Pt();
    _floatVars["M_Eta_LeptonJet1"] = (Jet.at(0)+Lepton).Eta();
    _floatVars["M_Phi_LeptonJet1"] = (Jet.at(0)+Lepton).Phi();
    
    _floatVars["M_Mass_Jet1MET"]     = (Jet.at(0)+Miss).M();
    _floatVars["M_E_Jet1MET"]     = (Jet.at(0)+Miss).E();
    _floatVars["M_Pt_Jet1MET"]     = (Jet.at(0)+Miss).Pt();
    _floatVars["M_Eta_Jet1MET"]     = (Jet.at(0)+Miss).Eta();
    _floatVars["M_Phi_Jet1MET"]     = (Jet.at(0)+Miss).Phi();
    


    _floatVars["M_Mass_Jet1LeptonMET"]     = (Jet.at(0) + Lepton +  Miss).M();
    _floatVars["M_E_Jet1LeptonMET"]     = (Jet.at(0) + Lepton +  Miss).E();
    _floatVars["M_Pt_Jet1LeptonMET"]     = (Jet.at(0) + Lepton +  Miss).Pt();
    _floatVars["M_Eta_Jet1LeptonMET"]     = (Jet.at(0) + Lepton +  Miss).Eta();
    _floatVars["M_Phi_Jet1LeptonMET"]     = (Jet.at(0) + Lepton +  Miss).Phi();
    

    _floatVars["M_Pt_RHT"] = _floatVars["M_Pt_Jet1LeptonMET"]/_floatVars["M_Pt_AllJetsLeptonMET"];
    _floatVars["M_DeltaEtaJet1MET"] =   fabs(Jet.at(0).Eta()-Miss.Eta());
    _floatVars["M_DeltaPhiJet1MET"] =   fabs(Jet.at(0).DeltaPhi(Miss));
    _floatVars["M_DeltaRJet1MET"] =   fabs(Jet.at(0).DeltaR(Miss));
    

    _floatVars["M_DeltaEtaLeptonJet1"] =  fabs(Lepton.Eta()-Jet.at(0).Eta());
    _floatVars["M_DeltaPhiLeptonJet1"] =  fabs(Jet.at(0).DeltaPhi(Lepton));
    _floatVars["M_DeltaRLeptonJet1"] =  fabs(Jet.at(0).DeltaR(Lepton));
    _floatVars["M_DeltaEtaWlvJet1"] =  fabs(Wlv.Eta()-Jet.at(0).Eta());
    _floatVars["M_DeltaPhiWlvJet1"] =  fabs(Jet.at(0).DeltaPhi(Wlv));
    _floatVars["M_DeltaRWlvJet1"] =  fabs(Jet.at(0).DeltaR(Wlv));
    
     if(Jet.size()>=2)
      {

	_floatVars["M_Mass_Jet1Jet2"]       = (Jet.at(0)+Jet.at(1)).M();
	_floatVars["M_E_Jet1Jet2"]       = (Jet.at(0)+Jet.at(1)).E();
	_floatVars["M_Pt_Jet1Jet2"]       = (Jet.at(0)+Jet.at(1)).Pt();
	_floatVars["M_Eta_Jet1Jet2"]       = (Jet.at(0)+Jet.at(1)).Eta();
	_floatVars["M_Phi_Jet1Jet2"]       = (Jet.at(0)+Jet.at(1)).Phi();
	
	_floatVars["M_Mass_LeptonJet2"] = (Jet.at(1)+Lepton).M();
	_floatVars["M_E_LeptonJet2"] = (Jet.at(1)+Lepton).E();
	_floatVars["M_Pt_LeptonJet2"] = (Jet.at(1)+Lepton).Pt();
	_floatVars["M_Eta_LeptonJet2"] = (Jet.at(1)+Lepton).Eta();
	_floatVars["M_Phi_LeptonJet2"] = (Jet.at(1)+Lepton).Phi();
	

	_floatVars["M_Mass_Jet2MET"] = (Jet.at(1) +  Miss).M();
	_floatVars["M_E_Jet2MET"] = (Jet.at(1) +  Miss).E();
	_floatVars["M_Pt_Jet2MET"] = (Jet.at(1) +  Miss).Pt();
	_floatVars["M_Eta_Jet2MET"] = (Jet.at(1) +  Miss).Eta();
	_floatVars["M_Phi_Jet2MET"] = (Jet.at(1) +  Miss).Phi();
	
	_floatVars["M_Mass_Jet1Jet2Lepton"] = (Jet.at(0)+Jet.at(1)+Lepton).M();
	_floatVars["M_E_Jet1Jet2Lepton"] = (Jet.at(0)+Jet.at(1)+Lepton).E();
	_floatVars["M_Pt_Jet1Jet2Lepton"] = (Jet.at(0)+Jet.at(1)+Lepton).Pt();
	_floatVars["M_Eta_Jet1Jet2Lepton"] = (Jet.at(0)+Jet.at(1)+Lepton).Eta();
	_floatVars["M_Phi_Jet1Jet2Lepton"] = (Jet.at(0)+Jet.at(1)+Lepton).Phi();
	
	_floatVars["M_Mass_Jet1Jet2MET"] = (Jet.at(0) + Jet.at(1) +  Miss).M();
	_floatVars["M_E_Jet1Jet2MET"] = (Jet.at(0) + Jet.at(1) +  Miss).E();
	_floatVars["M_Pt_Jet1Jet2MET"] = (Jet.at(0) + Jet.at(1) +  Miss).Pt();
	_floatVars["M_Eta_Jet1Jet2MET"] = (Jet.at(0) + Jet.at(1) +  Miss).Eta();
	_floatVars["M_Phi_Jet1Jet2MET"] = (Jet.at(0) + Jet.at(1) +  Miss).Phi();
	
	_floatVars["M_Mass_LeptonJet2MET"] = (Lepton + Jet.at(1) + Miss).M();
	_floatVars["M_E_LeptonJet2MET"] = (Lepton + Jet.at(1) + Miss).E();
	_floatVars["M_Pt_LeptonJet2MET"] = (Lepton + Jet.at(1) + Miss).Pt();
	_floatVars["M_Eta_LeptonJet2MET"] = (Lepton + Jet.at(1) + Miss).Eta();
	_floatVars["M_Phi_LeptonJet2MET"] = (Lepton + Jet.at(1) + Miss).Phi();
	

	_floatVars["M_Mass_LeptonJet1Jet2MET"] = (Lepton + Jet.at(0)  + Jet.at(1) + Miss).M();
	_floatVars["M_E_LeptonJet1Jet2MET"] = (Lepton + Jet.at(0)  + Jet.at(1) + Miss).E();
	_floatVars["M_Pt_LeptonJet1Jet2MET"] = (Lepton + Jet.at(0)  + Jet.at(1) + Miss).Pt();
	_floatVars["M_Eta_LeptonJet1Jet2MET"] = (Lepton + Jet.at(0)  + Jet.at(1) + Miss).Eta();
	_floatVars["M_Phi_LeptonJet1Jet2MET"] = (Lepton + Jet.at(0)  + Jet.at(1) + Miss).Phi();


	Float_t tempd12=-999;
	tempd12= pow(Jet.at(0).Et() + Jet.at(1).Et(),2) - pow(Jet.at(0).Px()+Jet.at(1).Px(),2) - pow(Jet.at(0).Py()+Jet.at(1).Py(),2);
	if(tempd12>0) tempd12=TMath::Sqrt(tempd12);
	_floatVars["M_TMass_Jet1Jet2"] = tempd12;
	

	
	_floatVars["M_DeltaRJet1Jet2"]  = fabs(Jet.at(0).DeltaR(Jet.at(1)));
	_floatVars["M_DeltaEtaJet1Jet2"]= fabs(Jet.at(0).Eta()-Jet.at(1).Eta());
	_floatVars["M_DeltaPhiJet1Jet2"]= fabs(Jet.at(0).DeltaPhi(Jet.at(1)));
	_floatVars["M_DeltaPhiJet2MET"] =   fabs(Jet.at(1).DeltaPhi(Miss));
	_floatVars["M_DeltaRJet2MET"] =   fabs(Jet.at(1).DeltaR(Miss));
	_floatVars["M_DeltaEtaJet2MET"] =   fabs(Jet.at(1).Eta()-Miss.Eta());
	

	_floatVars["M_DeltaEtaLeptonJet2"] =  fabs(Lepton.Eta()-Jet.at(1).Eta());
	_floatVars["M_DeltaPhiLeptonJet2"] =  fabs(Jet.at(1).DeltaPhi(Lepton));
	_floatVars["M_DeltaRLeptonJet2"] =  fabs(Jet.at(1).DeltaR(Lepton));
	_floatVars["M_DeltaEtaWlvJet2"] =  fabs(Wlv.Eta()-Jet.at(1).Eta());
	_floatVars["M_DeltaPhiWlvJet2"] =  fabs(Jet.at(1).DeltaPhi(Wlv));
	_floatVars["M_DeltaRWlvJet2"] =  fabs(Jet.at(1).DeltaR(Wlv));
	_floatVars["M_Jet1Jet2_Centrality"] = (Jet.at(0).Pt() + Jet.at(1).Pt())/ (Jet.at(0).P() + Jet.at(1).P());
	

	if(Jet.size()>=3)
	  {
	    _floatVars["M_Mass_Jet1Jet3"]   = (Jet.at(0)+Jet.at(2)).M();
	    _floatVars["M_E_Jet1Jet3"]   = (Jet.at(0)+Jet.at(2)).E();
	    _floatVars["M_Pt_Jet1Jet3"]   = (Jet.at(0)+Jet.at(2)).Pt();
	    _floatVars["M_Eta_Jet1Jet3"]   = (Jet.at(0)+Jet.at(2)).Eta();
	    _floatVars["M_Phi_Jet1Jet3"]   = (Jet.at(0)+Jet.at(2)).Phi();
	    
	    _floatVars["M_Mass_Jet2Jet3"]      = (Jet.at(1)+Jet.at(2)).M();
	    _floatVars["M_E_Jet2Jet3"]      = (Jet.at(1)+Jet.at(2)).E();
	    _floatVars["M_Pt_Jet2Jet3"]      = (Jet.at(1)+Jet.at(2)).Pt();
	    _floatVars["M_Eta_Jet2Jet3"]      = (Jet.at(1)+Jet.at(2)).Eta();
	    _floatVars["M_Phi_Jet2Jet3"]      = (Jet.at(1)+Jet.at(2)).Phi();
	    
	    _floatVars["M_Mass_LeptonJet3"] = (Jet.at(2)+Lepton).M();
	    _floatVars["M_E_LeptonJet3"] = (Jet.at(2)+Lepton).E();
	    _floatVars["M_Pt_LeptonJet3"] = (Jet.at(2)+Lepton).Pt();
	    _floatVars["M_Eta_LeptonJet3"] = (Jet.at(2)+Lepton).Eta();
	    _floatVars["M_Phi_LeptonJet3"] = (Jet.at(2)+Lepton).Phi();
	    

	    _floatVars["M_Mass_Jet3MET"] = (Jet.at(2)+Miss).M();
	    _floatVars["M_E_Jet3MET"] = (Jet.at(2)+Miss).E();
	    _floatVars["M_Pt_Jet3MET"] = (Jet.at(2)+Miss).Pt();
	    _floatVars["M_Eta_Jet3MET"] = (Jet.at(2)+Miss).Eta();
	    _floatVars["M_Phi_Jet3MET"] = (Jet.at(2)+Miss).Phi();
	    

	    _floatVars["M_Mass_Jet1Jet2Jet3"]       = (Jet.at(0)+Jet.at(1)+Jet.at(2)).M();
	    _floatVars["M_E_Jet1Jet2Jet3"]       = (Jet.at(0)+Jet.at(1)+Jet.at(2)).E();
	    _floatVars["M_Pt_Jet1Jet2Jet3"]       = (Jet.at(0)+Jet.at(1)+Jet.at(2)).Pt();
	    _floatVars["M_Eta_Jet1Jet2Jet3"]       = (Jet.at(0)+Jet.at(1)+Jet.at(2)).Eta();
	    _floatVars["M_Phi_Jet1Jet2Jet3"]       = (Jet.at(0)+Jet.at(1)+Jet.at(2)).Phi();
	    _floatVars["M_Mass_Jet1Jet3Lepton"] = (Jet.at(0)+Jet.at(2)+Lepton).M();
	    _floatVars["M_E_Jet1Jet3Lepton"] = (Jet.at(0)+Jet.at(2)+Lepton).E();
	    _floatVars["M_Pt_Jet1Jet3Lepton"] = (Jet.at(0)+Jet.at(2)+Lepton).Pt();
	    _floatVars["M_Eta_Jet1Jet3Lepton"] = (Jet.at(0)+Jet.at(2)+Lepton).Eta();
	    _floatVars["M_Phi_Jet1Jet3Lepton"] = (Jet.at(0)+Jet.at(2)+Lepton).Phi();
	    
	    _floatVars["M_Mass_Jet2Jet3Lepton"] = (Jet.at(1)+Jet.at(2)+Lepton).M();
	    _floatVars["M_E_Jet2Jet3Lepton"] = (Jet.at(1)+Jet.at(2)+Lepton).E();
	    _floatVars["M_Pt_Jet2Jet3Lepton"] = (Jet.at(1)+Jet.at(2)+Lepton).Pt();
	    _floatVars["M_Eta_Jet2Jet3Lepton"] = (Jet.at(1)+Jet.at(2)+Lepton).Eta();
	    _floatVars["M_Phi_Jet2Jet3Lepton"] = (Jet.at(1)+Jet.at(2)+Lepton).Phi();
	    
	    _floatVars["M_Mass_Jet1Jet3MET"] = (Jet.at(0)+Jet.at(2) + Miss).M();
	    _floatVars["M_E_Jet1Jet3MET"] = (Jet.at(0)+Jet.at(2) + Miss).E();
	    _floatVars["M_Pt_Jet1Jet3MET"] = (Jet.at(0)+Jet.at(2) + Miss).Pt();
	    _floatVars["M_Eta_Jet1Jet3MET"] = (Jet.at(0)+Jet.at(2) + Miss).Eta();
	    _floatVars["M_Phi_Jet1Jet3MET"] = (Jet.at(0)+Jet.at(2) + Miss).Phi();
	    _floatVars["M_Mass_Jet2Jet3MET"] = (Jet.at(1)+Jet.at(2) + Miss).M();
	    _floatVars["M_E_Jet2Jet3MET"] = (Jet.at(1)+Jet.at(2) + Miss).E();
	    _floatVars["M_Pt_Jet2Jet3MET"] = (Jet.at(1)+Jet.at(2) + Miss).Pt();
	    _floatVars["M_Eta_Jet2Jet3MET"] = (Jet.at(1)+Jet.at(2) + Miss).Eta();
	    _floatVars["M_Phi_Jet2Jet3MET"] = (Jet.at(1)+Jet.at(2) + Miss).Phi();
	    _floatVars["M_Mass_LeptonJet2Jet3MET"] = (Lepton + Jet.at(1) +Jet.at(2) + Miss).M();
	    _floatVars["M_E_LeptonJet2Jet3MET"] = (Lepton + Jet.at(1) +Jet.at(2) + Miss).E();
	    _floatVars["M_Pt_LeptonJet2Jet3MET"] = (Lepton + Jet.at(1) +Jet.at(2) + Miss).Pt();
	    _floatVars["M_Eta_LeptonJet2Jet3MET"] = (Lepton + Jet.at(1) +Jet.at(2) + Miss).Eta();
	    _floatVars["M_Phi_LeptonJet2Jet3MET"] = (Lepton + Jet.at(1) +Jet.at(2) + Miss).Phi();
	    

	    _floatVars["M_Mass_LeptonJet1Jet3MET"] = (Lepton + Jet.at(0) + Jet.at(2) + Miss).M();
	    _floatVars["M_E_LeptonJet1Jet3MET"] = (Lepton + Jet.at(0) + Jet.at(2) + Miss).E();
	    _floatVars["M_Pt_LeptonJet1Jet3MET"] = (Lepton + Jet.at(0) + Jet.at(2) + Miss).Pt();
	    _floatVars["M_Eta_LeptonJet1Jet3MET"] = (Lepton + Jet.at(0) + Jet.at(2) + Miss).Eta();
	    _floatVars["M_Phi_LeptonJet1Jet3MET"] = (Lepton + Jet.at(0) + Jet.at(2) + Miss).Phi();
	    _floatVars["M_Mass_Jet1Jet2Jet3MET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Miss).M();
	    _floatVars["M_E_Jet1Jet2Jet3MET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Miss).E();
	    _floatVars["M_Pt_Jet1Jet2Jet3MET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Miss).Pt();
	    _floatVars["M_Eta_Jet1Jet2Jet3MET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Miss).Eta();
	    _floatVars["M_Phi_Jet1Jet2Jet3MET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Miss).Phi();
	    

	    _floatVars["M_Mass_Jet1Jet2Jet3Lepton"] = (Jet.at(0)+Jet.at(1)+Jet.at(2)+Lepton).M();
	    _floatVars["M_E_Jet1Jet2Jet3Lepton"] = (Jet.at(0)+Jet.at(1)+Jet.at(2)+Lepton).E();
	    _floatVars["M_Pt_Jet1Jet2Jet3Lepton"] = (Jet.at(0)+Jet.at(1)+Jet.at(2)+Lepton).Pt();
	    _floatVars["M_Eta_Jet1Jet2Jet3Lepton"] = (Jet.at(0)+Jet.at(1)+Jet.at(2)+Lepton).Eta();
	    _floatVars["M_Phi_Jet1Jet2Jet3Lepton"] = (Jet.at(0)+Jet.at(1)+Jet.at(2)+Lepton).Phi();
	    
	    _floatVars["M_Mass_Jet1Jet2Jet3LeptonMET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Lepton + Miss).M();
	    _floatVars["M_E_Jet1Jet2Jet3LeptonMET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Lepton + Miss).E();
	    _floatVars["M_Pt_Jet1Jet2Jet3LeptonMET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Lepton + Miss).Pt();
	    _floatVars["M_Eta_Jet1Jet2Jet3LeptonMET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Lepton + Miss).Eta();
	    _floatVars["M_Phi_Jet1Jet2Jet3LeptonMET"] = (Jet.at(0)+Jet.at(1)+Jet.at(2) + Lepton + Miss).Phi();
	    
	    Float_t tempd23=-999;
	    tempd23= pow(Jet.at(1).Et() + Jet.at(2).Et(),2) - pow(Jet.at(1).Px()+Jet.at(2).Px(),2) - pow(Jet.at(1).Py()+Jet.at(2).Py(),2);
	    if(tempd23>0) tempd23=TMath::Sqrt(tempd23);
	    _floatVars["M_TMass_Jet2Jet3"] = tempd23;
	    Float_t tempd13=-999;
	    tempd13= pow(Jet.at(0).Et() + Jet.at(2).Et(),2) - pow(Jet.at(0).Px()+Jet.at(2).Px(),2) - pow(Jet.at(0).Py()+Jet.at(2).Py(),2);
	    if(tempd13>0) tempd13=TMath::Sqrt(tempd13);
	    _floatVars["M_TMass_Jet1Jet3"] = tempd13;
	    

	    Float_t tempd123=-999;
	    tempd123= pow(Jet.at(0).Et() + Jet.at(1).Et()+ Jet.at(2).Et(),2) - pow(Jet.at(0).Px()+Jet.at(1).Px()+Jet.at(2).Px(),2) - pow(Jet.at(0).Py()+Jet.at(1).Py()+Jet.at(2).Py(),2);
	    if(tempd123>0) tempd123=TMath::Sqrt(tempd123);
	    _floatVars["M_TMass_Jet1Jet2Jet3"] = tempd123;
	    

	    _floatVars["M_DeltaRJet1Jet3"]  = fabs(Jet.at(0).DeltaR(Jet.at(2)));
	    _floatVars["M_DeltaEtaJet1Jet3"]= fabs(Jet.at(0).Eta()-Jet.at(2).Eta());
	    _floatVars["M_DeltaPhiJet1Jet3"]= fabs(Jet.at(0).DeltaPhi(Jet.at(2)));
	    
	    _floatVars["M_DeltaRJet2Jet3"]  = fabs(Jet.at(1).DeltaR(Jet.at(2)));
	    _floatVars["M_DeltaEtaJet2Jet3"]= fabs(Jet.at(1).Eta()-Jet.at(2).Eta());
	    _floatVars["M_DeltaPhiJet2Jet3"]= fabs(Jet.at(1).DeltaPhi(Jet.at(2)));
	    
	    _floatVars["M_DeltaEtaJet3MET"] =   fabs(Jet.at(2).Eta()-Miss.Eta());
	    _floatVars["M_DeltaPhiJet3MET"] =   fabs(Jet.at(2).DeltaPhi(Miss));
	    _floatVars["M_DeltaRJet3MET"] =   fabs(Jet.at(2).DeltaR(Miss));
	    

	    _floatVars["M_DeltaEtaLeptonJet3"] =  fabs(Lepton.Eta()-Jet.at(2).Eta());
	    _floatVars["M_DeltaPhiLeptonJet3"] =  fabs(Jet.at(2).DeltaPhi(Lepton));
	    _floatVars["M_DeltaRLeptonJet3"] =  fabs(Jet.at(2).DeltaR(Lepton));
	    
	    _floatVars["M_DeltaEtaWlvJet3"] =  fabs(Wlv.Eta()-Jet.at(2).Eta());
	    _floatVars["M_DeltaPhiWlvJet3"] =  fabs(Jet.at(2).DeltaPhi(Wlv));
	    _floatVars["M_DeltaRWlvJet3"] =  fabs(Jet.at(2).DeltaR(Wlv));
	    

	    _floatVars["M_Jet2Jet3_Centrality"] = (Jet.at(1).Pt() + Jet.at(2).Pt())/ (Jet.at(1).P() + Jet.at(2).P());
	    _floatVars["M_Jet1Jet3_Centrality"] = (Jet.at(0).Pt() + Jet.at(2).Pt())/ (Jet.at(0).P() + Jet.at(2).P());
	    _floatVars["M_Jet1Jet2Jet3_Centrality"] = (Jet.at(0).Pt() + Jet.at(1).Pt()+ Jet.at(2).Pt())/ (Jet.at(0).P() + Jet.at(1).P()+ Jet.at(2).P());

	    _floatVars["M_HT"] = (Jet.at(0).Pt() + Jet.at(1).Pt()+ Jet.at(2).Pt());
	  }
      }
  }

  if (evtObj->jets.size() == 3){ //It really should always be 3 by this point, but just to make sure.
    int selectedBJetIndex = -1;
    _floatVars["largestCSV"] = -2.;
    for (unsigned int i = 0; i < 3; i++){
      if (evtObj->jets[i].GetbDiscriminator() > _floatVars["largestCSV"]){
	_floatVars["largestCSV"] = evtObj->jets[i].GetbDiscriminator();
	selectedBJetIndex = i;
      }
    } //end loop finding largest b-tag
    _floatVars["firstJetCSV"] = evtObj->jets[0].GetbDiscriminator();
    _floatVars["secondJetCSV"] = -1.;
    _floatVars["thirdJetCSV"] = -1.;
    if (evtObj->jets.size() > 1) _floatVars["secondJetCSV"] = evtObj->jets[1].GetbDiscriminator();
    if (evtObj->jets.size() > 2) _floatVars["thirdJetCSV"] = evtObj->jets[2].GetbDiscriminator();
    
    _floatVars["taggedJetCSV"] = -1;
    _floatVars["lightJet1CSV"] = -1;
    _floatVars["lightJet2CSV"] = -1;
    if (evtObj->taggedJets.size() > 0) _floatVars["taggedJetCSV"] = evtObj->taggedJets[0].GetbDiscriminator();
    if (evtObj->unTaggedJets.size() > 0) _floatVars["lightJet1CSV"] = evtObj->unTaggedJets[0].GetbDiscriminator();
    if (evtObj->unTaggedJets.size() > 0) _floatVars["lightJet2CSV"] = evtObj->unTaggedJets[1].GetbDiscriminator();

    std::vector<TLorentzVector> lightJets;
    TLorentzVector W(0,0,0,0), Top(0,0,0,0);
    for (unsigned int i = 0; i < 3; i++){
      if (i == selectedBJetIndex) continue; //this is the b jet
      else {
	W = W + Jet.at(i);
	TLorentzVector lighttemp(0,0,0,0);
	lighttemp.SetPtEtaPhiE(Jet.at(i).Pt(),Jet.at(i).Eta(),Jet.at(i).Phi(),Jet.at(i).E());
	lightJets.push_back(lighttemp);
      }
    } //end second loop over jets to add up lightjets
    
    for(int ilightjet=0; ilightjet<lightJets.size();ilightjet++){
      if(fabs(Jet.at(selectedBJetIndex).DeltaR(lightJets.at(ilightjet)))>_floatVars["M_maxDeltaRBJetLightJet"]){
	_floatVars["M_maxDeltaRBJetLightJet"]=fabs(Jet.at(selectedBJetIndex).DeltaR(lightJets.at(ilightjet)));
      }
    }

    _floatVars["M_hadronicWmass"] = W.M();
    _floatVars["M_hadronicWPt"] = W.Pt();
    _floatVars["M_hadronicWEta"] = W.Eta();
    _floatVars["M_hadronicWPhi"] = W.Phi();
    _floatVars["M_DeltaRBJethadronicW"] = fabs(Jet.at(selectedBJetIndex).DeltaR(W));
    _floatVars["M_DeltaRBJetLepton"] = fabs(Jet.at(selectedBJetIndex).DeltaR(Lepton));
    _floatVars["M_DeltaRlightjets"] = fabs(lightJets.at(0).DeltaR(lightJets.at(1)));

    _floatVars["M_cosThetaStar_lepOnly"] = -5.;
    _floatVars["M_cosThetaStar_hadOnly"] = -5.;

    if(Jet.at(selectedBJetIndex).DeltaR(Lepton)>Jet.at(selectedBJetIndex).DeltaR(W)) {
      //Top is hadronic
      
      //Do the wrong top things first
      Top = Lepton + Jet.at(selectedBJetIndex) + Miss;
      _floatVars["M_topMass2_lep"] = -1;
      _floatVars["M_topMassLep_had"] = Top.M();
      //Then apply the correct one.
      Top = W + Jet.at(selectedBJetIndex);
      _floatVars["M_topMass2"] = Top.M();
      _floatVars["M_topMassLep"] = -1;
    }
    else { //Top is leptonic
      Top = W + Jet.at(selectedBJetIndex); //fill the hadronic ones first
      _floatVars["M_topMass2_lep"] = Top.M();
      _floatVars["M_topMassLep_had"] = -1;

      //Now fill the leptonic top variables
      Top = Lepton + Jet.at(selectedBJetIndex) + Miss;
      _floatVars["M_topMass2"] = -1;
      _floatVars["M_topMassLep"] = Top.M();
    }
    TLorentzVector W_topRF = W;
    TLorentzVector lep_wRF = Lepton;
    
    W_topRF.Boost( TVector3(-1.0*Top.BoostVector().Px(), -1.0*Top.BoostVector().Py(), -1.0*Top.BoostVector().Pz()));
    lep_wRF.Boost( TVector3(-1.0*W.BoostVector().Px(), -1.0*W.BoostVector().Py(), -1.0*W.BoostVector().Pz()));
    _floatVars["M_cosThetaStar"] = cos(lep_wRF.Vect().Angle(W_topRF.Vect()));
    if(Jet.at(selectedBJetIndex).DeltaR(Lepton)>Jet.at(selectedBJetIndex).DeltaR(W)) _floatVars["M_cosThetaStar_hadOnly"] = cos(lep_wRF.Vect().Angle(W_topRF.Vect()));
    else _floatVars["M_cosThetaStar_lepOnly"] = cos(lep_wRF.Vect().Angle(W_topRF.Vect()));
    
  } //end 3 jet loop

  if (DoHists()) FillHistograms(evtObj->GetEventWeight());

}
