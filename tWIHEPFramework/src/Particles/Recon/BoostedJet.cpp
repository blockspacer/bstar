/******************************************************************************
 * Jet.cpp                                                                    *
 *                                                                            *
 * Store information about final state jets                                   *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of Jet class                                       *
 *    Jet()                              -- Default Constructor               *
 *    Jet()                              -- Parameterized Constructor         *
 *    ~Jet()                             -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill                               -- Fill jet from event tree          *
 *    FillFastSim                        -- Fill jet from FastSim tree        *
 *    SetEtM0                            -- Set ET for EM layer 0             *
 *    GetEtM0                            -- Set ET for EM layer 0             *
 *    EtEM0                              -- Set ET for EM layer 0             *
 *    ApplyJESCorrection                 -- Shift jet energy by some amount   *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const Jet              *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of Jet class                                          *
 *    Double_t _EtEm0                     -- Et for ME layer 0                *
 *                                                                            *
 * History                                                                    *
 *      16 Nov 2006 - Created by R. Schwienhorst                              *
 *      17 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      20 Mar 2007 - RS: Add filling from fastsim tree                       *
 *      29 Nov 2007 - RS: Add b-tagging choice when filling                   *
 *      11 Jan 2008 - RS: Fill only if object passes ID cuts                  *
 *      14 Jan 2008 - RS: For TRFIP2D, add systematics                        *
 *      08 May 2009 - J. Holzbauer: added quarkFlavorLabel tag                *
 *****************************************************************************/
#include "SingleTopRootAnalysis/Particles/Recon/BoostedJet.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Electron.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Muon.hpp"

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <TRandom3.h>
#include<TString.h>


using namespace std;
//using namespace Analysis;

//
// give me a random number generator.
// put it here and make it static so that I don't create multiple copies
static TRandom3 boostedjet_myran;
static TRandom3 boostedjet_jer_myran;

// Integrate classes into the Root system
ClassImp(BoostedJet)

/******************************************************************************
 * Jet::Jet()                                                                 *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  BoostedJet::BoostedJet() : Particle::Particle(),
  _prunedmass(0.0),_tau1(0.0),_tau2(0.0),_tau3(0.0),_numberOfConstituents(0), _chargedMultiplicity(0),  _bDiscriminator ( -999.0), _mass ( 0.0), _uncorrPt ( 0.0), _neutralHadEnergyFraction(0.0), _neutralEmEmEnergyFraction ( 0.0), _chargedHadronEnergyFraction (0.0), _chargedEmEnergyFraction(0.0), _muonEnergyFraction(0.0), _electronEnergy(0.0), _photonEnergy(0.0), _jesUp(false), _jesDown(false), _jerUp(false), _jerDown(false), _hadronFlavour(-1), _tagged(0)
{
} //Jet()

/******************************************************************************
 * Jet::~Jet()                                                                *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
 BoostedJet::~BoostedJet()
{
} //~Jet()

/******************************************************************************
 * Jet::Jet(const Jet& other)                                                 *
 *                                                                            *
 * Copy Constructor Jet const                                                 *
 *                                                                            *
 * Input:  Jet                                                                *
 * Output: None                                                               *
 ******************************************************************************/
BoostedJet::BoostedJet(const BoostedJet& other): Particle(other),
_prunedmass			(other.Getprunedmass()),
_tau1				(other.Gettau1()),
_tau2				(other.Gettau2()),
_tau3				(other.Gettau3()),
_numberOfConstituents		(other.GetnumberOfConstituents()), 
_hadronFlavour		        (other.GethadronFlavour()), 
_tagged				(other.IsTagged()),
_chargedMultiplicity		(other.GetchargedMultiplicity()),  
_bDiscriminator 		(other.GetbDiscriminator()), 
_mass 				(other.Getmass()), 
_uncorrPt 			(other.GetuncorrPt()), 
_neutralHadEnergyFraction	(other.GetneutralHadEnergyFraction()), 
_neutralEmEmEnergyFraction 	(other.GetneutralEmEmEnergyFraction()),
 _chargedHadronEnergyFraction 	(other.GetchargedHadronEnergyFraction()), 
_chargedEmEnergyFraction	(other.GetchargedEmEnergyFraction()), 
_muonEnergyFraction		(other.GetmuonEnergyFraction()), 
_electronEnergy			(other.GetelectronEnergy()), 
_photonEnergy			(other.GetphotonEnergy())
{
} //Jet()

/******************************************************************************
 * Jet::Jet(const Particle& other)                                            *
 *                                                                            *
 * Copy Constructor Particle const                                            *
 *                                                                            *
 * Input:  Particle                                                           *
 * Output: None                                                               *
 ******************************************************************************/
BoostedJet::BoostedJet(const Particle& other): Particle(other),
				 _prunedmass(0.0),_tau1(0.0),_tau2(0.0),_tau3(0.0),_numberOfConstituents(0), _hadronFlavour(-1), _chargedMultiplicity(0),  _bDiscriminator ( -999.0),  _mass ( 0.0), _uncorrPt ( 0.0), _neutralHadEnergyFraction(0.0), _neutralEmEmEnergyFraction ( 0.0), _chargedHadronEnergyFraction (0.0), _chargedEmEnergyFraction(0.0), _muonEnergyFraction(0.0), _electronEnergy(0.0), _photonEnergy(0.0), _tagged(0)
{
 
} //Jet()


/******************************************************************************
 * Jet& Jet::operator+=(const Jet& other)                                     *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  Jet                                                                *
 * Output: Jet                                                                *
 ******************************************************************************/
BoostedJet& BoostedJet::operator+=(const BoostedJet& other)
{
  
  Particle::operator+=(other);
  return *this;
} //+=

/******************************************************************************
 * Jet& Jet::operator+(const Jet& other)                                      *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  Jet                                                                *
 * Output: Jet                                                                *
 ******************************************************************************/
BoostedJet BoostedJet::operator+(const BoostedJet& other)
{
  Particle particleTemp = Particle::operator+(other);
  BoostedJet boostedjetTemp = particleTemp;
  return boostedjetTemp;

} //+

/******************************************************************************
 * Jet& Jet::operator=(const Jet& other)                                      *
 *                                                                            *
 * Overload = Particle                                                        *
 *                                                                            *
 * Input:  Particle                                                           *
 * Output: Jet                                                                *
 ******************************************************************************/
BoostedJet& BoostedJet::operator=(const Particle& other)
{
  
  Particle::operator=(other);
  Setprunedmass(0.0),
  Settau1(0.0),
  Settau2(0.0),
  Settau3(0.0),
  SetnumberOfConstituents(0), 
  SethadronFlavour(-1),
  SetTagged(0),
  SetchargedMultiplicity(0),  
  SetbDiscriminator ( -999.0), 
  Setmass ( 0.0), 
  SetuncorrPt ( 0.0), 
  SetneutralHadEnergyFraction(0.0), 
  SetneutralEmEmEnergyFraction ( 0.0), 
  SetchargedHadronEnergyFraction (0.0), 
  SetchargedEmEnergyFraction(0.0), 
  SetmuonEnergyFraction(0.0), 
  SetelectronEnergy(0.0), 
  SetphotonEnergy(0.0);

  return *this;
} //= Particle

/******************************************************************************
 * Jet& Jet::operator=(const Jet& other)                                      *
 *                                                                            *
 * Overload = const Jet                                                       *
 *                                                                            *
 * Input:  Jet                                                                *
 * Output: Jet                                                                *
 ******************************************************************************/
BoostedJet& BoostedJet::operator=(const BoostedJet& other)
{
  
  Particle::operator=(other);
  Setprunedmass				(other.Getprunedmass());
  Settau1				(other.Gettau1());
  Settau2				(other.Gettau2());
  Settau3				(other.Gettau3());
  SetnumberOfConstituents		(other.GetnumberOfConstituents());
  SethadronFlavour	         	(other.GethadronFlavour());
  SetTagged				(other.IsTagged());
  SetchargedMultiplicity		(other.GetchargedMultiplicity()); 
  SetbDiscriminator 			(other.GetbDiscriminator());
  Setmass 				(other.Getmass());
  SetuncorrPt 				(other.GetuncorrPt());
  SetneutralHadEnergyFraction		(other.GetneutralHadEnergyFraction());
  SetneutralEmEmEnergyFraction 		(other.GetneutralEmEmEnergyFraction());
  SetchargedHadronEnergyFraction 	(other.GetchargedHadronEnergyFraction());
  SetchargedEmEnergyFraction		(other.GetchargedEmEnergyFraction());
  SetmuonEnergyFraction			(other.GetmuonEnergyFraction());
  SetelectronEnergy			(other.GetelectronEnergy());
  SetphotonEnergy			(other.GetphotonEnergy());
  return *this;
} //= const

/******************************************************************************
 * Jet& Jet::operator=(const Jet& other)                                      *
 *                                                                            *
 * Overload = non-const Jet                                                   *
 *                                                                            *
 * Input:  Jet                                                                *
 * Output: Jet                                                                *
 ******************************************************************************/
BoostedJet& BoostedJet::operator=(BoostedJet& other)
{
  Particle::operator=(other);
  Setprunedmass                         (other.Getprunedmass());
  Settau1                               (other.Gettau1());
  Settau2                               (other.Gettau2());
  Settau3                               (other.Gettau3());

  SetnumberOfConstituents		(other.GetnumberOfConstituents());
  SethadronFlavour	         	(other.GethadronFlavour());
  SetTagged			        (other.IsTagged());
  SetchargedMultiplicity		(other.GetchargedMultiplicity()); 
  SetbDiscriminator 			(other.GetbDiscriminator());
  Setmass 				(other.Getmass());
  SetuncorrPt 				(other.GetuncorrPt());
  SetneutralHadEnergyFraction		(other.GetneutralHadEnergyFraction());
  SetneutralEmEmEnergyFraction 		(other.GetneutralEmEmEnergyFraction());
  SetchargedHadronEnergyFraction 	(other.GetchargedHadronEnergyFraction());
  SetchargedEmEnergyFraction		(other.GetchargedEmEnergyFraction());
  SetmuonEnergyFraction			(other.GetmuonEnergyFraction());
  SetelectronEnergy			(other.GetelectronEnergy());
  SetphotonEnergy			(other.GetphotonEnergy());
  return *this;
} //= non-const

/******************************************************************************         
 * void Jet::SetCuts(TEnv* config)                                            *
 *                                                                            *         
 * Set the cuts to be used on the jet objects                                 *
 *                                                                            *         
 * Input:  TEnv* config                                                       * 
 *         - Contains the cuts to be applied to each jet                      *
 * Output: None                                                               *
 ******************************************************************************/

void BoostedJet::SetCuts(TEnv * config)
{
  _maxEtaCut = 			config -> GetValue("ObjectID.BoostedJet.MaxEta",100.);
  _minPtCut = 			config -> GetValue("ObjectID.BoostedJet.MinPt",0.);
  _minCRPtCut = 		config -> GetValue("ObjectID.BoostedJet.CR.MinPt",0.);
  _closestLeptonCut = 		config -> GetValue("ObjectID.BoostedJet.LepCleanR",0.0);
  _jesUp = 			config -> GetValue("Systs.doJESUp",0);
  _jesDown = 			config -> GetValue("Systs.doJESDown",0);
  _jerUp = 			config -> GetValue("Systs.doJERUp",0);
  _jerDown = 			config -> GetValue("Systs.doJERDown",0);
  _bstar = 			config -> GetValue("ifbstar",0);
  _TT_CR = 			config -> GetValue("TT_CR",0);
  _QCD_CR = 			config -> GetValue("QCD_CR",0);
  _minMuonBoostedJetDetaR =     config -> GetValue("ObjectID.BoostedJet.MuonDeltaRMin",0.0);
  _minprunedmass = 		config -> GetValue("ObjectID.BoostedJet.MinPrunedMass",0.0);
  _maxprunedmass = 		config -> GetValue("ObjectID.BoostedJet.MaxPrunedMass",999.9);
  _minsoftdropmass = 		config -> GetValue("ObjectID.BoostedJet.MinSoftdropMass",0.0);
  _maxsoftdropmass = 		config -> GetValue("ObjectID.BoostedJet.MaxSoftdropMass",999.9);
  _tau21 = 			config -> GetValue("ObjectID.BoostedJet.Tau21",0.0);
  _tau32 = 			config -> GetValue("ObjectID.BoostedJet.Tau32",0.0);
  _closestLeptonCut =           config -> GetValue("ObjectID.BoostedJet.LepCleanR",0.0); 


}

/******************************************************************************         
 * void Jet::Fill(EventTree *evtr, Int_t iE,const TString& tagName)           *         
 *                                                                            *         
 * Fill jet vector from tree                                                  *         
 *                                                                            *         
 * Input:  Event Tree pointer                                                 *         
 *         - which jet                                                        *
 * Output: True if this jet passes jet ID cuts                                *         
 ******************************************************************************/

Bool_t BoostedJet::Fill( double myJESCorr, double myJERCorr, std::vector<Muon>& selectedMuons, std::vector<Electron>& selectedElectrons, EventTree *evtr, Int_t iE, TLorentzVector * met, bool isMC)
{

  Double_t boostedjetPt, boostedjetEta,boostedjetPhi,boostedjetE, boostedjetCharge, boostedjetM,prunedmass,softdropmass,tau1,tau2,tau3,boostedjet_Uncorr_pt,boostedjet_JesSF,boostedjet_JerSF,boostedjet_pt;
  boostedjetPt     = evtr -> BoostedJet_pt     -> operator[](iE);
  boostedjetEta    = evtr -> BoostedJet_eta    -> operator[](iE);
  boostedjetPhi    = evtr -> BoostedJet_phi    -> operator[](iE);
  boostedjetE      = evtr -> BoostedJet_energy -> operator[](iE);
  prunedmass       = evtr -> BoostedJet_pruned_mass -> operator[](iE);  
  softdropmass     = evtr -> BoostedJet_softdrop_mass-> operator[](iE);
  tau1     	   = evtr -> BoostedJet_tau1-> operator[](iE);
  tau2     	   = evtr -> BoostedJet_tau2-> operator[](iE);
  tau3     	   = evtr -> BoostedJet_tau3-> operator[](iE);
  boostedjet_Uncorr_pt     = evtr -> BoostedJet_Uncorr_pt     -> operator[](iE);
  boostedjet_JesSF     = evtr -> BoostedJet_JesSF     -> operator[](iE);
  boostedjet_JerSF     = evtr -> BoostedJet_JerSF     -> operator[](iE);
  boostedjet_pt = boostedjet_Uncorr_pt*boostedjet_JesSF*boostedjet_JerSF;

  if(boostedjetE > 0){
    //SetPtEtaPhiE(boostedjetPt,boostedjetEta,boostedjetPhi,boostedjetE);
    SetPtEtaPhiE(boostedjet_pt,boostedjetEta,boostedjetPhi,boostedjetE);
  }else{
    SetPtEtaPhiE(1,1,1,3);
  }
  Setprunedmass                         (evtr -> BoostedJet_pruned_mass -> operator[](iE));
  Settau1				(evtr -> BoostedJet_tau1-> operator[](iE));
  Settau2				(evtr -> BoostedJet_tau2-> operator[](iE));
  Settau3				(evtr -> BoostedJet_tau3-> operator[](iE));
  SetnumberOfConstituents		(evtr -> BoostedJet_numberOfConstituents     	-> operator[](iE));
  SetchargedMultiplicity		(evtr -> BoostedJet_chargedMultiplicity     	-> operator[](iE)); 
  SetbDiscriminator 			(evtr -> BoostedJet_pfCombinedInclusiveSecondaryVertexV2BJetTags     	-> operator[](iE));
  Setmass 				(evtr -> BoostedJet_mass     			-> operator[](iE));
  SetuncorrPt 				(evtr -> BoostedJet_Uncorr_pt     		-> operator[](iE));
  SetneutralHadEnergyFraction		(evtr -> BoostedJet_neutralHadEnergyFraction	-> operator[](iE));
  SetchargedHadronEnergyFraction 	(evtr -> BoostedJet_chargedHadronEnergyFraction-> operator[](iE));
  SetchargedEmEnergyFraction		(evtr -> BoostedJet_chargedEmEnergyFraction    -> operator[](iE));
  SetmuonEnergyFraction			(evtr -> BoostedJet_muonEnergyFraction     	-> operator[](iE));
  SetelectronEnergy			(evtr -> BoostedJet_electronEnergy     	-> operator[](iE));
  SetphotonEnergy			(evtr -> BoostedJet_photonEnergy     		-> operator[](iE));

  // Now we want to do the JER and JES systematic adjustments to the jet. This also requires correcting the MET.
  if (_jesUp || _jesDown || _jerUp || _jerDown) SystematicPtShift(evtr, iE, met);  


  /////////////////////////////////////////////////////////////////////////////
  // Jet Pt and Eta Cuts
  /////////////////////////////////////////////////////////////////////////////

  Bool_t passPt = boostedjet_pt > _minPtCut;
  Bool_t passEta = TMath::Abs(Eta()) < _maxEtaCut;


  /////////////////////////////////////////////////////////////////////////
  //  // Jet Cleaning
  //    /////////////////////////////////////////////////////////////////////////
  //    //  cout<<"we are going to jet cleaning "<<endl;
  //
  Bool_t passesCleaning = kTRUE;

  Bool_t selectMuon = kFALSE;
  Double_t closestLepton = 999.;

  for (auto const & ele : selectedElectrons){
	  if (ele.DeltaR(*this) < closestLepton) closestLepton = ele.DeltaR(*this);
  }
  for (auto const & mu : selectedMuons){
	  selectMuon = kTRUE;
	  if (mu.DeltaR(*this) < closestLepton) closestLepton = mu.DeltaR(*this);
  }
  if (closestLepton < _closestLeptonCut) passesCleaning = kFALSE;
//  SetClosestLep(closestLepton);





 /////////////////////////////////////////////////////////////////////////
 //   overlap related cuts
 ///////////////////////////////////////////////////////////////////////////
 Bool_t passMuonBoostedJetDetaR = kFALSE;
 Double_t muonDeltaR = 0.0;
 std::vector<Muon>::iterator muonIterator; 
 for(muonIterator = selectedMuons.begin(); muonIterator != selectedMuons.end(); muonIterator++){
    muonDeltaR = muonIterator->DeltaR(*this);
    if(muonDeltaR > _minMuonBoostedJetDetaR){
	passMuonBoostedJetDetaR = kTRUE ;
	break;
   }
// cout<<"if is pass  "<<passMuonJetDetaR<<endl<<"the mag of detaR = "<<muonDeltaR<<endl;

  }
Double_t tau21 = tau2/tau1;
Double_t tau32 = tau3/tau2;
//cout<<"tau21 = "<<tau21 <<endl;
Bool_t passmasswindow;
if(_TT_CR){
	passmasswindow =_minsoftdropmass<softdropmass && softdropmass< _maxsoftdropmass;}
else{	
	passmasswindow =_minprunedmass <prunedmass && prunedmass< _maxprunedmass;}
Bool_t passtau;
if(_TT_CR){
	passtau = tau32<_tau32;}
else if(_QCD_CR){
	passtau =tau21>_tau21;}
else{
        passtau	= tau21<_tau21;}


 //if (passPt && passEta  && passMuonBoostedJetDetaR && passmasswindow && passtau){ cout<<"pruned mass :"<<prunedmass<<endl;return kTRUE;}
if (passPt && passEta && passMuonBoostedJetDetaR && passmasswindow && passtau) return kTRUE;

  return kFALSE;

} //Fill()


/******************************************************************************         
 * void Jet::SystematicPtShift(EventTree * evtr)                         *
 *                                                                            *         
 * Apply systematic shifts in jet pt                                          *
 *                                                                            *         
 * Input:  - the event tree (to access systematic SFs and MET info)      *
 * Output: -                                                                  *
 ******************************************************************************/
void BoostedJet::SystematicPtShift(EventTree * evtr, Int_t iE, TLorentzVector * met){


  //  std::cout << "syst correct" << std::endl;
  float ptSF = 1.0;
  if (_jesUp){
    ptSF = evtr->BoostedJet_JesSFup->operator[](iE)/evtr->BoostedJet_JesSF->operator[](iE);
  }
  if (_jesDown){
    ptSF = evtr->BoostedJet_JesSFdown->operator[](iE)/evtr->BoostedJet_JesSF->operator[](iE);
  }
  if (_jerUp){
    ptSF = evtr->BoostedJet_JerSFup->operator[](iE)/evtr->BoostedJet_JerSF->operator[](iE);
  }
  if (_jerDown){
    ptSF = evtr->BoostedJet_JerSFdown->operator[](iE)/evtr->BoostedJet_JerSF->operator[](iE);
  }
  //  float ptBefore = Pt();
  //Remove jet from MET
  met->SetPx(met->Px() + Px());
  met->SetPy(met->Py() + Py());
  //Apply the correction
  //  std::cout << Px() << " " << Py() << " " << Pt() << std::endl;
  SetPx(Px()*ptSF);
  SetPy(Py()*ptSF);
  SetPz(Pz()*ptSF);
  //std::cout << Px() << " " << Py() << " " << Pt() << std::endl << std::endl;
  
  //Propagate to MET
  met->SetPx(met->Px() - Px());
  met->SetPy(met->Py() - Py());

}
