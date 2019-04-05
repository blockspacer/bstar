/******************************************************************************
 * Muon.cpp                                                                   *
 *                                                                            *
 * Store information about final state muons                                  *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of Muon class                                      *
 *    Muon()                             -- Default Constructor               *
 *    Muon()                             -- Parameterized Constructor         *
 *    Muon()                             -- Constructor with Particle         *
 *    ~Muon()                            -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill                               -- Fill muon from event tree         *
 *    FillFastSim                        -- Fill muon from FastSim tree       *
 *    SetEtCone                          -- Set ET for EM layer 0             *
 *    GetEtCone                          -- Set ET for EM layer 0             *
 *    EtCone                             -- Set ET for EM layer 0             *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const Muon             *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of Muon class                                         *
 *    Double_t _EtCone20                     -- CAL Energy in cone around muon*
 *    Double_t _EtCone30                     -- CAL Energy in cone around muon*
 *                                                                            *
 * History                                                                    *
 *      16 June 2015 - Created by Huaqiao ZHANG                               *
 *****************************************************************************/
#include "SingleTopRootAnalysis/Particles/Recon/Jet.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Muon.hpp"
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

// Integrate classes into the Root system
ClassImp(Muon)

/******************************************************************************
 * Muon::Muon()                                                               *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  Muon::Muon() : Particle::Particle(), 
  _passTightId(kFALSE), 
  _passLooseId(kFALSE),
  _isSoft     (kFALSE),
  _isHighPt   (kFALSE),
  _isPf       (kFALSE),
  _isoCharged (0.0),
  _isoSum     (0.0),
  _isoCharParPt		(0.0),
  _isoNeutralHadron	(0.0),
  _isoPhoton		(0.0),
  _isoPU		(0.0),
  _chi2			(0.0),
  _dxy			(0.0),
  _dz			(0.0),
  _validHits		(0.0),
  _validHitsInner	(0.0),
  _matchedStat		(0.0),
  _TLayers		(0.0),
  _relIsoR04		(0.0),
  _ndof			(0.0),
  _charge		(0.0)
{

} //Muon

/******************************************************************************
 * Muon::~Muon()                                                              *
 *                                                                            *
 * Destructor                                                                 *  
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
Muon::~Muon()
{
  //cout<<"In muon destructor"<<endl;
} //~Muon

/******************************************************************************
 * Muon::Muon(const Muon& other)                                              *
 *                                                                            *
 * Copy constructor Muon const                                                *
 *                                                                            *
 * Input:  Muon                                                               *
 * Output: None                                                               *
 ******************************************************************************/
Muon::Muon(const Muon& other): Particle(other),
  _passTightId(other.GetpassTightId()),
  _passLooseId(other.GetpassLooseId()),
  _isSoft(other.GetisSoft()),
  _isHighPt(other.GetisHighPt()),
  _isPf(other.GetisPf()),
  _isoCharged(other.GetisoCharged()),
  _isoSum(other.GetisoSum()),
  _isoCharParPt(other.GetisoCharParPt()),
  _isoNeutralHadron(other.GetisoNeutralHadron()),
  _isoPhoton(other.GetisoPhoton()),
  _isoPU(other.GetisoPU()),
  _chi2(other.Getchi2()),
  _dxy(other.Getdxy()),
  _dz(other.Getdz()),
  _validHits(other.GetvalidHits()),
  _validHitsInner(other.GetvalidHitsInner()),
  _matchedStat(other.GetmatchedStat()),
  _TLayers(other.GetTLayers()),
  _relIsoR04(other.GetrelIsoR04()),
			       _ndof(other.Getndof()),
			       _charge(other.GetCharge())			       
			       
{
} //Muon()

/******************************************************************************
 * Muon::Muon(const Muon& other)                                              *
 *                                                                            *
 * Copy constructor Particle const                                            *
 *                                                                            *
 * Input:  Muon                                                               *
 * Output: None                                                               *
 ******************************************************************************/
Muon::Muon(const Particle& other): Particle(other),
  _passTightId(kFALSE), 
  _passLooseId(kFALSE),
  _isSoft     (kFALSE),
  _isHighPt   (kFALSE),
  _isPf       (kFALSE),
  _isoCharged (0.0),
  _isoSum     (0.0),
  _isoCharParPt		(0.0),
  _isoNeutralHadron	(0.0),
  _isoPhoton		(0.0),
  _isoPU		(0.0),
  _chi2			(0.0),
  _dxy			(0.0),
  _dz			(0.0),
  _validHits		(0.0),
  _validHitsInner	(0.0),
  _matchedStat		(0.0),
  _TLayers		(0.0),
				   _relIsoR04            (0.0),
				   _ndof (0.0),
				   _charge(0.0)
{
} //Muon

/******************************************************************************
 * Muon& Muon::operator+=(const Muon& other)                                  *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  Muon                                                               *
 * Output: Muon                                                               *
 ******************************************************************************/
Muon& Muon::operator+=(const Muon& other)
{
  
  Particle::operator+=(other);
  return *this;
}

/******************************************************************************
 * Muon& Muon::operator+(const Muon& other)                                   *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  Muon                                                               *
 * Output: Muon                                                               *
 ******************************************************************************/
Muon Muon::operator+(const Muon& other)
{

  Particle particleTemp = Particle::operator+(other);
  return other;
} //+

/******************************************************************************
 * Muon& Muon::operator=(const Muon& other)                                   *
 *                                                                            *
 * Overload = Particle                                                        *
 *                                                                            *
 * Input:  Particle                                                           *
 * Output: Muon                                                               *
 ******************************************************************************/
Muon& Muon::operator=(const Particle& other)
{
  Particle::operator=(other);
  SetpassTightId(kFALSE);
  SetpassLooseId(kFALSE);
  SetisSoft     (kFALSE);
  SetisHighPt   (kFALSE);
  SetisPf       (kFALSE);
  SetisoCharged (0.0);
  SetisoSum     (0.0);
  SetisoCharParPt	(0.0);
  SetisoNeutralHadron	(0.0);
  SetisoPhoton		(0.0);
  SetisoPU		(0.0);
  Setchi2		(0.0);
  Setdxy		(0.0);
  Setdz			(0.0);
  SetvalidHits		(0.0);
  SetvalidHitsInner	(0.0);
  SetmatchedStat	(0.0);
  SetTLayers		(0.0);
  SetrelIsoR04		(0.0);
  Setndof		(0.0);
  SetCharge		(0.0);
  return *this;
} //= Particle

/******************************************************************************
 * Muon& Muon::operator=(const Muon& other)                                   *
 *                                                                            *
 * Overload = const Muon                                                      *
 *                                                                            *
 * Input:  Muon                                                               *
 * Output: Muon                                                               *
 ******************************************************************************/
Muon& Muon::operator=(const Muon& other)
{
  Particle::operator=(other);
  SetpassTightId(other.GetpassTightId());
  SetpassLooseId(other.GetpassLooseId());
  SetisSoft(other.GetisSoft());
  SetisHighPt(other.GetisHighPt());
  SetisPf(other.GetisPf());
  SetisoCharged(other.GetisoCharged());
  SetisoSum(other.GetisoSum());
  SetisoCharParPt(other.GetisoCharParPt());
  SetisoNeutralHadron(other.GetisoNeutralHadron());
  SetisoPhoton(other.GetisoPhoton());
  SetisoPU(other.GetisoPU());
  Setchi2(other.Getchi2());
  Setdxy(other.Getdxy());
  Setdz(other.Getdz());
  SetvalidHits(other.GetvalidHits());
  SetvalidHitsInner(other.GetvalidHitsInner());
  SetmatchedStat(other.GetmatchedStat());
  SetTLayers(other.GetTLayers());
  SetrelIsoR04(other.GetrelIsoR04());
  Setndof(other.Getndof());
  SetCharge(other.GetCharge());
  return *this;
} //= const muon

/******************************************************************************
 * Muon& Muon::operator=(const Muon& other)                                   *
 *                                                                            *
 * Overload = non-const Muon                                                  *
 *                                                                            *
 * Input:  Muon                                                               *
 * Output: Muon                                                               *
 ******************************************************************************/
Muon& Muon::operator=(Muon& other)
{
  Particle::operator=(other);
  SetpassTightId(other.GetpassTightId());
  SetpassLooseId(other.GetpassLooseId());
  SetisSoft(other.GetisSoft());
  SetisHighPt(other.GetisHighPt());
  SetisPf(other.GetisPf());
  SetisoCharged(other.GetisoCharged());
  SetisoSum(other.GetisoSum());
  SetisoCharParPt(other.GetisoCharParPt());
  SetisoNeutralHadron(other.GetisoNeutralHadron());
  SetisoPhoton(other.GetisoPhoton());
  SetisoPU(other.GetisoPU());
  Setchi2(other.Getchi2());
  Setdxy(other.Getdxy());
  Setdz(other.Getdz());
  SetvalidHits(other.GetvalidHits());
  SetvalidHitsInner(other.GetvalidHitsInner());
  SetmatchedStat(other.GetmatchedStat());
  SetTLayers(other.GetTLayers());
  SetrelIsoR04(other.GetrelIsoR04());
  Setndof(other.Getndof());
  SetCharge(other.GetCharge());
  return *this;
} //= non-const muon

/******************************************************************************         
 * void Muon::SetCuts(TEnv* config, TString muonType)                         * 
 *                                                                            *         
 * Set up the cuts for muonType                                               *
 *                                                                            *         
 * Input:  TEnv* config, TString muonType                                     *
 * Output: None                                                               *
 ******************************************************************************/
void Muon::SetCuts(TEnv* config, TString muonType)
{
	_minPtCuts[muonType] = config -> GetValue("ObjectID.Muon."+muonType+".MinPt", 100.0);
	_maxEtaCuts[muonType] = config -> GetValue("ObjectID.Muon."+muonType+".MaxEta", 0.0);
	_maxRelIsoCuts[muonType] = config -> GetValue("ObjectID.Muon."+muonType+".MaxRelIso", 100.0);
	_minLeptonJetDetaR = config -> GetValue("ObjectID.Lepton.ClostJetDeltaRMin", 100.0);
	_relPtcut = config -> GetValue("ObjectID.Lepton.RelPt", 100.0);
	_TT_CR =              config -> GetValue("TT_CR",0);
	_QCD_CR =             config -> GetValue("QCD_CR",0);
	_bstar = config -> GetValue("ifbstar",0);
}

/******************************************************************************         
 * void Muon::Fill(EventTree *evtr, Int_t iE)                                 *         
 *                                                                            *         
 * Fill Muon vector from tree                                                 *         
 *                                                                            *         
 * Input:  Event Tree                                                         *         
 * Output: kTRUE if the muon passes object ID cuts                            *         
 ******************************************************************************/
Bool_t Muon::Fill(EventTree *evtr,std::vector<Jet>& selectedjets,int iE,TString muonType, Bool_t isSimulation)
{
  // **************************************************************
  // Check muon type
  // **************************************************************
  if( (muonType != "All") && (muonType != "Veto") && (muonType != "Tight") && (muonType != "Isolated") && (muonType != "UnIsolated") && (muonType != "PtEtaCut") ){
    std::cout << "ERROR: <Muon::Fill()> " << "Passed variable muonType of value " << muonType << " is not valid.  "
	      << "Must be All, Veto, Tight, Isolated, UnIsolated, PtEtaCut" << std::endl;
  } //if
  
  // **************************************************************
  // Fill muon
  // **************************************************************
  Double_t muPx     = evtr -> Muon_px       -> operator[](iE);
  Double_t muPy     = evtr -> Muon_py       -> operator[](iE);
  Double_t muPz     = evtr -> Muon_pz       -> operator[](iE);
  Double_t muPt     = evtr -> Muon_pt       -> operator[](iE);
  Double_t muEta    = evtr -> Muon_eta      -> operator[](iE);
  Double_t muPhi    = evtr -> Muon_phi      -> operator[](iE);
  Double_t muE      = evtr -> Muon_energy   -> operator[](iE);
  Double_t muCharge = evtr -> Muon_charge   -> operator[](iE);
  //Double_t murelpt= evtr -> Muon_ptrel-> operator[](iE);
  int trigger       = evtr -> HLT_Mu50    ; 
  SetPtEtaPhiE(muPt, muEta, muPhi, muE);

  SetpassTightId	(evtr -> Muon_tight   		-> operator[](iE));
  SetpassLooseId	(evtr -> Muon_loose   		-> operator[](iE));
  SetisSoft     	(evtr -> Muon_soft    		-> operator[](iE));
  SetisHighPt   	(evtr -> Muon_isHighPt		-> operator[](iE));
  SetisPf       	(evtr -> Muon_pf      		-> operator[](iE));
  SetisoCharged 	(evtr -> Muon_isoR04Charged 	-> operator[](iE));
  SetisoSum     	(evtr -> Muon_isoSum  		-> operator[](iE));
  SetisoCharParPt	(evtr -> Muon_isoR04CharParPt  	-> operator[](iE));
  SetisoNeutralHadron	(evtr -> Muon_isoR04NeutralHadron  -> operator[](iE));
  SetisoPhoton		(evtr -> Muon_isoR04Photon   	-> operator[](iE));
  SetisoPU		(evtr -> Muon_isoR04PU   		-> operator[](iE));
  Setchi2		(evtr -> Muon_chi2   		-> operator[](iE));
  Setdxy		(evtr -> Muon_dxy_pv   		-> operator[](iE));
  Setdz			(evtr -> Muon_dz_pv   		-> operator[](iE));
  SetvalidHits		(evtr -> Muon_validHits   	-> operator[](iE));
  SetvalidHitsInner	(evtr -> Muon_validHitsInner   	-> operator[](iE));
  SetmatchedStat	(evtr -> Muon_matchedStat   	-> operator[](iE));
  SetTLayers		(evtr -> Muon_TLayers   	-> operator[](iE));
  SetrelIsoR04		(evtr -> Muon_relIsoDeltaBetaR04-> operator[](iE));
  Setndof		(evtr -> Muon_ndof		-> operator[](iE));
  SetCharge		(evtr -> Muon_charge		-> operator[](iE));
 






  // **************************************************************
  // Run 2 relative isolation cuts
  // **************************************************************
  Bool_t passRelIso = kTRUE;

  //  std::cout << relIsoR04() << " " << _maxRelIsoCuts[muonType] << " " << muonType << " ";

  if (relIsoR04() > _maxRelIsoCuts[muonType]) {passRelIso = kFALSE;
    //    std::cout << "false";
  }
  //  else std::cout << "true";

   // std::cout << " " << passRelIso << std::endl;
  
  // **************************************************************
  // Pt and Eta Cuts
  // **************************************************************
  // If event passes or fails requirements
  Double_t closestJetDetaR = 999.;
  Double_t murelpt= 999.;
  Bool_t passMinPt  = kTRUE;
  Bool_t passMaxEta = kTRUE;

  Bool_t passCustomID = kTRUE;
  
  Bool_t passtrigger= kTRUE;
  Bool_t passMuonJet2Dcut = kFALSE;
  Bool_t passMuonclosetJetDeatR = kTRUE;
  TVector3 Muon(muPx,muPy,muPz);
  Double_t muonMag = 0.0;
   Double_t muonAngle = 0.0 ;
   muonMag = Muon.Mag();
  for (auto const & jet : selectedjets){
	      if (jet.Pt()>30&&jet.DeltaR(*this) < closestJetDetaR) {
		      closestJetDetaR = jet.DeltaR(*this);
		      muonAngle = jet.Angle(Muon);
		      murelpt=muonMag*sin(muonAngle);
		      }
	        }
	if (closestJetDetaR<_minLeptonJetDetaR) passMuonclosetJetDeatR = kFALSE;
  if (_QCD_CR) {
	  if(passMuonclosetJetDeatR||murelpt<_relPtcut)passMuonJet2Dcut =kTRUE;
	  }
  else {
	 if(passMuonclosetJetDeatR||murelpt>_relPtcut)passMuonJet2Dcut =kTRUE;
	}

if (chi2()/ndof() 	>= 10  ||
      TLayers()		<= 5   ||
      validHits()	<  1   ||
      dxy()		>= 0.2 ||
      dz()		>= 0.5 ||
      validHitsInner()	<  1   ||
      matchedStat()	<  2)
    passCustomID = kFALSE;
  
  // Test Requirements
  if(muPt <= _minPtCuts[muonType])               passMinPt  = kFALSE;
  if(trigger == 0)                		 passtrigger  = kFALSE;
   //cout<<"if it pass Pt cut :"<<passMinPt<<endl;
  if(TMath::Abs(muEta) >= _maxEtaCuts[muonType]) passMaxEta = kFALSE;
  
  
  
  
  if(_bstar == 1){
	//if(     "Tight"      == muonType) return (passMinPt && passMaxEta  && passTightId() &&passRelIso);
	if(     "Tight"      == muonType) return (passMinPt && passMaxEta && passTightId());
        else if("Veto"       == muonType) return (passMinPt && passMaxEta);//no isolation req. or inner det or jet overlap.
   	else if("UnIsolated" == muonType) return (passMinPt && passMaxEta  && passTightId() && passCustomID && ! passRelIso); //The same as tight muons, but with an inverted isolation requirement
}
  else{ 

  //  if(     "Tight"      == muonType) return( passMinPt && passMaxEta  && IsTight() && Isolation() && !GetOverlapWithJet() && IsCombinedMuon());
  	if(     "Tight"      == muonType) return (passMinPt && passMaxEta  && passTightId() && passCustomID && passRelIso);
  //if(     "Tight"      == muonType) return (passMinPt && passMaxEta  && passTightId() );
  	else if("Veto"       == muonType) return (passMinPt && passMaxEta);//no isolation req. or inner det or jet overlap.
 	else if("UnIsolated" == muonType) return (passMinPt && passMaxEta  && passTightId() && passCustomID && ! passRelIso); //The same as tight muons, but with an inverted isolation requirement
    //    std::cout << muPt << " " << minPt << " " << muEta << " " << maxEta << std::end;
  //else if("Isolated"   == muonType) return( GetIsolation()  && !GetOverlapWithJet() && IsCombinedMuon() && OverlapUse());
  //else if("UnIsolated" == muonType) return( !GetIsolation()  && passMinPt && passMaxEta && IsTight() && !GetOverlapWithJet()&& IsCombinedMuon());
  //else if("All"        == muonType) return( kTRUE );
  //else return kTRUE;
}
  return kTRUE;

} //Fill()
 
/******************************************************************************         
 * void Muon::FillFastSim(FastSimTree *tr, Int_t iE)                          *         
 *                                                                            *         
 * Fill Muon vector from tree                                                 *         
 *                                                                            *         
 * Input:  FastSim Tree                                                       *         
 * Output: kTRUE if the muon passes object ID cuts                            *         
 ******************************************************************************/
Bool_t Muon::FillFastSim(TEnv *config,FastSimTree *tr,int iE,TString muonType)
{

} //FillFastSim()


