/******************************************************************************
 * Electron.hpp                                                               *
 *                                                                            *
 * Store information about final state electrons                              *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of Electron class                                  *
 *    Electron()                         -- Default Constructor               *
 *    Electron()                         -- Parameterized Constructor         *
 *    Clear()                            -- Set contents to default           *
 *    Fill                               -- Fill electron from event tree     *
 *    FillFastSim                        -- Fill electron from FastSim tree   *
 *    SetEoverP                          -- Set E over P for electron         *
 *    GetEoverP                          -- GetEoverP                         *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const Electron         *
 *    operator =                         -- Overload =                        *
 *                                                                            *
 * Private Data Members of Electron class                                     *
 *                                                                            *
 * History                                                                    *
 *      16 June 2015 - Created by Huaqiao ZHANG for CMS                       *
 *****************************************************************************/

#ifndef electron_h
#define electron_h

#include "SingleTopRootAnalysis/Particles/Recon/Jet.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Particle.hpp"
#include "SingleTopRootAnalysis/Trees/EventTree.hpp"
#include "SingleTopRootAnalysis/Trees/FastSimTree.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/Top_MET.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Muon.hpp"

class Muon;
class Jet;
class Electron: public Particle
{
 public:

  // Default constructor
  Electron();

  // Copy constructor for Electrons
  Electron(const Electron& other);

  // Copy constructor for Particles
  Electron(const Particle& other);

  // Destructor
  ~Electron();

  // Overloaded operators:
  // +=
  Electron& operator+=(const Electron& other);
  // +
  Electron operator+(const Electron& other);
  // = Particle
  Electron& operator=(const Particle& other);
  // = Electron const
  Electron& operator=(const Electron& other);
  // = Electron non-const
  Electron& operator=(Electron& other);

  // Set all contents to their defaults
  inline void Clear() { Particle::Clear(); _passVetoId = 0.0; _passLooseId = 0.0; _passMediumId = 0.0; _passTightId = 0.0; _passHEEPId = 0.0; _passConversionVeto = 0.0;
  _patElectron_d0 = 0.0; _patElectron_dz = 0.0; _expectedMissingInnerHits = 0.0;
  _isoChargedHadrons = 0.0; _isoNeutralHadrons = 0.0; _isoPhotons = 0.0; _isoPU = 0.0;
  }

  //Set up the cuts
  void SetCuts(TEnv* config, TString electronType);

  // Fill the electron from an EventTree
  //Bool_t Fill(EventTree *evtr, Int_t iE, TString electronType, Bool_t isSimulation=false , std::vector<Muon>& selectedMuons);
  Bool_t Fill(EventTree *evtr, Int_t iE, TString electronType, Bool_t isSimulation, std::vector<Muon>& selectedMuons,std::vector<Jet>& selectedjets);
  //Bool_t Fill(EventTree *evtr, Int_t iE, TString electronType, Bool_t isSimulation=false);

  // also fill from a fastsim tree
  Bool_t FillFastSim(FastSimTree *tr, Int_t iE,TEnv* config,TString electronType);

  //inline void SetthetaID(Double_t thetaID){_thetaID = thetaID;};
  //inline Double_t GetthetaID() const {return _thetaID;};
  //inline Double_t thetaID() const {return _thetaID;};

  inline void SetpassVetoId(Int_t passVetoId){_passVetoId = passVetoId;};
  inline Int_t GetpassVetoId() const {return _passVetoId;};
  inline Int_t passVetoId() const {return _passVetoId;};

  inline void SetpassLooseId(Int_t passLooseId){_passLooseId = passLooseId;};
  inline Int_t GetpassLooseId() const {return _passLooseId;};
  inline Int_t passLooseId() const {return _passLooseId;};

  inline void SetpassMediumId(Int_t passMediumId){_passMediumId = passMediumId;};
  inline Int_t GetpassMediumId() const {return _passMediumId;};
  inline Int_t passMediumId() const {return _passMediumId;};

  inline void SetpassTightId(Int_t passTightId){_passTightId = passTightId;};
  inline Int_t GetpassTightId() const {return _passTightId;};
  inline Int_t passTightId() const {return _passTightId;};

  inline void SetpassHEEPId(Int_t passHEEPId){_passHEEPId = passHEEPId;};
  inline Int_t GetpassHEEPId() const {return _passHEEPId;};
  inline Int_t passHEEPId() const {return _passHEEPId;};

  inline void SetpassConversionVeto(Int_t passConversionVeto){_passConversionVeto = passConversionVeto;};
  inline Int_t GetpassConversionVeto() const {return _passConversionVeto;};
  inline Int_t passConversionVeto() const {return _passConversionVeto;};

  inline void SetexpectedMissingInnerHits(Int_t expectedMissingInnerHits){_expectedMissingInnerHits = expectedMissingInnerHits;};
  inline Int_t GetexpectedMissingInnerHits() const {return _expectedMissingInnerHits;};
  inline Int_t expectedMissingInnerHits() const {return _expectedMissingInnerHits;};

  inline void SetpatElectron_d0(Double_t patElectron_d0){_patElectron_d0 = patElectron_d0;};
  inline Double_t GetpatElectron_d0() const {return _patElectron_d0;};
  inline Double_t patElectron_d0() const {return _patElectron_d0;};

  inline void SetpatElectron_dz(Double_t patElectron_dz){_patElectron_dz = patElectron_dz;};
  inline Double_t GetpatElectron_dz() const {return _patElectron_dz;};
  inline Double_t patElectron_dz() const {return _patElectron_dz;};

  inline void SetisoChargedHadrons(Double_t isoChargedHadrons){_isoChargedHadrons = isoChargedHadrons;};
  inline Double_t GetisoChargedHadrons() const {return _isoChargedHadrons;};
  inline Double_t isoChargedHadrons() const {return _isoChargedHadrons;};

  inline void SetisoNeutralHadrons(Double_t isoNeutralHadrons){_isoNeutralHadrons = isoNeutralHadrons;};
  inline Double_t GetisoNeutralHadrons() const {return _isoNeutralHadrons;};
  inline Double_t isoNeutralHadrons() const {return _isoNeutralHadrons;};

  inline void SetisoPhotons(Double_t isoPhotons){_isoPhotons = isoPhotons;};
  inline Double_t GetisoPhotons() const {return _isoPhotons;};
  inline Double_t isoPhotons() const {return _isoPhotons;};

  inline void SetisoPU(Double_t isoPU){_isoPU = isoPU;};
  inline Double_t GetisoPU() const {return _isoPU;};
  inline Double_t isoPU() const {return _isoPU;};

  inline void SetCharge(Double_t charge){_charge = charge;};
  inline Double_t GetCharge() const {return _charge;};
  inline Double_t charge() const {return _charge;};

  inline void SetSCeta(Double_t scEta){_scEta = scEta;};
  inline Double_t GetSCEta() const {return _scEta;};
  inline Double_t scEta() const {return _scEta;};

  inline void SetInCrack(Int_t inCrack){_inCrack =  inCrack;};
  inline Int_t GetInCrack() const {return _inCrack;};
  inline Int_t inCrack() const {return _inCrack;};

  inline void SetSigmaEtaEta(Double_t sigmaEtaEta){_sigmaEtaEta = sigmaEtaEta;};
  inline Double_t GetSigmaEtaEta() const {return _sigmaEtaEta;};
  inline Double_t sigmaEtaEta() const {return _sigmaEtaEta;};

  inline void SetDEtaInSeed(Double_t dEtaInSeed){_dEtaInSeed = dEtaInSeed;};
  inline Double_t GetDEtaInSeed() const {return _dEtaInSeed;};
  inline Double_t dEtaInSeed() const {return _dEtaInSeed;};

  inline void SetDPhiIn(Double_t dPhiIn){_dPhiIn = dPhiIn;};
  inline Double_t GetDPhiIn() const {return _dPhiIn;};
  inline Double_t dPhiIn() const {return _dPhiIn;};

  inline void SetHOverE(Double_t hOverE){_hOverE = hOverE;};
  inline Double_t GetHOverE() const {return _hOverE;};
  inline Double_t hOverE() const {return _hOverE;};

  inline void SetOoEmooP(Double_t ooEmooP){_ooEmooP = ooEmooP;};
  inline Double_t GetOoEmooP() const {return _ooEmooP;};
  inline Double_t ooEmooP() const {return _ooEmooP;};

  inline void SetRelIsoPFRhoEA(Double_t relIsoPFRhoEA){_relIsoPFRhoEA = relIsoPFRhoEA;};
  inline Double_t GetRelIsoPFRhoEA() const {return _relIsoPFRhoEA;};
  inline Double_t relIsoPFRhoEA() const {return _relIsoPFRhoEA;};

  inline void SetMissingHits(Int_t missingHits){_missingHits = missingHits;};
  inline Int_t GetMissingHits() const {return _missingHits;};
  inline Int_t missingHits() const {return _missingHits;};

//inline void SetRelPt(Double_t relpt){_relPt = relPt;};
  //inline Double_t GetRelPt() const {return _relPt;};
   // inline Double_t relPt() const {return _relPt;};
 private:

  Int_t _passVetoId;
  Int_t _passLooseId;
  Int_t _passMediumId;
  Int_t _passTightId;
  Int_t _passHEEPId;
  Int_t _passConversionVeto;

  Int_t _expectedMissingInnerHits;

  Double_t _patElectron_d0;
  Double_t _patElectron_dz;
  Double_t _isoChargedHadrons;
  Double_t _isoNeutralHadrons;
  Double_t _isoPhotons;
  Double_t _isoPU;
  Double_t _charge;
  Double_t _scEta;
  Int_t _inCrack;

  //The ID values
  Double_t _sigmaEtaEta;
  Double_t _dEtaInSeed;
  Double_t _dPhiIn;
  Double_t _hOverE;
  Double_t _ooEmooP;
  Int_t _missingHits;
  
  //Isolation
  Double_t _relIsoPFRhoEA;
  Double_t _minLeptonJetDetaR;
  Double_t _relPtcut;
  Double_t  _minMuonEleDetaR;



  ///////////////////////////////////////////
  // Maps containing the cut values to be placed on the different type of selected leptons
  std::map<TString,Double_t> _maxEtaCuts;
  std::map<TString,Double_t> _minPtCuts;
  std::map<TString,Double_t> _minEtaGapCuts;
  std::map<TString,Double_t> _maxEtaGapCuts;
  std::map<TString,Double_t> _dZCutEndcap;
  std::map<TString,Double_t> _dZCutBarrel;
  std::map<TString,Double_t> _d0CutEndcap;
  std::map<TString,Double_t> _d0CutBarrel;
  Int_t _bstar;  
  Int_t _Elechannel;  
  Int_t _QCD_CR;  
  Int_t _TT_CR;  

  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(Electron,0)
  ////////////////////////////////////////////////////////////////////////////////
};


#endif
