/******************************************************************************
 * Jet.hpp                                                                    *
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
 *                                                                            *
 * History                                                                    *
 *      17 June 2015 - Created by Huaqiao ZHANG                               *
 *****************************************************************************/

#ifndef boostedjet_h
#define boostedjet_h

//#include "SingleTopRootAnalysis/Particles/Truth/MCJet.hpp"
#include "SingleTopRootAnalysis/Particles/Truth/MCTau.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Particle.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Electron.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Muon.hpp"
#include "SingleTopRootAnalysis/Trees/EventTree.hpp"
#include "SingleTopRootAnalysis/Trees/FastSimTree.hpp"
//#include "SingleTopRootAnalysis/Base/Dictionary/MultijetJESUncertaintyProvider.hpp"
//#include "SingleTopRootAnalysis/Base/Dictionary/JetEfficiencyEstimator.hpp"
//#include "SingleTopRootAnalysis/Base/Dictionary/JERProvider.hpp"
//
#include <TString.h>

//using namespace Analysis;

class Electron;
class Muon;

class BoostedJet: public Particle
{
 public:
  // Default Constructor
  BoostedJet();

  // Copy Constructor
  BoostedJet(const BoostedJet& other);

  // Construtor using particle
  BoostedJet(const Particle& other);

  // Destructor
  ~BoostedJet();

  // Set all contents to their defaults
  inline void Clear() { Particle::Clear();_prunedmass=0.0;_softdropmass=0.0;_tau1=0.0;_tau2=0.0;_tau3=0.0; _numberOfConstituents=0; _chargedMultiplicity=0;  _bDiscriminator = -999.0; _mass = 0.0; _uncorrPt = 0.0; _neutralHadEnergyFraction=0.0; _neutralEmEmEnergyFraction = 0.0; _chargedHadronEnergyFraction =0.0; _chargedEmEnergyFraction=0.0; _muonEnergyFraction=0.0; _electronEnergy=0.0; _photonEnergy=0.0; _tagged = kFALSE; _hadronFlavour=-1;
}

  void SetCuts(TEnv* config);

  // Fill the boostedjet from an EventTree 
  Bool_t Fill( double myJESCorr, double myJERCorr, std::vector<Muon>& selectedMuons, std::vector<Electron>& selectedElectrons, EventTree *evtr, Int_t iE, TLorentzVector * met, bool isMC);
  //  Bool_t Fill( double myJESCorr, double myJERCorr, std::vector<Electron>& selectedElectrons, EventTree *evtr, Int_t iE);

  inline void Setprunedmass(Double_t prunedmass){_prunedmass =prunedmass;};
  inline Double_t Getprunedmass() const {return _prunedmass;};
  inline Double_t prunedmass() const {return _prunedmass;};

  inline void Setsoftdropmass(Double_t softdropmass){_softdropmass =softdropmass;};
  inline Double_t Getsoftdropmass() const {return _softdropmass;};
  inline Double_t softdropmass() const {return _softdropmass;};


  inline void Settau1(Double_t tau1){_tau1 = tau1;};
  inline Double_t Gettau1() const {return _tau1;};
  inline Double_t tau1() const {return _tau1;};

  inline void Settau2(Double_t tau2){_tau2 = tau2;};
  inline Double_t Gettau2() const {return _tau2;};
  inline Double_t tau2() const {return _tau2;};

  inline void Settau3(Double_t tau3){_tau3 = tau3;};
  inline Double_t Gettau3() const {return _tau3;};
  inline Double_t tau3() const {return _tau3;};


  inline void SetnumberOfConstituents(Int_t numberOfConstituents){_numberOfConstituents = numberOfConstituents;};
  inline Int_t GetnumberOfConstituents() const {return _numberOfConstituents;};
  inline Int_t numberOfConstituents() const {return _numberOfConstituents;};

  inline void SethadronFlavour(Int_t hadronFlavour){_hadronFlavour = hadronFlavour;};
  inline Int_t GethadronFlavour() const {return _hadronFlavour;};
  inline Int_t hadronFlavour() const {return _hadronFlavour;};

  inline void SetchargedMultiplicity(Int_t chargedMultiplicity){_chargedMultiplicity = chargedMultiplicity;};
  inline Int_t GetchargedMultiplicity() const {return _chargedMultiplicity;};
  inline Int_t chargedMultiplicity() const {return _chargedMultiplicity;};

  inline void SetbDiscriminator(Double_t bDiscriminator){_bDiscriminator = bDiscriminator;};
  inline Double_t GetbDiscriminator() const {return _bDiscriminator;};
  inline Double_t bDiscriminator() const {return _bDiscriminator;};


  inline void Setmass(Double_t mass){_mass = mass;};
  inline Double_t Getmass() const {return _mass;};
  inline Double_t mass() const {return _mass;};

  inline void SetneutralHadEnergyFraction(Double_t neutralHadEnergyFraction){_neutralHadEnergyFraction = neutralHadEnergyFraction;};
  inline Double_t GetneutralHadEnergyFraction() const {return _neutralHadEnergyFraction;};
  inline Double_t neutralHadEnergyFraction() const {return _neutralHadEnergyFraction;};

  inline void SetneutralEmEmEnergyFraction(Double_t neutralEmEmEnergyFraction){_neutralEmEmEnergyFraction = neutralEmEmEnergyFraction;};
  inline Double_t GetneutralEmEmEnergyFraction() const {return _neutralEmEmEnergyFraction;};
  inline Double_t neutralEmEmEnergyFraction() const {return _neutralEmEmEnergyFraction;};

  inline void SetchargedHadronEnergyFraction(Double_t chargedHadronEnergyFraction){_chargedHadronEnergyFraction = chargedHadronEnergyFraction;};
  inline Double_t GetchargedHadronEnergyFraction() const {return _chargedHadronEnergyFraction;};
  inline Double_t chargedHadronEnergyFraction() const {return _chargedHadronEnergyFraction;};

  inline void SetchargedEmEnergyFraction(Double_t chargedEmEnergyFraction){_chargedEmEnergyFraction = chargedEmEnergyFraction;};
  inline Double_t GetchargedEmEnergyFraction() const {return _chargedEmEnergyFraction;};
  inline Double_t chargedEmEnergyFraction() const {return _chargedEmEnergyFraction;};

  inline void SetmuonEnergyFraction(Double_t muonEnergyFraction){_muonEnergyFraction = muonEnergyFraction;};
  inline Double_t GetmuonEnergyFraction() const {return _muonEnergyFraction;};
  inline Double_t muonEnergyFraction() const {return _muonEnergyFraction;};

  inline void SetelectronEnergy(Double_t electronEnergy){_electronEnergy = electronEnergy;};
  inline Double_t GetelectronEnergy() const {return _electronEnergy;};
  inline Double_t electronEnergy() const {return _electronEnergy;};

  inline void SetphotonEnergy(Double_t photonEnergy){_photonEnergy = photonEnergy;};
  inline Double_t GetphotonEnergy() const {return _photonEnergy;};
  inline Double_t photonEnergy() const {return _photonEnergy;};

  inline void SetuncorrPt(Double_t uncorrPt){_uncorrPt = uncorrPt;};
  inline Double_t GetuncorrPt() const {return _uncorrPt;};
  inline Double_t uncorrPt() const {return _uncorrPt;};

  inline void SetTagged(Int_t isTagged){_tagged = isTagged;};
  inline Int_t IsTagged() const {return _tagged;};
  inline Int_t tagged() const {return _tagged;};

  inline void SetClosestLep(Double_t closestLep){_closestLep = closestLep;};
  inline Double_t GetClosestLep() const {return _closestLep;};
  inline Double_t closestLep() const {return _closestLep;};

  // Overloaded Operators
  // +=
  BoostedJet& operator+=(const BoostedJet& other);
  // +
  BoostedJet operator+(const BoostedJet& other);
  // = BoostedJet const
  BoostedJet& operator=(const Particle& other);
  // = Particle const
  BoostedJet& operator=(const BoostedJet& other);
  // = BoostedJet non-const
  BoostedJet& operator=(BoostedJet& other);
    

 private:
  Int_t _numberOfConstituents;   
  Int_t _chargedMultiplicity;   
  Int_t _hadronFlavour;
  Double_t _bDiscriminator;   
  Double_t _mass;   
  Double_t _prunedmass;   
  Double_t _softdropmass;   
  Double_t _tau1;   
  Double_t _tau2;   
  Double_t _tau3;   
  Double_t _neutralHadEnergyFraction;   
  Double_t _neutralEmEmEnergyFraction;   
  Double_t _chargedHadronEnergyFraction;   
  Double_t _chargedEmEnergyFraction;   
  Double_t _muonEnergyFraction;   
  Double_t _electronEnergy;   
  Double_t _photonEnergy;   
  Double_t _uncorrPt;  
  Int_t _tagged;
  Double_t _closestLep;

  // Cuts applied to the boostedjet objects
  Double_t _maxEtaCut;
  Double_t _minPtCut;
  Double_t _minCRPtCut;
  Double_t _bMaxEtaCut;
  Double_t _bMinPtCut;
  Double_t _bTagCut;
  Double_t _closestLeptonCut;
  Double_t _minMuonBoostedJetDetaR;
  Double_t _minprunedmass;
  Double_t _maxprunedmass;
  Double_t _minQCDprunedmass;
  Double_t _maxQCDprunedmass;
  Double_t _minsoftdropmass;
  Double_t _maxsoftdropmass;
  Double_t _tau21;
  Double_t _tau32;

  Int_t    _bstar;
  Int_t    _Elechannel;
  Int_t    _TT_CR;
  Int_t    _QCD_CR;
  // Are we running systematic variations?
  Int_t _jesUp;
  Int_t _jesDown;
  Int_t _jerUp;
  Int_t _jerDown;

  // Apply the boostedjet correction systematics
  void SystematicPtShift(EventTree * evtr, Int_t iE, TLorentzVector * met);
 
  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(BoostedJet,0)
  ////////////////////////////////////////////////////////////////////////////////
};


#endif
