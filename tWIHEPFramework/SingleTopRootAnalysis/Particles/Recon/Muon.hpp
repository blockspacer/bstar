/******************************************************************************
 * Muon.hpp                                                                   *
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
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const Muon             *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of Muon class                                         *
 *                                                                            *
 * History                                                                    *
 *      17 June 2015 - Created by Huaqiao ZHANG                               *
 *****************************************************************************/

#ifndef muon_h
#define muon_h

#include "SingleTopRootAnalysis/Particles/Recon/Jet.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Particle.hpp"
#include "SingleTopRootAnalysis/Trees/EventTree.hpp"
#include "SingleTopRootAnalysis/Trees/FastSimTree.hpp"
#include <TEnv.h>

using namespace std;
class Muon: public Particle
{
public:

  // Default Constructor
  Muon();

  // Copy constructor from Muon
  Muon(const Muon& other);

  // Copy constructor from Particle
  Muon(const Particle& other);

  // Dectructor
  ~Muon();

  // Set all contents to their defaults
  inline void Clear() { Particle::Clear(); _passTightId = kFALSE; _passLooseId = kFALSE;_isSoft = kFALSE; _isHighPt = kFALSE; _isoCharged = 0.0; _isoSum = 0.0; _isoCharParPt = 0.0; _isoNeutralHadron = 0.0;  _isoPhoton = 0.0;  _isoPU = 0.0; _chi2=0.0; _dxy=0.0; _validHits = 0.0; _validHitsInner = 0.0; _matchedStat=0.0; _TLayers=0.0; _relIsoR04=0.0;};

  // Fill the muon from an EventTree
  Bool_t Fill(EventTree *evtr,int iE,TString muonType, Bool_t isSimulation);
  // also fill from a fastsim tree
  Bool_t FillFastSim(TEnv *config, FastSimTree *tr, Int_t iE,TString muonType);

 // Overloaded operators:
  // +=
  Muon& operator+=(const Muon& other);
  // +
  Muon operator+(const Muon& other);
  // = Particle
  Muon& operator=(const Particle& other);
  // = Muon const
  Muon& operator=(const Muon& other);
  // = Muon non-const
  Muon& operator=(Muon& other);

  void SetCuts(TEnv *config, TString muonType);
    
  inline void SetpassTightId(Bool_t passTightId){_passTightId = passTightId;};
  inline Bool_t GetpassTightId() const {return _passTightId;};
  inline Bool_t passTightId() const {return _passTightId;};

  inline void SetpassLooseId(Bool_t passLooseId){_passLooseId = passLooseId;};
  inline Bool_t GetpassLooseId() const {return _passLooseId;};
  inline Bool_t passLooseId() const {return _passLooseId;};

  inline void SetisSoft(Bool_t isSoft){_isSoft = isSoft;};
  inline Bool_t GetisSoft() const {return _isSoft;};
  inline Bool_t isSoft() const {return _isSoft;};

  inline void SetisHighPt(Bool_t isHighPt){_isHighPt = isHighPt;};
  inline Bool_t GetisHighPt() const {return _isHighPt;};
  inline Bool_t isHighPt() const {return _isHighPt;};

  inline void SetisPf(Bool_t isPf){_isPf = isPf;};
  inline Bool_t GetisPf() const {return _isPf;};
  inline Bool_t isPf() const {return _isPf;};

  inline void SetisoCharged(Double_t isoCharged){_isoCharged = isoCharged;};
  inline Double_t GetisoCharged() const {return _isoCharged;};
  inline Double_t isoCharged() const {return _isoCharged;};

  inline void SetisoSum(Double_t isoSum){_isoSum = isoSum;};
  inline Double_t GetisoSum() const {return _isoSum;};
  inline Double_t isoSum() const {return _isoSum;};

  inline void SetisoCharParPt(Double_t isoCharParPt){_isoCharParPt = isoCharParPt;};
  inline Double_t GetisoCharParPt() const {return _isoCharParPt;};
  inline Double_t isoCharParPt() const {return _isoCharParPt;};

  inline void SetisoNeutralHadron(Double_t isoNeutralHadron){_isoNeutralHadron = isoNeutralHadron;};
  inline Double_t GetisoNeutralHadron() const {return _isoNeutralHadron;};
  inline Double_t isoNeutralHadron() const {return _isoNeutralHadron;};

  inline void SetisoPhoton(Double_t isoPhoton){_isoPhoton = isoPhoton;};
  inline Double_t GetisoPhoton() const {return _isoPhoton;};
  inline Double_t isoPhoton() const {return _isoPhoton;};

  inline void SetisoPU(Double_t isoPU){_isoPU = isoPU;};
  inline Double_t GetisoPU() const {return _isoPU;};
  inline Double_t isoPU() const {return _isoPU;};

  inline void Setchi2(Double_t chi2){_chi2 = chi2;};
  inline Double_t Getchi2() const {return _chi2;};
  inline Double_t chi2() const {return _chi2;};

  inline void SetvalidHits(Double_t validHits){_validHits = validHits;};
  inline Double_t GetvalidHits() const {return _validHits;};
  inline Double_t validHits() const {return _validHits;};

  inline void SetvalidHitsInner(Double_t validHitsInner){_validHitsInner = validHitsInner;};
  inline Double_t GetvalidHitsInner() const {return _validHitsInner;};
  inline Double_t validHitsInner() const {return _validHitsInner;};

  inline void SetmatchedStat(Double_t matchedStat){_matchedStat = matchedStat;};
  inline Double_t GetmatchedStat() const {return _matchedStat;};
  inline Double_t matchedStat() const {return _matchedStat;};

  inline void Setdxy(Double_t dxy){_dxy = dxy;};
  inline Double_t Getdxy() const {return _dxy;};
  inline Double_t dxy() const {return _dxy;};

  inline void Setdz(Double_t dz){_dz = dz;};
  inline Double_t Getdz() const {return _dz;};
  inline Double_t dz() const {return _dz;};

  inline void SetTLayers(Double_t TLayers){_TLayers = TLayers;};
  inline Double_t GetTLayers() const {return _TLayers;};
  inline Double_t TLayers() const {return _TLayers;};

  inline void SetrelIsoR04(Double_t relIsoR04){_relIsoR04 = relIsoR04;};
  inline Double_t GetrelIsoR04() const {return _relIsoR04;};
  inline Double_t relIsoR04() const {return _relIsoR04;};

  inline void Setndof(Double_t ndof){_ndof = ndof;};
  inline Double_t Getndof() const {return _ndof;};
  inline Double_t ndof() const {return _ndof;};

  inline void SetCharge(Double_t charge){_charge = charge;};
  inline Double_t GetCharge() const {return _charge;};
  inline Double_t charge() const {return _charge;};

private:

  Bool_t _passTightId;
  Bool_t _passLooseId;
  Bool_t _isSoft;
  Bool_t _isHighPt;
  Bool_t _isPf;
  Double_t _isoCharged;
  Double_t _isoSum;
  Double_t _isoCharParPt;
  Double_t _isoNeutralHadron;
  Double_t _isoPhoton;
  Double_t _isoPU;
  Double_t _chi2;
  Double_t _validHits;
  Double_t _validHitsInner;
  Double_t _matchedStat;
  Double_t _dxy;
  Double_t _dz;
  Double_t _TLayers;
  Double_t _relIsoR04;
  Double_t _ndof;
  Double_t _charge;

  //////////////////////////////////
  // Definitions of the objects go here. This way we don't need to access the configuration file for every particle
  map<TString,Double_t> _minPtCuts;
  map<TString,Double_t> _maxEtaCuts;
  map<TString,Double_t> _maxRelIsoCuts;
  
  Int_t _bstar;

  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(Muon,0)
  ////////////////////////////////////////////////////////////////////////////////

};


#endif
