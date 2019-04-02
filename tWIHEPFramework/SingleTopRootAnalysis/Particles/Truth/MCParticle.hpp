/******************************************************************************
 * MCParticle.hpp                                                             *
 *                                                                            *
 * Container for a MCParticle                                                 *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCParticle class                                *
 *    MCParticle()                       -- Default Constructor               *
 *    MCParticle()                       -- Parameterized Constructor         *
 *    MCParticle()                       -- Constructor with Particle         *
 *    ~MCParticle()                      -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill                               -- Fill electron from event tree     *
 *    SetPdgId                           -- Set Particle Data Group ID        *
 *    SetStatus                          -- Set Status                        *
 *    SetBarCode                         -- Set Bar Code                      *
 *    GetPdgId                           -- Get Particle Data Group ID        *
 *    PdgId                              -- Get Particle Data Group ID        *
 *    GetAbsPdgId                        -- Get abs Particle Data Group ID    *
 *    AbsPdgId                           -- Get abs Particle Data Group ID    *
 *    GetStatus                          -- Get Status                        *
 *    Status                             -- Get Status                        *
 *    GetBarCode                         -- Get Bar Code                      *
 *    BarCode                            -- Get Bar Code                      *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCParticle       *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCParticle class                                   *
 *    Int_t _PdgId                       -- PDG ID for this particle          *
 *    Int_t _Status                      -- Status flag for this particle     *
 *    Int_t _BarCode                     -- BarCode flag for this particle    *
 *                                                                            *
 * History                                                                    *
 *      15 Nov 2006 - Created by R. Schwienhorst                              *
 *      20 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      23 Apr 2007 - Modified by J. Holzbauer - Added more PDGID requirements*
 *      08 May 2007 - Modified by J. Holzbauer- added FillDeltaR              *
 *      17 Jul 2007 - Modified by J. Holzbauer- added more jet class tags     *
 *      19 Jul 2007 - RS: Make all "is..." methods return a bool.             *
 *****************************************************************************/

#ifndef MCParticle_h
#define MCParticle_h

#include "SingleTopRootAnalysis/Particles/Recon/Particle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"
//

class MCParticle: public Particle
{
public:
  // Default constructor
  MCParticle();

  // Copy constructor
  MCParticle(const MCParticle& other);

  // Copy Constructor - Particle
  MCParticle(const Particle& other);

  // Destructor
  ~MCParticle();

  // Set all contents to their defaults
  inline void Clear() { Particle::Clear(); _PdgId = 0; _Status = 0; _BarCode = 0;  _nDecay = 0;};

  // Fill the MCParticle from a TruthTree
  void Fill(TruthTree *trtr, Int_t iE);

 
  // Setters PDG ID
  inline void SetPdgId(Double_t pdgId)  { _PdgId   = pdgId;   };
  inline void SetStatus(Int_t status)   { _Status  = status;  };
  inline void SetnDecay(Int_t nparents)   { _nDecay  = nparents;  };
  inline void SetBarCode(Int_t barCode) { _BarCode = barCode; };


  // Getters PDG ID
  inline Double_t GetPdgId() const { return _PdgId; };
  inline Double_t PdgId() const { return GetPdgId(); };
  inline Double_t GetAbsPdgId() const { return(_PdgId>0 ? _PdgId : -_PdgId) ; };
  inline Double_t AbsPdgId() const { return GetAbsPdgId(); };
  


  // Particles ID check functions
  inline bool isMuon() const {return(AbsPdgId()==13); };
  inline bool isElectron() const {return(AbsPdgId()==11); };
  inline bool isTau() const {return(AbsPdgId()==15); };
  inline bool isTop() const {return(AbsPdgId()==6); };
  inline bool isW() const {return(AbsPdgId()==24); };
  inline bool isChargedPion() const {return(AbsPdgId()==211); };
  inline bool isChargedK() const {return(AbsPdgId()==321); };
  inline bool isBJet() const {return(AbsPdgId()==5); };
  inline bool isCJet() const {return(AbsPdgId()==4); };
  inline bool isLightQuarkJet() const {return(AbsPdgId()>=1 && AbsPdgId()<=3); };
  inline bool isJet() const {return(AbsPdgId()>=1 && AbsPdgId()<=5); };
  inline bool isQuark() const {return(AbsPdgId()>=1 && AbsPdgId()<=6); };
  inline bool isNeutrino() const {return(AbsPdgId()==12 || AbsPdgId()==14|| AbsPdgId()==16); };

  // Status flag
  inline Int_t GetStatus() const { return _Status; };
  inline Int_t Status() const { return GetStatus(); };
 
  // nDecay flag
  inline Int_t GetnDecay() const { return _nDecay; };
  inline Int_t nDecay() const { return GetnDecay(); };

  // Stable flag
  inline Int_t GetBarCode() const { return _BarCode; };
  inline Int_t BarCode() const { return GetBarCode(); };


  // Overloaded operators:
  // +=
  MCParticle& operator+=(const MCParticle& other);
  // +
  MCParticle operator+(const MCParticle& other);
  // = const Particle
  MCParticle& operator=(const Particle& other);
  // = MCParticle const
  MCParticle& operator=(const MCParticle& other);
  // = MC Particle const
  MCParticle& operator=(MCParticle& other);

protected:
  Double_t _PdgId;   // PDG ID for this particle
  Int_t _Status;     // Status flag for this particle
  Int_t _nDecay;     // Status flag for this particle
  Int_t _BarCode;    // BarCode flag for this particle
 
  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCParticle,0)
  ////////////////////////////////////////////////////////////////////////////////

};

#endif
