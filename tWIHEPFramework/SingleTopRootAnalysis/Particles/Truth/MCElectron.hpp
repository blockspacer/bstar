/******************************************************************************
 * MCElectron.hpp                                                             *
 *                                                                            *
 * Container for a MCElectron                                                 *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCElectron class                                *
 *    MCElectron()                           -- Default Constructor           *
 *    MCElectron()                           -- Parameterized Constructor     *
 *    MCElectron()                           -- Constructor with Particle     *
 *    ~MCElectron()                          -- Destructor                    *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill muon from event tree         *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCElectron       *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCElectron class                                   *
 *    Double_t _Matched                   -- If Electron is Matched           *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      13 Apr 2006 - Created by P. Ryan / J. Holzbauer                       *
 *      25 Apr 2007 - Modified by J. Holzbauer, added DeltaR from tree        *
 *                    set default deltaR to -1                                *
 *****************************************************************************/
#ifndef MCElectron_h
#define MCElectron_h

#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

class MCElectron: public MCParticle
{
public:
  // Default Constructor
  MCElectron();

  // Copy constructor from MCElectron
  MCElectron(const MCElectron& other);

  // Copy constructor from MCParticle
  MCElectron(const MCParticle& other);

  // Destructor
  ~MCElectron();

  // Set all contents to their defaults
  inline void Clear() { MCParticle::Clear(); _DeltaR = -1; _Matched = kFALSE;};

  // Fill the MCElectron from a TruthTree
  void FillElectron(TruthTree *trtr, Int_t iE);

  // Getters and Setters
  // DeltaR
  inline void SetDeltaR(const Double_t& deltaR) { _DeltaR = deltaR; };
  inline Double_t GetDeltaR() const { return _DeltaR; };
  //Is no DeltaR() call because want to enable TLorentz DeltaR as the default method

  // Matched
  inline void SetMatched(const Bool_t Matched) { _Matched = Matched; };
  inline Bool_t GetMatched() const { return _Matched; };
  inline Bool_t Matched() const { return GetMatched(); };

  // Overloaded operators:
  // +=
  MCElectron& operator+=(const MCElectron& other);
  // +
  MCElectron operator+(const MCElectron& other);
  // = MCParticle const
  MCElectron& operator=(const MCParticle& other);
  // = MCElectron const
  MCElectron& operator=(const MCElectron& other);
  // = MCElectron non-const
  MCElectron& operator=(MCElectron& other);
    
private:
    Double_t _DeltaR;       // Distance in eta-phi space between reconstructed and truth particles
  Bool_t _Matched;        // If particle is matched at reconstructed level


  ///////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCElectron,0)
  ///////////////////////////////////////////////////////////////////////////////

};


#endif
