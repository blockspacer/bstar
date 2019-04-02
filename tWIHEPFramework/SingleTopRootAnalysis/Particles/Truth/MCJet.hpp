/******************************************************************************
 * MCJet.hpp                                                                  *
 *                                                                            *
 * Container for a MCJet                                                      *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCJet class                                     *
 *    MCJet()                           -- Default Constructor                *
 *    MCJet()                           -- Parameterized Constructor          *
 *    MCJet()                           -- Constructor with Particle          *
 *    ~MCJet()                          -- Destructor                         *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill muon from event tree         *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCJet            *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCJet class                                        *
 *    Double_t _Matched                   -- If Jet is Matched                *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      13 Apr 2007 - Created by P. Ryan / J. Holzbauer                       *
 *      25 Apr 2007 - Modified by J. Holzbauer, added DeltaR from tree        *
 *                    set default deltaR to -1                                *
 *****************************************************************************/
#ifndef MCJet_h
#define MCJet_h

#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

class MCJet: public MCParticle
{
public:
  // Default Constructor
  MCJet();

  // Copy constructor from MCJet
  MCJet(const MCJet& other);

  // Copy constructor from MCParticle
  MCJet(const MCParticle& other);

  // Destructor
  ~MCJet();

  // Set all contents to their defaults
  inline void Clear() { MCParticle::Clear(); _DeltaR = -1; _Matched = kFALSE;};

  // Fill the MCJet from a TruthTree
  void FillJet(TruthTree *trtr, Int_t iE);

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
  MCJet& operator+=(const MCJet& other);
  // +
  MCJet operator+(const MCJet& other);
  // = MCParticle const
  MCJet& operator=(const MCParticle& other);
  // = MCJet const
  MCJet& operator=(const MCJet& other);
  // = MCJet non-const
  MCJet& operator=(MCJet& other);
    
private:
    Double_t _DeltaR;       // Distance in eta-phi space between reconstructed and truth particles
  Bool_t _Matched;        // If particle is matched at reconstructed level


  ///////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCJet,0)
  ///////////////////////////////////////////////////////////////////////////////

};


#endif
