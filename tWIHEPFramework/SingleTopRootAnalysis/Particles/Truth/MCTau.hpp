/******************************************************************************
 * MCTau.hpp                                                                  *
 *                                                                            *
 * Container for a MCTau                                                      *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCTau class                                     *
 *    MCTau()                           -- Default Constructor                *
 *    MCTau()                           -- Parameterized Constructor          *
 *    MCTau()                           -- Constructor with Particle          *
 *    ~MCTau()                          -- Destructor                         *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill muon from event tree         *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCTau            *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCTau class                                        *
 *    Double_t _Matched                   -- If Tau is Matched                *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      13 Apr 2006 - Created by P. Ryan / J. Holzbauer                       *
 *      25 Apr 2007 - Modified by J. Holzbauer, added DeltaR from tree        *
 *                    set default deltaR to -1                                *
 *****************************************************************************/
#ifndef MCTau_h
#define MCTau_h

#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

class MCTau: public MCParticle
{
public:
  // Default Constructor
  MCTau();

  // Copy constructor from MCTau
  MCTau(const MCTau& other);

  // Copy constructor from MCParticle
  MCTau(const MCParticle& other);

  // Destructor
  ~MCTau();

  // Set all contents to their defaults
  inline void Clear() { MCParticle::Clear(); _DeltaR = -1; _Matched = kFALSE;};

  // Fill the MCTau from a TruthTree
  void FillTau(TruthTree *trtr, Int_t iE);

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
  MCTau& operator+=(const MCTau& other);
  // +
  MCTau operator+(const MCTau& other);
  // = MCParticle const
  MCTau& operator=(const MCParticle& other);
  // = MCTau const
  MCTau& operator=(const MCTau& other);
  // = MCTau non-const
  MCTau& operator=(MCTau& other);
    
private:
    Double_t _DeltaR;       // Distance in eta-phi space between reconstructed and truth particles
  Bool_t _Matched;        // If particle is matched at reconstructed level


  ///////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCTau,0)
  ///////////////////////////////////////////////////////////////////////////////

};


#endif
