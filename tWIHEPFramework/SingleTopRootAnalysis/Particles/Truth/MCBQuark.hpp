/******************************************************************************
 * MCBQuark.hpp                                                               *
 *                                                                            *
 * Container for a MCBQuark                                                   *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCBQuark class                                  *
 *    MCBQuark()                         -- Default Constructor               *
 *    MCBQuark()                         -- Parameterized Constructor         *
 *    MCBQuark()                         -- Constructor with Particle         *
 *    ~MCBQuark()                        -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill muon from event tree         *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCBQuark         *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCBQuark class                                     *
 *    Double_t _Matched                   -- If BQuark is Matched             *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      30 Jun 2009 - Created by P. Ryan                                      *
 *                                                                            *
 *****************************************************************************/
#ifndef MCBQuark_h
#define MCBQuark_h

#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

class MCBQuark: public MCParticle
{
public:
  // Default Constructor
  MCBQuark();

  // Copy constructor from MCBQuark
  MCBQuark(const MCBQuark& other);

  // Copy constructor from MCParticle
  MCBQuark(const MCParticle& other);

  // Destructor
  ~MCBQuark();

  // Set all contents to their defaults
  inline void Clear() { MCParticle::Clear(); _DeltaR = -1; _Matched = kFALSE;};

  // Fill the MCBQuark from a TruthTree
  void FillBQuark(TruthTree *trtr, Int_t iE);

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
  MCBQuark& operator+=(const MCBQuark& other);
  // +
  MCBQuark operator+(const MCBQuark& other);
  // = MCParticle const
  MCBQuark& operator=(const MCParticle& other);
  // = MCBQuark const
  MCBQuark& operator=(const MCBQuark& other);
  // = MCBQuark non-const
  MCBQuark& operator=(MCBQuark& other);
    
private:
    Double_t _DeltaR;       // Distance in eta-phi space between reconstructed and truth particles
  Bool_t _Matched;        // If particle is matched at reconstructed level


  ///////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCBQuark,0)
  ///////////////////////////////////////////////////////////////////////////////

};


#endif
