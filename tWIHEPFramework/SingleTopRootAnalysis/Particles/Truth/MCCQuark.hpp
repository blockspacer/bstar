/******************************************************************************
 * MCCQuark.hpp                                                               *
 *                                                                            *
 * Container for a MCCQuark                                                   *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCCQuark class                                  *
 *    MCCQuark()                         -- Default Constructor               *
 *    MCCQuark()                         -- Parameterized Constructor         *
 *    MCCQuark()                         -- Constructor with Particle         *
 *    ~MCCQuark()                        -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill muon from event tree         *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCCQuark         *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCCQuark class                                     *
 *    Double_t _Matched                   -- If BQuark is Matched             *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      30 Jun 2009 - Created by P. Ryan                                      *
 *      22 Jul 2009 - Created by J. Holzbauer, based on MCBQuark              *
 *****************************************************************************/
#ifndef MCCQuark_h
#define MCCQuark_h

#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

class MCCQuark: public MCParticle
{
public:
  // Default Constructor
  MCCQuark();

  // Copy constructor from MCCQuark
  MCCQuark(const MCCQuark& other);

  // Copy constructor from MCParticle
  MCCQuark(const MCParticle& other);

  // Destructor
  ~MCCQuark();

  // Set all contents to their defaults
  inline void Clear() { MCParticle::Clear(); _DeltaR = -1; _Matched = kFALSE;};

  // Fill the MCCQuark from a TruthTree
  void FillCQuark(TruthTree *trtr, Int_t iE);

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
  MCCQuark& operator+=(const MCCQuark& other);
  // +
  MCCQuark operator+(const MCCQuark& other);
  // = MCParticle const
  MCCQuark& operator=(const MCParticle& other);
  // = MCCQuark const
  MCCQuark& operator=(const MCCQuark& other);
  // = MCCQuark non-const
  MCCQuark& operator=(MCCQuark& other);
    
private:
    Double_t _DeltaR;       // Distance in eta-phi space between reconstructed and truth particles
  Bool_t _Matched;        // If particle is matched at reconstructed level


  ///////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCCQuark,0)
  ///////////////////////////////////////////////////////////////////////////////

};


#endif
