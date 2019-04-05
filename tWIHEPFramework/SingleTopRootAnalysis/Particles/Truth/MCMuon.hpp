/******************************************************************************
 * MCMuon.hpp                                                                 *
 *                                                                            *
 * Container for a MCMuon                                                     *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCMuon class                                    *
 *    MCMuon()                           -- Default Constructor               *
 *    MCMuon()                           -- Parameterized Constructor         *
 *    MCMuon()                           -- Constructor with Particle         *
 *    ~MCMuon()                          -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill muon from event tree         *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCMuon           *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCMuon class                                       *
 *    Double_t _Matched                   -- If Muon is Matched               *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      13 Apr 2006 - Created by P. Ryan / J. Holzbauer                       *
 *      25 Apr 2007 - Modified by J. Holzbauer, added back DeltaR from tree   *
 *                    set default deltaR to -1                                *
 *****************************************************************************/
#ifndef MCMuon_h
#define MCMuon_h

#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

class MCMuon: public MCParticle
{
public:
  // Default Constructor
  MCMuon();

  // Copy constructor from MCMuon
  MCMuon(const MCMuon& other);

  // Copy constructor from MCParticle
  MCMuon(const MCParticle& other);

  // Destructor
  ~MCMuon();

  // Set all contents to their defaults
  inline void Clear() { MCParticle::Clear(); _DeltaR = -1; _Matched = kFALSE;};

  // Fill the MCMuon from a TruthTree
  void FillMuon(TruthTree *trtr, Int_t iE);

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
  MCMuon& operator+=(const MCMuon& other);
  // +
  MCMuon operator+(const MCMuon& other);
  // = MCParticle const
  MCMuon& operator=(const MCParticle& other);
  // = MCMuon const
  MCMuon& operator=(const MCMuon& other);
  // = MCMuon non-const
  MCMuon& operator=(MCMuon& other);
    
private:
  
  Double_t _DeltaR;       // Distance in eta-phi space between reconstructed and truth particles

  Bool_t _Matched;        // If particle is matched at reconstructed level
  //  Double_t _PgdId;
  //Int_t _Status;
  //Int_t _BarCode;

  ///////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCMuon,0)
  ///////////////////////////////////////////////////////////////////////////////

};


#endif
