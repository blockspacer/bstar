/******************************************************************************
 * MCW.hpp                                                                    *
 *                                                                            *
 * Container for a MCW                                                        *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCW class                                       *
 *    MCW()                           -- Default Constructor                  *
 *    MCW()                           -- Parameterized Constructor            *
 *    MCW()                           -- Constructor with Particle            *
 *    ~MCW()                          -- Destructor                           *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill w from event tree            *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCW              *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCW class                                          *
 *                                                                            *
 * History                                                                    *
 *      15 May 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/
#ifndef MCW_h
#define MCW_h

#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

class MCW: public MCParticle
{
public:
  // Default Constructor
  MCW();

  // Copy constructor from MCW
  MCW(const MCW& other);

  // Copy constructor from MCParticle
  MCW(const MCParticle& other);

  // Destructor
  ~MCW();

  // Set all contents to their defaults
  inline void Clear() { MCParticle::Clear();};

  // Fill the MCW from a TruthTree
  void FillW(TruthTree *trtr, Int_t iE);

  // Getters and Setters


  // Overloaded operators:
  // +=
  MCW& operator+=(const MCW& other);
  // +
  MCW operator+(const MCW& other);
  // = MCParticle const
  MCW& operator=(const MCParticle& other);
  // = MCW const
  MCW& operator=(const MCW& other);
  // = MCW non-const
  MCW& operator=(MCW& other);
    
private:

  ///////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCW,0)
  ///////////////////////////////////////////////////////////////////////////////

};


#endif
