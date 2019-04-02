/******************************************************************************
 * MCNeutrino.hpp                                                                  *
 *                                                                            *
 * Container for a MCNeutrino                                                      *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCNeutrino class                                     *
 *    MCNeutrino()                           -- Default Constructor                *
 *    MCNeutrino()                           -- Parameterized Constructor          *
 *    MCNeutrino()                           -- Constructor with Particle          *
 *    ~MCNeutrino()                          -- Destructor                         *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill top from event tree          *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCNeutrino            *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCNeutrino class                                        *
 *                                                                            *
 * History                                                                    *
 *      15 May 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/
#ifndef MCNeutrino_h
#define MCNeutrino_h

#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

class MCNeutrino: public MCParticle
{
public:
  // Default Constructor
  MCNeutrino();

  // Copy constructor from MCNeutrino
  MCNeutrino(const MCNeutrino& other);

  // Copy constructor from MCParticle
  MCNeutrino(const MCParticle& other);

  // Destructor
  ~MCNeutrino();

  // Set all contents to their defaults
  inline void Clear() { MCParticle::Clear();};

  // Fill the MCNeutrino from a TruthTree
  void FillNeutrino(TruthTree *trtr, Int_t iE);

  // Getters and Setters


  // Overloaded operators:
  // +=
  MCNeutrino& operator+=(const MCNeutrino& other);
  // +
  MCNeutrino operator+(const MCNeutrino& other);
  // = MCParticle const
  MCNeutrino& operator=(const MCParticle& other);
  // = MCNeutrino const
  MCNeutrino& operator=(const MCNeutrino& other);
  // = MCNeutrino non-const
  MCNeutrino& operator=(MCNeutrino& other);
    
private:

  ///////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCNeutrino,0)
  ///////////////////////////////////////////////////////////////////////////////

};


#endif
