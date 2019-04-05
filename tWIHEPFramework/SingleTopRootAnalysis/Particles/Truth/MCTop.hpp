/******************************************************************************
 * MCTop.hpp                                                                  *
 *                                                                            *
 * Container for a MCTop                                                      *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCTop class                                     *
 *    MCTop()                           -- Default Constructor                *
 *    MCTop()                           -- Parameterized Constructor          *
 *    MCTop()                           -- Constructor with Particle          *
 *    ~MCTop()                          -- Destructor                         *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill top from event tree          *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCTop            *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCTop class                                        *
 *                                                                            *
 * History                                                                    *
 *      15 May 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/
#ifndef MCTop_h
#define MCTop_h

#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include "SingleTopRootAnalysis/Trees/TruthTree.hpp"

class MCTop: public MCParticle
{
public:
  // Default Constructor
  MCTop();

  // Copy constructor from MCTop
  MCTop(const MCTop& other);

  // Copy constructor from MCParticle
  MCTop(const MCParticle& other);

  // Destructor
  ~MCTop();

  // Set all contents to their defaults
  inline void Clear() { MCParticle::Clear();};

  // Fill the MCTop from a TruthTree
  void FillTop(TruthTree *trtr, Int_t iE);

  // Getters and Setters


  // Overloaded operators:
  // +=
  MCTop& operator+=(const MCTop& other);
  // +
  MCTop operator+(const MCTop& other);
  // = MCParticle const
  MCTop& operator=(const MCParticle& other);
  // = MCTop const
  MCTop& operator=(const MCTop& other);
  // = MCTop non-const
  MCTop& operator=(MCTop& other);
    
private:

  ///////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(MCTop,0)
  ///////////////////////////////////////////////////////////////////////////////

};


#endif
