/******************************************************************************
 * Neutrino.hpp                                                               *
 *                                                                            *
 * Store information about final state neutrinos                              *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of Neutrino class                                  *
 *    Neutrino()                         -- Default Constructor               *
 *    Neutrino()                         -- Parameterized Constructor         *
 *    Clear()                            -- Set contents to default           *
 *    Fill                               -- Fill neutrino from event tree     *
 *    FillFastSim                        -- Fill neutrino from FastSim tree   *
 *    SetEoverP                          -- Set E over P for neutrino         *
 *    GetEoverP                          -- GetEoverP                         *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const Neutrino         *
 *    operator =                         -- Overload =                        *
 *                                                                            *
 * Private Data Members of Neutrino class                                     *
 *    Double_t _EoverP;                  -- E over P                          *
 *                                                                            *
 * History                                                                    *
 *      17 June 2015 - Port from STAR by Huaqiao ZHANG                        *
 *****************************************************************************/

#ifndef neutrino_h
#define neutrino_h

#include "SingleTopRootAnalysis/Particles/Recon/Particle.hpp"
#include "SingleTopRootAnalysis/Trees/EventTree.hpp"
#include "SingleTopRootAnalysis/Trees/FastSimTree.hpp"

//#include "SingleTopRootAnalysis/utils.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Jet.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Electron.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Muon.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Tau.hpp"


class Neutrino: public Particle
{
 public:

  // Default constructor
  Neutrino();

  // Copy constructor for Neutrinos
  Neutrino(const Neutrino& other);

  // Copy constructor for Particles
  Neutrino(const Particle& other);

  // Destructor
  ~Neutrino();

  // Set all contents to their defaults
  inline void Clear() { Particle::Clear();  }

  // Fill the neutrino from an EventTree
  Bool_t Fill(const std::vector<Muon>& tightMuons,const std::vector<Electron>& tightElectrons,const Double_t& missingEx, const Double_t& missingEy);


  // Overloaded operators:
  // +=
  Neutrino& operator+=(const Neutrino& other);
  // +
  Neutrino operator+(const Neutrino& other);
  // = Particle
  Neutrino& operator=(const Particle& other);
  // = Neutrino const
  Neutrino& operator=(const Neutrino& other);
  // = Neutrino non-const
  Neutrino& operator=(Neutrino& other);

 private:

 
  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(Neutrino,0)
  ////////////////////////////////////////////////////////////////////////////////
};


#endif
