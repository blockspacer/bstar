/******************************************************************************
 * Tau.hpp                                                                    *
 *                                                                            *
 * Store information about final state Taus                                   *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of Tau class                                       *
 *    Tau()                              -- Default Constructor               *
 *    Tau()                              -- Parameterized Constructor         *
 *    Tau()                              -- Constructor with Particle         *
 *    ~Tau()                             -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill                               -- Fill tau from event tree          *
 *    FillFastSim                        -- Fill tau from FastSim tree        *
 *    SetSetNumTracks                    -- Set Num Tracks                    *
 *    GetSetNumTracks                    -- Set Num Tracks                    *
 *    NumTracks                          -- Set Num Tracks                    *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const Tau              *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of Tau class                                          *
 *    Double_t _NumTracks                -- CAL Energy in cone around Tau     *
 *                                                                            *
 *****************************************************************************/

#ifndef tau_h
#define tau_h

#include "SingleTopRootAnalysis/Particles/Recon/Particle.hpp"
#include "SingleTopRootAnalysis/Trees/EventTree.hpp"
#include "SingleTopRootAnalysis/Trees/FastSimTree.hpp"

class Tau: public Particle
{
public:
  // Default Constructor
  Tau();

  // Copy constructor  - Tau
  Tau(const Tau& other);

  // Copy Constructor - Particle
  Tau(const Particle& other);

  // Destructor
  ~Tau();

  // Set all contents to their defaults
  inline void Clear() { Particle::Clear(); _NumTracks = 0.0; };

  // Fill the tau from an EventTree
  Bool_t Fill(EventTree *evtr, Int_t iE, TEnv* config);
  // also fill from a fastsim tree
  Bool_t FillFastSim(FastSimTree *tr, Int_t iE, TEnv* config);

  // Get and Set Number of tracks for this tau
  void SetNumTracks(const Double_t& eop) { _NumTracks = eop; };
  Double_t GetNumTracks() const { return _NumTracks; };
  Double_t NumTracks() const { return _NumTracks; };

  // Overloaded operators:
  // +=
  Tau& operator+=(const Tau& other);
  // +
  Tau operator+(const Tau& other);
  // = const Particle
  Tau& operator=(const Particle& other);
  // = const Tau
  Tau& operator=(const Tau& other);
  // = non-const Tau
  Tau& operator=(Tau& other);

private:
  Double_t _NumTracks;   // number of tracks

  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef(Tau,0)
  ////////////////////////////////////////////////////////////////////////////////    

};


#endif
