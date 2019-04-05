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
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst                              *
 *      20 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      20 Mar 2007 - RS: Add filling from fastsim tree                       *
 *      11 Jan 2008 - RS: Fill only if object passes ID cuts                  *
 *****************************************************************************/
#include "SingleTopRootAnalysis/Particles/Recon/Tau.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(Tau)
/******************************************************************************
 * Tau::Tau()                                                                 *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
Tau::Tau() : Particle::Particle(),_NumTracks(0.0)
{

} //Tau()

/******************************************************************************
 * Tau::~Tau()                                                              *
 *                                                                            *
 * Destructor                                                                 *  
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
Tau::~Tau()
{

} //~Tau()


/******************************************************************************
 * Tau::Tau(const Tau& other)                                              *
 *                                                                            *
 * Copy constructor Tau const                                                *
 *                                                                            *
 * Input:  Tau                                                               *
 * Output: None                                                               *
 ******************************************************************************/
Tau::Tau(const Tau& other): Particle(other),_NumTracks(other.GetNumTracks())
{

} //Tau()

/******************************************************************************
 * Tau::Tau(const Tau& other)                                                 *
 *                                                                            *
 * Copy constructor Particle const                                            *
 *                                                                            *
 * Input:  Tau                                                                *
 * Output: None                                                               *
 ******************************************************************************/
Tau::Tau(const Particle& other): Particle(other),_NumTracks(0.0)
{

} //Tau() particle

/******************************************************************************
 * Tau& Tau::operator+=(const Tau& other)                                     *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  Tau                                                                *
 * Output: Tau                                                                *
 ******************************************************************************/
Tau& Tau::operator+=(const Tau& other)
{
  Particle::operator+=(other);
  // Also add NumTracks
  SetNumTracks(GetNumTracks() + other.GetNumTracks());
  return *this;
} //+=

/******************************************************************************
 * Tau& Tau::operator+(const Tau& other)                                      *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  Tau                                                                *
 * Output: Tau                                                                *
 ******************************************************************************/
Tau Tau::operator+(const Tau& other)
{

  Particle particleTemp = Particle::operator+(other);
  Tau tauTemp = particleTemp;
  
  // I don't know how to add the NumTracks for two different taus.
  // For now, add them
  _NumTracks += other.GetNumTracks();
  tauTemp.SetNumTracks( _NumTracks );

  // Return the temp particle
  return tauTemp;
} //+

/******************************************************************************
 * Tau& Tau::operator=(const Particle& other)                                 *
 *                                                                            *
 * Overload = Particle                                                        *
 *                                                                            *
 * Input:  Particle                                                           *
 * Output: Tau                                                                *
 ******************************************************************************/
Tau& Tau::operator=(const Particle& other)
{
  Particle::operator=(other);
  SetNumTracks(0.0);
  return *this;
} // = Particle

/******************************************************************************
 * Tau& Tau::operator=(const Tau& other)                                      *
 *                                                                            *
 * Overload = const Tau                                                       *
 *                                                                            *
 * Input:  Tau                                                                *
 * Output: Tau                                                                *
 ******************************************************************************/
Tau& Tau::operator=(const Tau& other)
{
  
  Particle::operator=(other);
  SetNumTracks(other.GetNumTracks());
  return *this;
} //= const tau

/******************************************************************************
 * Tau& Tau::operator=(const Tau& other)                                      *
 *                                                                            *
 * Overload = non-const Tau                                                   *
 *                                                                            *
 * Input:  Tau                                                                *
 * Output: Tau                                                                *
 ******************************************************************************/
Tau& Tau::operator=(Tau& other)
{
  
  Particle::operator=(other);
  SetNumTracks(other.GetNumTracks());
  return *this;
} //= tau non-const

/******************************************************************************         
 * void Tau::Fill(EventTree *evtr, Int_t iE)                                  *         
 *                                                                            *         
 * Fill Tau vector from tree                                                  *         
 *                                                                            *         
 * Input:  Event Tree                                                         *         
 * Output: kTRUE if object passes ID cuts                                     *         
 ******************************************************************************/
Bool_t Tau::Fill(EventTree *evtr,int iE,TEnv* config)
{
  /*
    Double_t tauPt     = evtr -> TauJet_pt       -> operator[](iE) / 1000.0;
    Double_t tauEta    = evtr -> TauJet_eta      -> operator[](iE);
    Double_t tauPhi    = evtr -> TauJet_phi      -> operator[](iE);
    Double_t tauE      = evtr -> TauJet_e        -> operator[](iE) / 1000.0;
    _NumTracks         = evtr -> TauJet_numTrack -> operator[](iE);
    Double_t tauCharge = evtr -> TauJet_charge   -> operator[](iE);

    SetPtEtaPhiE(tauPt, tauEta, tauPhi, tauE);
    SetCharge(tauCharge);
    //
    // now check electron object ID cuts
    double maxEta=config->GetValue("ObjectID.Tau.MaxEta",100.);
    double minPt=config->GetValue("ObjectID.Tau.MinPt",0.);

    if(tauPt<minPt) return kFALSE;
    if(TMath::Abs(tauEta)>maxEta) return kFALSE;

    return kTRUE;
  */
} //Fill()

/******************************************************************************         
 * void Tau::FillFastSim(FastSimTree *evtr, Int_t iE)                         *         
 *                                                                            *         
 * Fill Tau vector from tree                                                  *         
 *                                                                            *         
 * Input:  FastSim Tree                                                       *         
 * Output: kTRUE if object passes ID cuts                                     *         
 ******************************************************************************/
Bool_t Tau::FillFastSim(FastSimTree *tr,int iE,TEnv* config)
{

    Double_t tauPt     = tr -> Tau_p_T      -> operator[](iE) / 1000.0;
    Double_t tauEta    = tr -> Tau_eta      -> operator[](iE);
    Double_t tauPhi    = tr -> Tau_phi      -> operator[](iE);
    Double_t tauE      = tr -> Tau_E        -> operator[](iE) / 1000.0;
    _NumTracks         = 0; // not available in fast sim
    Double_t tauCharge = tr -> Tau_charge   -> operator[](iE);

    SetPtEtaPhiE(tauPt, tauEta, tauPhi, tauE);
    SetCharge(tauCharge);
    //
    // now check electron object ID cuts
    double maxEta=config->GetValue("ObjectID.Tau.MaxEta",100.);
    double minPt=config->GetValue("ObjectID.Tau.MinPt",0.);

    if(tauPt<minPt) return kFALSE;
    if(TMath::Abs(tauEta)>maxEta) return kFALSE;

    return kTRUE;
} //FillFastSim()


