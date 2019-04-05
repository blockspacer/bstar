/******************************************************************************
 * MCElectron.cpp                                                             *
 *                                                                            *
 * Container for a MCElectron                                                 *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCElectron class                                *
 *    MCElectron()                           -- Default Constructor           *
 *    MCElectron()                           -- Parameterized Constructor     *
 *    MCElectron()                           -- Constructor with Particle     *
 *    ~MCElectron()                          -- Destructor                    *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill electron from event tree     *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCElectron       *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCElectron class                                   *
 *    Double_t _Matched                   -- If Electron is Matched           *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      13 Apr 2006 - Created by P. Ryan / J. Holzbauer                       *
 *      25 Apr 2007 - Modified by J. Holzbauer, added back DeltaR from tree   *
 *                    set default deltaR to -1                                *
 *      11 Sep 2009 - Removed overloaded DeltaR() so can use TLorentz verison *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Particles/Truth/MCElectron.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(MCElectron)

/******************************************************************************
 * MCElectron::MCElectron()                                                   *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  MCElectron::MCElectron() : _DeltaR(0.0), _Matched(kFALSE)
{

} //MCElectron()

/******************************************************************************
 * MCElectron::~MCElectron()                                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
MCElectron::~MCElectron()
{

} //~MCElectron()


/******************************************************************************
 * MCElectron::MCElectron(const MCElectron& other)                            *
 *                                                                            *
 * Copy constructor MCElectron const                                          *
 *                                                                            *
 * Input:  MCElectron                                                         *
 * Output: None                                                               *
 ******************************************************************************/
MCElectron::MCElectron(const MCElectron& other): MCParticle(other) ,_DeltaR(other.GetDeltaR()),  _Matched(other.Matched())
{

} //MCElectron()

/******************************************************************************
 * MCElectron::MCElectron(const MCElectron& other)                            *
 *                                                                            *
 * Copy constructor MCParticle const                                          *
 *                                                                            *
 * Input:  MCElectron                                                         *
 * Output: None                                                               *
 ******************************************************************************/
MCElectron::MCElectron(const MCParticle& other): MCParticle(other),_DeltaR(-1), _Matched(kFALSE)
{

} //MCElectron

/******************************************************************************
 * MCElectron& MCElectron::operator+=(const MCElectron& other)                *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  MCElectron                                                         *
 * Output: MCElectron                                                         *
 ******************************************************************************/
MCElectron& MCElectron::operator+=(const MCElectron& other)
{
  
  MCParticle::operator+=(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
  if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+=

/******************************************************************************
 * MCElectron& MCElectron::operator+(const MCElectron& other)                 *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  MCElectron                                                         *
 * Output: MCElectron                                                         *
 ******************************************************************************/
MCElectron MCElectron::operator+(const MCElectron& other)
{

  MCParticle::operator+(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
   if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+

/******************************************************************************
 * MCElectron& MCElectron::operator=(const Electron& other)                   *
 *                                                                            *
 * Overload = MCParticle                                                      *
 *                                                                            *
 * Input:  MCParticle                                                         *
 * Output: MCElectron                                                         *
 ******************************************************************************/
MCElectron& MCElectron::operator=(const MCParticle& other)
{
  
  MCParticle::operator=(other);
  SetDeltaR(-1);
  SetMatched(kFALSE);
  return *this;
} // MCParticle

/******************************************************************************
 * MCElectron& MCElectron::operator=(const MCElectron& other)                 *
 *                                                                            *
 * Overload = const MCElectron                                                *
 *                                                                            *
 * Input:  MCElectron                                                         *
 * Output: MCElectron                                                         *
 ******************************************************************************/
MCElectron& MCElectron::operator=(const MCElectron& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
   SetMatched(other.GetMatched());
  return *this;
} //MCElectron

/******************************************************************************
 * MCElectron& MCElectron::operator=(MCElectron& other)                       *
 *                                                                            *
 * Overload = non-const MCElectron                                            *
 *                                                                            *
 * Input:  MCElectron                                                         *
 * Output: MCElectron                                                         *
 ******************************************************************************/
MCElectron& MCElectron::operator=(MCElectron& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
  SetMatched(other.GetMatched());
  return *this;
} //MCElectron()

/******************************************************************************
 * void MCElectron::Fill(TruthTree *truthTree, Int_t iE)                      *
 *                                                                            *
 * Fill MCElectron vector from tree                                           *
 *                                                                            *
 * Input:  Event Tree                                                         *
 * Output: None                                                               *
 ******************************************************************************/
void MCElectron::FillElectron(TruthTree *trtr, int iE)
{
  /*
  //MCParticle::Fill(trtr, iE);

  ///////////////////////////////////////////////////////////////
  // Version 13.0.30 of ntuples doesn't have these variables.
  // Set to -999 for now
  ///////////////////////////////////////////////////////////////
  _DeltaR =  -999;
  _Matched  = -999.0;

  //_DeltaR = trtr -> Tru_El_DeltaR -> operator[](iE);
  //_Matched  = trtr -> Tru_El_Matched -> operator[](iE);
  //std::cout<<trtr -> Electron_pt     -> operator[](iE)<<endl;
  Double_t pPt    = trtr -> Electron_pt     -> operator[](iE)/1000.;
  Double_t pEta    = trtr -> Electron_eta    -> operator[](iE);
  Double_t pPhi    = trtr -> Electron_phi    -> operator[](iE);
  Double_t pM      = trtr -> Electron_m      -> operator[](iE)/1000.;
  Double_t pCharge = trtr -> Electron_charge -> operator[](iE);
  
  _PdgId   = static_cast<Int_t>(trtr -> Electron_pdgId  -> operator[](iE));
  _Status  = static_cast<Int_t>(trtr -> Electron_status -> operator[](iE));
  _BarCode = trtr -> Electron_barcode->operator[](iE);

  SetPtEtaPhiM(pPt,pEta,pPhi,pM);
  SetCharge(pCharge);
  */
} //Fill()
