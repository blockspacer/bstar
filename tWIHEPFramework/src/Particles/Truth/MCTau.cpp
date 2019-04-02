/******************************************************************************
 * MCTau.cpp                                                                  *
 *                                                                            *
 * Container for a MCTau                                                      *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCTau class                                     *
 *    MCTau()                           -- Default Constructor                *
 *    MCTau()                           -- Parameterized Constructor          *
 *    MCTau()                           -- Constructor with Particle          *
 *    ~MCTau()                          -- Destructor                         *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill tau from event tree          *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCTau            *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCTau class                                        *
 *    Double_t _Matched                   -- If Tau is Matched                *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      13 Apr 2006 - Created by P. Ryan / J. Holzbauer                       *
 *      25 Apr 2007 - Modified by J. Holzbauer, added back DeltaR from tree   *
 *                    set default deltaR to -1                                *
 *      11 Sep 2009 - Removed overloaded DeltaR() so can use TLorentz verison *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Particles/Truth/MCTau.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(MCTau)

/******************************************************************************
 * MCTau::MCTau()                                                             *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  MCTau::MCTau() : _DeltaR(-1), _Matched(kFALSE)
{

} //MCTau()

/******************************************************************************
 * MCTau::~MCTau()                                                            *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
MCTau::~MCTau()
{

} //~MCTau()


/******************************************************************************
 * MCTau::MCTau(const MCTau& other)                                           *
 *                                                                            *
 * Copy constructor MCTau const                                               *
 *                                                                            *
 * Input:  MCTau                                                              *
 * Output: None                                                               *
 ******************************************************************************/
MCTau::MCTau(const MCTau& other): MCParticle(other) ,_DeltaR(other.GetDeltaR()),  _Matched(other.Matched())
{

} //MCTau()

/******************************************************************************
 * MCTau::MCTau(const MCTau& other)                                           *
 *                                                                            *
 * Copy constructor MCParticle const                                          *
 *                                                                            *
 * Input:  MCTau                                                              *
 * Output: None                                                               *
 ******************************************************************************/
MCTau::MCTau(const MCParticle& other): MCParticle(other),_DeltaR(-1), _Matched(kFALSE)
{

} //MCTau

/******************************************************************************
 * MCTau& MCTau::operator+=(const MCTau& other)                               *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  MCTau                                                              *
 * Output: MCTau                                                              *
 ******************************************************************************/
MCTau& MCTau::operator+=(const MCTau& other)
{
  
  MCParticle::operator+=(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
  if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+=

/******************************************************************************
 * MCTau& MCTau::operator+(const MCTau& other)                                *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  MCTau                                                              *
 * Output: MCTau                                                              *
 ******************************************************************************/
MCTau MCTau::operator+(const MCTau& other)
{

  MCParticle::operator+(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
   if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+

/******************************************************************************
 * MCTau& MCTau::operator=(const Tau& other)                                  *
 *                                                                            *
 * Overload = MCParticle                                                      *
 *                                                                            *
 * Input:  MCParticle                                                         *
 * Output: MCTau                                                              *
 ******************************************************************************/
MCTau& MCTau::operator=(const MCParticle& other)
{
  
  MCParticle::operator=(other);
  SetDeltaR(-1);
  SetMatched(kFALSE);
  return *this;
} // MCParticle

/******************************************************************************
 * MCTau& MCTau::operator=(const MCTau& other)                                *
 *                                                                            *
 * Overload = const MCTau                                                     *
 *                                                                            *
 * Input:  MCTau                                                              *
 * Output: MCTau                                                              *
 ******************************************************************************/
MCTau& MCTau::operator=(const MCTau& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
   SetMatched(other.GetMatched());
  return *this;
} //MCTau

/******************************************************************************
 * MCTau& MCTau::operator=(MCTau& other)                                      *
 *                                                                            *
 * Overload = non-const MCTau                                                 *
 *                                                                            *
 * Input:  MCTau                                                              *
 * Output: MCTau                                                              *
 ******************************************************************************/
MCTau& MCTau::operator=(MCTau& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
  SetMatched(other.GetMatched());
  return *this;
} //MCTau()

/******************************************************************************
 * void MCTau::Fill(TruthTree *truthTree, Int_t iE)                           *
 *                                                                            *
 * Fill MCTau vector from tree                                                *
 *                                                                            *
 * Input:  Event Tree                                                         *
 * Output: None                                                               *
 ******************************************************************************/
void MCTau::FillTau(TruthTree *trtr, int iE)
{
  /*
  //MCParticle::Fill(trtr, iE);

  ///////////////////////////////////////////////////////////////
  // Version 13.0.30 of ntuples doesn't have these variables.
  // Set to -999 for now
  ///////////////////////////////////////////////////////////////

  _DeltaR = -999.0;
  _Matched  = -999;

  //_DeltaR = trtr -> Tru_Tau_DeltaR -> operator[](iE);
  //_Matched  = trtr -> Tru_Tau_Matched -> operator[](iE);

  Double_t pPt     = trtr -> Tau_pt    -> operator[](iE)/1000.;
  Double_t pEta    = trtr -> Tau_eta    -> operator[](iE);
  Double_t pPhi    = trtr -> Tau_phi    -> operator[](iE);
  Double_t pM      = trtr -> Tau_m      -> operator[](iE)/1000.;
  Double_t pCharge = trtr -> Tau_charge -> operator[](iE);
  
  _PdgId   = static_cast<Int_t>(trtr -> Tau_pdgId  -> operator[](iE));
  _Status  = static_cast<Int_t>(trtr -> Tau_status -> operator[](iE));
  _BarCode = trtr -> Tau_barcode->operator[](iE);

 SetPtEtaPhiM(pPt,pEta,pPhi,pM);
    SetCharge(pCharge);
  */
} //Fill()
