/******************************************************************************
 * MCJet.cpp                                                                  *
 *                                                                            *
 * Container for a MCJet                                                      *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCJet class                                     *
 *    MCJet()                           -- Default Constructor                *
 *    MCJet()                           -- Parameterized Constructor          *
 *    MCJet()                           -- Constructor with Particle          *
 *    ~MCJet()                          -- Destructor                         *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill jet from event tree          *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCJet            *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCJet class                                        *
 *    Double_t _Matched                   -- If Jet is Matched                *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      13 Apr 2006 - Created by P. Ryan / J. Holzbauer                       *
 *      25 Apr 2007 - Modified by J. Holzbauer, added back DeltaR from tree   *
 *                    set default deltaR to -1                                *
 *      11 Sep 2009 - Removed overloaded DeltaR() so can use TLorentz verison *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Particles/Truth/MCJet.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(MCJet)

/******************************************************************************
 * MCJet::MCJet()                                                             *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  MCJet::MCJet() : _DeltaR(-1), _Matched(kFALSE)
{

} //MCJet()

/******************************************************************************
 * MCJet::~MCJet()                                                            *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
MCJet::~MCJet()
{

} //~MCJet()


/******************************************************************************
 * MCJet::MCJet(const MCJet& other)                                           *
 *                                                                            *
 * Copy constructor MCJet const                                               *
 *                                                                            *
 * Input:  MCJet                                                              *
 * Output: None                                                               *
 ******************************************************************************/
MCJet::MCJet(const MCJet& other): MCParticle(other) ,_DeltaR(other.GetDeltaR()),  _Matched(other.Matched())
{

} //MCJet()

/******************************************************************************
 * MCJet::MCJet(const MCJet& other)                                           *
 *                                                                            *
 * Copy constructor MCParticle const                                          *
 *                                                                            *
 * Input:  MCJet                                                              *
 * Output: None                                                               *
 ******************************************************************************/
MCJet::MCJet(const MCParticle& other): MCParticle(other),_DeltaR(-1), _Matched(kFALSE)
{

} //MCJet

/******************************************************************************
 * MCJet& MCJet::operator+=(const MCJet& other)                               *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  MCJet                                                              *
 * Output: MCJet                                                              *
 ******************************************************************************/
MCJet& MCJet::operator+=(const MCJet& other)
{
  
  MCParticle::operator+=(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
  if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+=

/******************************************************************************
 * MCJet& MCJet::operator+(const MCJet& other)                                *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  MCJet                                                              *
 * Output: MCJet                                                              *
 ******************************************************************************/
MCJet MCJet::operator+(const MCJet& other)
{

  MCParticle::operator+(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
  if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+

/******************************************************************************
 * MCJet& MCJet::operator=(const Jet& other)                                  *
 *                                                                            *
 * Overload = MCParticle                                                      *
 *                                                                            *
 * Input:  MCParticle                                                         *
 * Output: MCJet                                                              *
 ******************************************************************************/
MCJet& MCJet::operator=(const MCParticle& other)
{
  
  MCParticle::operator=(other);
  SetDeltaR(-1);
  SetMatched(kFALSE);
  return *this;
} // MCParticle

/******************************************************************************
 * MCJet& MCJet::operator=(const MCJet& other)                                *
 *                                                                            *
 * Overload = const MCJet                                                     *
 *                                                                            *
 * Input:  MCJet                                                              *
 * Output: MCJet                                                              *
 ******************************************************************************/
MCJet& MCJet::operator=(const MCJet& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
  SetMatched(other.GetMatched());
  return *this;
} //MCJet

/******************************************************************************
 * MCJet& MCJet::operator=(MCJet& other)                                      *
 *                                                                            *
 * Overload = non-const MCJet                                                 *
 *                                                                            *
 * Input:  MCJet                                                              *
 * Output: MCJet                                                              *
 ******************************************************************************/
MCJet& MCJet::operator=(MCJet& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
  SetMatched(other.GetMatched());
  return *this;
} //MCJet()

/******************************************************************************
 * void MCJet::Fill(TruthTree *truthTree, Int_t iE)                           *
 *                                                                            *
 * Fill MCJet vector from tree                                                *
 *                                                                            *
 * Input:  Event Tree                                                         *
 * Output: None                                                               *
 ******************************************************************************/
void MCJet::FillJet(TruthTree *trtr, int iE)
{
  /*
  // MCParticle::Fill(trtr, iE);
  // _Matched = kFALSE;
  
  ///////////////////////////////////////////////////////
  // Tru_PJet_Matched and Tru_PJet_DeltaR not in 
  // TruthTree for version 13.0.30 - set to -999
  /////////////////////////////////////////////////////
  _Matched  = -999;
  _DeltaR = -999.0;

  //_Matched  = trtr -> Tru_PJet_Matched -> operator[](iE);
  //_DeltaR = trtr -> Tru_PJet_DeltaR -> operator[](iE);
  
  Double_t pPt     = trtr -> Jet_pt    -> operator[](iE)/1000.;
  Double_t pEta    = trtr -> Jet_eta    -> operator[](iE);
  Double_t pPhi    = trtr -> Jet_phi    -> operator[](iE);
  Double_t pM      = trtr -> Jet_m      -> operator[](iE)/1000.;
  Double_t pCharge = trtr -> Jet_charge -> operator[](iE);
  
  _PdgId   = static_cast<Int_t>(trtr -> Jet_pdgId  -> operator[](iE));
  //  _Status  = static_cast<Int_t>(trtr -> Jet_status -> operator[](iE));
  _Status  = 0;
  _BarCode = trtr -> Jet_Barcode->operator[](iE);

 SetPtEtaPhiM(pPt,pEta,pPhi,pM);
    SetCharge(pCharge);
  */
} //Fill()
