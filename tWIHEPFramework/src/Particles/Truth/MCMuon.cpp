/******************************************************************************
 * MCMuon.cpp                                                                 *
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
 *      11 Sep 2009 - Removed overloaded DeltaR() so can use TLorentz verison *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Particles/Truth/MCMuon.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(MCMuon)

/******************************************************************************
 * MCMuon::MCMuon()                                                           *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  MCMuon::MCMuon() : _DeltaR(-1), _Matched(kFALSE)
{

} //MCMuon()

/******************************************************************************
 * MCMuon::~MCMuon()                                                          *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
MCMuon::~MCMuon()
{

} //~MCMuon()


/******************************************************************************
 * MCMuon::MCMuon(const MCMuon& other)                                        *
 *                                                                            *
 * Copy constructor MCMuon const                                              *
 *                                                                            *
 * Input:  MCMuon                                                             *
 * Output: None                                                               *
 ******************************************************************************/
MCMuon::MCMuon(const MCMuon& other): MCParticle(other) ,_DeltaR(other.GetDeltaR()),  _Matched(other.Matched())
{

} //MCMuon()

/******************************************************************************
 * MCMuon::MCMuon(const MCMuon& other)                                        *
 *                                                                            *
 * Copy constructor MCParticle const                                          *
 *                                                                            *
 * Input:  MCMuon                                                             *
 * Output: None                                                               *
 ******************************************************************************/
MCMuon::MCMuon(const MCParticle& other): MCParticle(other),_DeltaR(-1), _Matched(kFALSE)
{

} //MCMuon

/******************************************************************************
 * MCMuon& MCMuon::operator+=(const MCMuon& other)                            *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  MCMuon                                                             *
 * Output: MCMuon                                                             *
 ******************************************************************************/
MCMuon& MCMuon::operator+=(const MCMuon& other)
{
  
  MCParticle::operator+=(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
  if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+=

/******************************************************************************
 * MCMuon& MCMuon::operator+(const MCMuon& other)                             *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  MCMuon                                                             *
 * Output: MCMuon                                                             *
 ******************************************************************************/
MCMuon MCMuon::operator+(const MCMuon& other)
{

  MCParticle::operator+(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
   if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+

/******************************************************************************
 * MCMuon& MCMuon::operator=(const Muon& other)                               *
 *                                                                            *
 * Overload = MCParticle                                                      *
 *                                                                            *
 * Input:  MCParticle                                                         *
 * Output: MCMuon                                                             *
 ******************************************************************************/
MCMuon& MCMuon::operator=(const MCParticle& other)
{
  
  MCParticle::operator=(other);
  SetDeltaR(-1);
  SetMatched(kFALSE);
  return *this;
} // MCParticle

/******************************************************************************
 * MCMuon& MCMuon::operator=(const MCMuon& other)                             *
 *                                                                            *
 * Overload = const MCMuon                                                    *
 *                                                                            *
 * Input:  MCMuon                                                             *
 * Output: MCMuon                                                             *
 ******************************************************************************/
MCMuon& MCMuon::operator=(const MCMuon& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
   SetMatched(other.GetMatched());
  return *this;
} //MCMuon

/******************************************************************************
 * MCMuon& MCMuon::operator=(MCMuon& other)                                   *
 *                                                                            *
 * Overload = non-const MCMuon                                                *
 *                                                                            *
 * Input:  MCMuon                                                             *
 * Output: MCMuon                                                             *
 ******************************************************************************/
MCMuon& MCMuon::operator=(MCMuon& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
  SetMatched(other.GetMatched());
  return *this;
} //MCMuon()

/******************************************************************************
 * void MCMuon::Fill(TruthTree *truthTree, Int_t iE)                          *
 *                                                                            *
 * Fill MCMuon vector from tree                                               *
 *                                                                            *
 * Input:  Event Tree                                                         *
 * Output: None                                                               *
 ******************************************************************************/
void MCMuon::FillMuon(TruthTree *trtr, int iE)
{
  /*
  //  MCParticle::Fill(trtr, iE);

   ///////////////////////////////////////////////////////////////
  // Version 13.0.30 of ntuples doesn't have these variables.
  // Set to -999 for now
  ///////////////////////////////////////////////////////////////

  _DeltaR = -999.0;
  _Matched  = -999;


  //  _DeltaR = trtr -> Tru_Mu_DeltaR -> operator[](iE);
  //  _Matched  = trtr -> Tru_Mu_Matched -> operator[](iE);

  Double_t pPt     = trtr -> Muon_pt    -> operator[](iE)/1000.;
  Double_t pEta    = trtr -> Muon_eta    -> operator[](iE);
  Double_t pPhi    = trtr -> Muon_phi    -> operator[](iE);
  Double_t pM      = trtr -> Muon_m      -> operator[](iE)/1000.;
  Double_t pCharge = trtr -> Muon_charge -> operator[](iE);
  
  _PdgId   = static_cast<Int_t>(trtr -> Muon_pdgId  -> operator[](iE));
  _Status  = static_cast<Int_t>(trtr -> Muon_status -> operator[](iE));
  _BarCode = trtr -> Muon_barcode->operator[](iE);


 SetPtEtaPhiM(pPt,pEta,pPhi,pM);
    SetCharge(pCharge);
  */
} //Fill()
