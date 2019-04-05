/******************************************************************************
 * MCBQuark.hpp                                                               *
 *                                                                            *
 * Container for a MCBQuark                                                   *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCBQuark class                                  *
 *    MCBQuark()                         -- Default Constructor               *
 *    MCBQuark()                         -- Parameterized Constructor         *
 *    MCBQuark()                         -- Constructor with Particle         *
 *    ~MCBQuark()                        -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill muon from event tree         *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCBQuark         *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCBQuark class                                     *
 *    Double_t _Matched                   -- If BQuark is Matched             *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      30 Jun 2009 - Created by P. Ryan                                      *
 *                                                                            *
 *      11 Sep 2009 - Removed overloaded DeltaR() so can use TLorentz verison *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Particles/Truth/MCBQuark.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(MCBQuark)

/******************************************************************************
 * MCBQuark::MCBQuark()                                                       *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  MCBQuark::MCBQuark() : _DeltaR(-1), _Matched(kFALSE)
{

} //MCBQuark()

/******************************************************************************
 * MCBQuark::~MCBQuark()                                                      *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
MCBQuark::~MCBQuark()
{

} //~MCBQuark()


/******************************************************************************
 * MCBQuark::MCBQuark(const MCBQuark& other)                                  *
 *                                                                            *
 * Copy constructor MCBQuark const                                            *
 *                                                                            *
 * Input:  MCBQuark                                                           *
 * Output: None                                                               *
 ******************************************************************************/
MCBQuark::MCBQuark(const MCBQuark& other): MCParticle(other) ,_DeltaR(other.GetDeltaR()),  _Matched(other.Matched())
{

} //MCBQuark()

/******************************************************************************
 * MCBQuark::MCBQuark(const MCBQuark& other)                                  *
 *                                                                            *
 * Copy constructor MCParticle const                                          *
 *                                                                            *
 * Input:  MCBQuark                                                           *
 * Output: None                                                               *
 ******************************************************************************/
MCBQuark::MCBQuark(const MCParticle& other): MCParticle(other),_DeltaR(-1), _Matched(kFALSE)
{

} //MCBQuark

/******************************************************************************
 * MCBQuark& MCBQuark::operator+=(const MCBQuark& other)                      *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  MCBQuark                                                           *
 * Output: MCBQuark                                                           *
 ******************************************************************************/
MCBQuark& MCBQuark::operator+=(const MCBQuark& other)
{
  
  MCParticle::operator+=(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
  if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+=

/******************************************************************************
 * MCBQuark& MCBQuark::operator+(const MCBQuark& other)                       *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  MCBQuark                                                           *
 * Output: MCBQuark                                                           *
 ******************************************************************************/
MCBQuark MCBQuark::operator+(const MCBQuark& other)
{

  MCParticle::operator+(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
  if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+

/******************************************************************************
 * MCBQuark& MCBQuark::operator=(const BQuark& other)                         *
 *                                                                            *
 * Overload = MCParticle                                                      *
 *                                                                            *
 * Input:  MCParticle                                                         *
 * Output: MCBQuark                                                           *
 ******************************************************************************/
MCBQuark& MCBQuark::operator=(const MCParticle& other)
{
  
  MCParticle::operator=(other);
  SetDeltaR(-1);
  SetMatched(kFALSE);
  return *this;
} // MCParticle

/******************************************************************************
 * MCBQuark& MCBQuark::operator=(const MCBQuark& other)                       *
 *                                                                            *
 * Overload = const MCBQuark                                                  *
 *                                                                            *
 * Input:  MCBQuark                                                           *
 * Output: MCBQuark                                                           *
 ******************************************************************************/
MCBQuark& MCBQuark::operator=(const MCBQuark& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
  SetMatched(other.GetMatched());
  return *this;
} //MCBQuark

/******************************************************************************
 * MCBQuark& MCBQuark::operator=(MCBQuark& other)                             *
 *                                                                            *
 * Overload = non-const MCBQuark                                              *
 *                                                                            *
 * Input:  MCBQuark                                                           *
 * Output: MCBQuark                                                           *
 ******************************************************************************/
MCBQuark& MCBQuark::operator=(MCBQuark& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
  SetMatched(other.GetMatched());
  return *this;
} //MCBQuark()

/******************************************************************************
 * void MCBQuark::Fill(TruthTree *truthTree, Int_t iE)                        *
 *                                                                            *
 * Fill MCBQuark vector from tree                                             *
 *                                                                            *
 * Input:  Event Tree                                                         *
 * Output: None                                                               *
 ******************************************************************************/
void MCBQuark::FillBQuark(TruthTree *trtr, int iE)
{
  /*
  // MCParticle::Fill(trtr, iE);
  // _Matched = kFALSE;
  
  ///////////////////////////////////////////////////////
  // Tru_PBQuark_Matched and Tru_PBQuark_DeltaR not in 
  // TruthTree for version 13.0.30 - set to -999
  /////////////////////////////////////////////////////
  _Matched  = -999;
  _DeltaR = -999.0;

  //_Matched  = trtr -> Tru_PBQuark_Matched -> operator[](iE);
  //_DeltaR = trtr -> Tru_PBQuark_DeltaR -> operator[](iE);

  Double_t pPt     = trtr -> BCQuark_pt    -> operator[](iE)/1000.;
  Double_t pEta    = trtr -> BCQuark_eta    -> operator[](iE);
  Double_t pPhi    = trtr -> BCQuark_phi    -> operator[](iE);
  Double_t pM      = trtr -> BCQuark_m      -> operator[](iE)/1000.;
  Double_t pCharge = trtr -> BCQuark_charge -> operator[](iE);
  
  _PdgId   = static_cast<Int_t>(trtr -> BCQuark_pdgId  -> operator[](iE));
  _Status  = static_cast<Int_t>(trtr -> BCQuark_status -> operator[](iE));
  _BarCode = trtr -> BCQuark_barcode->operator[](iE);

   SetPtEtaPhiM(pPt,pEta,pPhi,pM);
    SetCharge(pCharge);
  */
} //Fill()
