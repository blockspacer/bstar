/******************************************************************************
 * MCCQuark.hpp                                                               *
 *                                                                            *
 * Container for a MCCQuark                                                   *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCCQuark class                                  *
 *    MCCQuark()                         -- Default Constructor               *
 *    MCCQuark()                         -- Parameterized Constructor         *
 *    MCCQuark()                         -- Constructor with Particle         *
 *    ~MCCQuark()                        -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill muon from event tree         *
 *    SetDeltaR()                        -- Set DeltaR                        *
 *    GetDeltaR()                        -- Get DeltaR                        *
 *    DeltaR()                           -- Get DeltaR                        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCCQuark         *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCCQuark class                                     *
 *    Double_t _Matched                   -- If BQuark is Matched             *
 *    Double_t _DeltaR                   -- Distance between truth & recon    *
 *                                                                            *
 * History                                                                    *
 *      30 Jun 2009 - Created by P. Ryan                                      *
 *      22 Jul 2009 - Created by J. Holzbauer, based on MCCQuark              *
 *      11 Sep 2009 - Removed overloaded DeltaR() so can use TLorentz verison *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Particles/Truth/MCCQuark.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(MCCQuark)

/******************************************************************************
 * MCCQuark::MCCQuark()                                                       *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  MCCQuark::MCCQuark() : _DeltaR(-1), _Matched(kFALSE)
{

} //MCCQuark()

/******************************************************************************
 * MCCQuark::~MCCQuark()                                                      *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
MCCQuark::~MCCQuark()
{

} //~MCCQuark()


/******************************************************************************
 * MCCQuark::MCCQuark(const MCCQuark& other)                                  *
 *                                                                            *
 * Copy constructor MCCQuark const                                            *
 *                                                                            *
 * Input:  MCCQuark                                                           *
 * Output: None                                                               *
 ******************************************************************************/
MCCQuark::MCCQuark(const MCCQuark& other): MCParticle(other) ,_DeltaR(other.GetDeltaR()),  _Matched(other.Matched())
{

} //MCCQuark()

/******************************************************************************
 * MCCQuark::MCCQuark(const MCCQuark& other)                                  *
 *                                                                            *
 * Copy constructor MCParticle const                                          *
 *                                                                            *
 * Input:  MCCQuark                                                           *
 * Output: None                                                               *
 ******************************************************************************/
MCCQuark::MCCQuark(const MCParticle& other): MCParticle(other),_DeltaR(-1), _Matched(kFALSE)
{

} //MCCQuark

/******************************************************************************
 * MCCQuark& MCCQuark::operator+=(const MCCQuark& other)                      *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  MCCQuark                                                           *
 * Output: MCCQuark                                                           *
 ******************************************************************************/
MCCQuark& MCCQuark::operator+=(const MCCQuark& other)
{
  
  MCParticle::operator+=(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
  if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+=

/******************************************************************************
 * MCCQuark& MCCQuark::operator+(const MCCQuark& other)                       *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  MCCQuark                                                           *
 * Output: MCCQuark                                                           *
 ******************************************************************************/
MCCQuark MCCQuark::operator+(const MCCQuark& other)
{

  MCParticle::operator+(other);
  // also add DeltaR
  SetDeltaR(GetDeltaR() + other.GetDeltaR());
  if( !GetMatched() || other.GetMatched() ) SetMatched(kFALSE);
  else SetMatched(kTRUE);
  return *this;
} //+

/******************************************************************************
 * MCCQuark& MCCQuark::operator=(const BQuark& other)                         *
 *                                                                            *
 * Overload = MCParticle                                                      *
 *                                                                            *
 * Input:  MCParticle                                                         *
 * Output: MCCQuark                                                           *
 ******************************************************************************/
MCCQuark& MCCQuark::operator=(const MCParticle& other)
{
  
  MCParticle::operator=(other);
  SetDeltaR(-1);
  SetMatched(kFALSE);
  return *this;
} // MCParticle

/******************************************************************************
 * MCCQuark& MCCQuark::operator=(const MCCQuark& other)                       *
 *                                                                            *
 * Overload = const MCCQuark                                                  *
 *                                                                            *
 * Input:  MCCQuark                                                           *
 * Output: MCCQuark                                                           *
 ******************************************************************************/
MCCQuark& MCCQuark::operator=(const MCCQuark& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
  SetMatched(other.GetMatched());
  return *this;
} //MCCQuark

/******************************************************************************
 * MCCQuark& MCCQuark::operator=(MCCQuark& other)                             *
 *                                                                            *
 * Overload = non-const MCCQuark                                              *
 *                                                                            *
 * Input:  MCCQuark                                                           *
 * Output: MCCQuark                                                           *
 ******************************************************************************/
MCCQuark& MCCQuark::operator=(MCCQuark& other)
{
  MCParticle::operator=(other);
  SetDeltaR(other.GetDeltaR());
  SetMatched(other.GetMatched());
  return *this;
} //MCCQuark()

/******************************************************************************
 * void MCCQuark::Fill(TruthTree *truthTree, Int_t iE)                        *
 *                                                                            *
 * Fill MCCQuark vector from tree                                             *
 *                                                                            *
 * Input:  Event Tree                                                         *
 * Output: None                                                               *
 ******************************************************************************/
void MCCQuark::FillCQuark(TruthTree *trtr, int iE)
{
  /*
  //MCParticle::Fill(trtr, iE);
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
