/******************************************************************************
 * MCTop.cpp                                                                  *
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

#include "SingleTopRootAnalysis/Particles/Truth/MCTop.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(MCTop)

/******************************************************************************
 * MCTop::MCTop()                                                             *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  MCTop::MCTop()
{

} //MCTop()

/******************************************************************************
 * MCTop::~MCTop()                                                            *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
MCTop::~MCTop()
{

} //~MCTop()


/******************************************************************************
 * MCTop::MCTop(const MCTop& other)                                           *
 *                                                                            *
 * Copy constructor MCTop const                                               *
 *                                                                            *
 * Input:  MCTop                                                              *
 * Output: None                                                               *
 ******************************************************************************/
MCTop::MCTop(const MCTop& other): MCParticle(other)
{

} //MCTop()

/******************************************************************************
 * MCTop::MCTop(const MCTop& other)                                           *
 *                                                                            *
 * Copy constructor MCParticle const                                          *
 *                                                                            *
 * Input:  MCTop                                                              *
 * Output: None                                                               *
 ******************************************************************************/
MCTop::MCTop(const MCParticle& other): MCParticle(other)
{ 

} //MCTop

/******************************************************************************
 * MCTop& MCTop::operator+=(const MCTop& other)                               *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  MCTop                                                              *
 * Output: MCTop                                                              *
 ******************************************************************************/
MCTop& MCTop::operator+=(const MCTop& other)
{
  
  MCParticle::operator+=(other);
   return *this;
} //+=

/******************************************************************************
 * MCTop& MCTop::operator+(const MCTop& other)                                *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  MCTop                                                              *
 * Output: MCTop                                                              *
 ******************************************************************************/
MCTop MCTop::operator+(const MCTop& other)
{

  MCParticle::operator+(other);
  return *this;
} //+

/******************************************************************************
 * MCTop& MCTop::operator=(const Top& other)                                  *
 *                                                                            *
 * Overload = MCParticle                                                      *
 *                                                                            *
 * Input:  MCParticleicle::Clear();                                           *
 * Output: MCTop                                                              *
 ******************************************************************************/
MCTop& MCTop::operator=(const MCParticle& other)
{
  
  MCParticle::operator=(other);
  return *this;
} // MCParticle

/******************************************************************************
 * MCTop& MCTop::operator=(const MCTop& other)                                *
 *                                                                            *
 * Overload = const MCTop                                                     *
 *                                                                            *
 * Input:  MCTop                                                              *
 * Output: MCTop                                                              *
 ******************************************************************************/
MCTop& MCTop::operator=(const MCTop& other)
{
  MCParticle::operator=(other);
  return *this;
} //MCTop

/******************************************************************************
 * MCTop& MCTop::operator=(MCTop& other)                                      *
 *                                                                            *
 * Overload = non-const MCTop                                                 *
 *                                                                            *
 * Input:  MCTop                                                              *
 * Output: MCTop                                                              *
 ******************************************************************************/
MCTop& MCTop::operator=(MCTop& other)
{
  MCParticle::operator=(other);
  return *this;
} //MCTop()

/******************************************************************************
 * void MCTop::Fill(TruthTree *truthTree, Int_t iE)                           *
 *                                                                            *
 * Fill MCTop vector from tree                                                *
 *                                                                            *
 * Input:  Event Tree                                                         *
 * Output: None                                                               *
 ******************************************************************************/
void MCTop::FillTop(TruthTree *trtr, int iE)
{
  /**
  //MCParticle::Fill(trtr, iE);

  Double_t pPt     = trtr -> TopDecay_pt    -> operator[](iE)/1000.;
  Double_t pEta    = trtr -> TopDecay_eta    -> operator[](iE);
  Double_t pPhi    = trtr -> TopDecay_phi    -> operator[](iE);
  Double_t pM      = trtr -> TopDecay_m      -> operator[](iE)/1000.;
  Double_t pCharge = trtr -> TopDecay_charge -> operator[](iE);
  
  _PdgId   = static_cast<Int_t>(trtr -> TopDecay_pdgId  -> operator[](iE));
  _Status  = static_cast<Int_t>(trtr -> TopDecay_status -> operator[](iE));
  _BarCode = trtr -> TopDecay_barcode->operator[](iE);
 SetPtEtaPhiM(pPt,pEta,pPhi,pM);
    SetCharge(pCharge);
  **/

} //Fill()
