/******************************************************************************
 * MCNeutrino.cpp                                                                  *
 *                                                                            *
 * Container for a MCNeutrino                                                      *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCNeutrino class                                     *
 *    MCNeutrino()                           -- Default Constructor                *
 *    MCNeutrino()                           -- Parameterized Constructor          *
 *    MCNeutrino()                           -- Constructor with Particle          *
 *    ~MCNeutrino()                          -- Destructor                         *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill top from event tree          *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCNeutrino            *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCNeutrino class                                        *
 *                                                                            *
 * History                                                                    *
 *      15 May 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Particles/Truth/MCNeutrino.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(MCNeutrino)

/******************************************************************************
 * MCNeutrino::MCNeutrino()                                                             *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  MCNeutrino::MCNeutrino()
{

} //MCNeutrino()

/******************************************************************************
 * MCNeutrino::~MCNeutrino()                                                            *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
MCNeutrino::~MCNeutrino()
{

} //~MCNeutrino()


/******************************************************************************
 * MCNeutrino::MCNeutrino(const MCNeutrino& other)                                           *
 *                                                                            *
 * Copy constructor MCNeutrino const                                               *
 *                                                                            *
 * Input:  MCNeutrino                                                              *
 * Output: None                                                               *
 ******************************************************************************/
MCNeutrino::MCNeutrino(const MCNeutrino& other): MCParticle(other)
{

} //MCNeutrino()

/******************************************************************************
 * MCNeutrino::MCNeutrino(const MCNeutrino& other)                                           *
 *                                                                            *
 * Copy constructor MCParticle const                                          *
 *                                                                            *
 * Input:  MCNeutrino                                                              *
 * Output: None                                                               *
 ******************************************************************************/
MCNeutrino::MCNeutrino(const MCParticle& other): MCParticle(other)
{ 

} //MCNeutrino

/******************************************************************************
 * MCNeutrino& MCNeutrino::operator+=(const MCNeutrino& other)                               *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  MCNeutrino                                                              *
 * Output: MCNeutrino                                                              *
 ******************************************************************************/
MCNeutrino& MCNeutrino::operator+=(const MCNeutrino& other)
{
  
  MCParticle::operator+=(other);
   return *this;
} //+=

/******************************************************************************
 * MCNeutrino& MCNeutrino::operator+(const MCNeutrino& other)                                *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  MCNeutrino                                                              *
 * Output: MCNeutrino                                                              *
 ******************************************************************************/
MCNeutrino MCNeutrino::operator+(const MCNeutrino& other)
{

  MCParticle::operator+(other);
  return *this;
} //+

/******************************************************************************
 * MCNeutrino& MCNeutrino::operator=(const Neutrino& other)                                  *
 *                                                                            *
 * Overload = MCParticle                                                      *
 *                                                                            *
 * Input:  MCParticleicle::Clear();                                           *
 * Output: MCNeutrino                                                              *
 ******************************************************************************/
MCNeutrino& MCNeutrino::operator=(const MCParticle& other)
{
  
  MCParticle::operator=(other);
  return *this;
} // MCParticle

/******************************************************************************
 * MCNeutrino& MCNeutrino::operator=(const MCNeutrino& other)                                *
 *                                                                            *
 * Overload = const MCNeutrino                                                     *
 *                                                                            *
 * Input:  MCNeutrino                                                              *
 * Output: MCNeutrino                                                              *
 ******************************************************************************/
MCNeutrino& MCNeutrino::operator=(const MCNeutrino& other)
{
  MCParticle::operator=(other);
  return *this;
} //MCNeutrino

/******************************************************************************
 * MCNeutrino& MCNeutrino::operator=(MCNeutrino& other)                                      *
 *                                                                            *
 * Overload = non-const MCNeutrino                                                 *
 *                                                                            *
 * Input:  MCNeutrino                                                              *
 * Output: MCNeutrino                                                              *
 ******************************************************************************/
MCNeutrino& MCNeutrino::operator=(MCNeutrino& other)
{
  MCParticle::operator=(other);
  return *this;
} //MCNeutrino()

/******************************************************************************
 * void MCNeutrino::Fill(TruthTree *truthTree, Int_t iE)                           *
 *                                                                            *
 * Fill MCNeutrino vector from tree                                                *
 *                                                                            *
 * Input:  Event Tree                                                         *
 * Output: None                                                               *
 ******************************************************************************/
void MCNeutrino::FillNeutrino(TruthTree *trtr, int iE)
{
  /*
  //  MCParticle::Fill(trtr, iE);

  Double_t pPt     = trtr -> Neutrino_pt    -> operator[](iE)/1000.;
  Double_t pEta    = trtr -> Neutrino_eta    -> operator[](iE);
  Double_t pPhi    = trtr -> Neutrino_phi    -> operator[](iE);
  Double_t pM      = trtr -> Neutrino_m      -> operator[](iE)/1000.;
  Double_t pCharge = trtr -> Neutrino_charge -> operator[](iE);
  
  _PdgId   = static_cast<Int_t>(trtr -> Neutrino_pdgId  -> operator[](iE));
  _Status  = static_cast<Int_t>(trtr -> Neutrino_status -> operator[](iE));
  _BarCode = trtr -> Neutrino_barcode->operator[](iE);

  */

} //Fill()
