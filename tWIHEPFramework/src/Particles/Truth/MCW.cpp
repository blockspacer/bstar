/******************************************************************************
 * MCW.cpp                                                                    *
 *                                                                            *
 * Container for a MCW                                                        *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of MCW class                                       *
 *    MCW()                           -- Default Constructor                  *
 *    MCW()                           -- Parameterized Constructor            *
 *    MCW()                           -- Constructor with Particle            *
 *    ~MCW()                          -- Destructor                           *
 *    Clear()                            -- Set contents to default           *
 *    Fill()                             -- Fill w from event tree            *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const MCW              *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of MCW class                                          *
 *                                                                            *
 * History                                                                    *
 *      15 May 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Particles/Truth/MCW.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// Integrate classes into the Root system
ClassImp(MCW)

/******************************************************************************
 * MCW::MCW()                                                                 *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
  MCW::MCW()
{

} //MCW()

/******************************************************************************
 * MCW::~MCW()                                                                *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
MCW::~MCW()
{

} //~MCW()


/******************************************************************************
 * MCW::MCW(const MCW& other)                                                 *
 *                                                                            *
 * Copy constructor MCW const                                                 *
 *                                                                            *
 * Input:  MCW                                                                *
 * Output: None                                                               *
 ******************************************************************************/
MCW::MCW(const MCW& other): MCParticle(other)
{

} //MCW()

/******************************************************************************
 * MCW::MCW(const MCW& other)                                                 *
 *                                                                            *
 * Copy constructor MCParticle const                                          *
 *                                                                            *
 * Input:  MCW                                                                *
 * Output: None                                                               *
 ******************************************************************************/
MCW::MCW(const MCParticle& other): MCParticle(other)
{ 

} //MCW

/******************************************************************************
 * MCW& MCW::operator+=(const MCW& other)                                     *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  MCW                                                                *
 * Output: MCW                                                                *
 ******************************************************************************/
MCW& MCW::operator+=(const MCW& other)
{
  
  MCParticle::operator+=(other);
   return *this;
} //+=

/******************************************************************************
 * MCW& MCW::operator+(const MCW& other)                                      *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  MCW                                                                *
 * Output: MCW                                                                *
 ******************************************************************************/
MCW MCW::operator+(const MCW& other)
{

  MCParticle::operator+(other);
  return *this;
} //+

/******************************************************************************
 * MCW& MCW::operator=(const W& other)                                        *
 *                                                                            *
 * Overload = MCParticle                                                      *
 *                                                                            *
 * Input:  MCParticleicle::Clear();                                           *
 * Output: MCW                                                                *
 ******************************************************************************/
MCW& MCW::operator=(const MCParticle& other)
{
  
  MCParticle::operator=(other);
  return *this;
} // MCParticle

/******************************************************************************
 * MCW& MCW::operator=(const MCW& other)                                      *
 *                                                                            *
 * Overload = const MCW                                                       *
 *                                                                            *
 * Input:  MCW                                                                *
 * Output: MCW                                                                *
 ******************************************************************************/
MCW& MCW::operator=(const MCW& other)
{
  MCParticle::operator=(other);
  return *this;
} //MCW

/******************************************************************************
 * MCW& MCW::operator=(MCW& other)                                            *
 *                                                                            *
 * Overload = non-const MCW                                                   *
 *                                                                            *
 * Input:  MCW                                                                *
 * Output: MCW                                                                *
 ******************************************************************************/
MCW& MCW::operator=(MCW& other)
{
  MCParticle::operator=(other);
  return *this;
} //MCW()

/******************************************************************************
 * void MCW::Fill(TruthTree *truthTree, Int_t iE)                             *
 *                                                                            *
 * Fill MCW vector from tree                                                  *
 *                                                                            *
 * Input:  Event Tree                                                         *
 * Output: None                                                               *
 ******************************************************************************/
void MCW::FillW(TruthTree *trtr, int iE)
{

  MCParticle::Fill(trtr, iE);

} //Fill()
