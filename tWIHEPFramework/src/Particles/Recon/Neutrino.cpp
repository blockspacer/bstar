/******************************************************************************
 * Neutrino.cpp                                                                    *
 *                                                                            *
 * Store information about final state Neutrinos                                   *
 *                                                                            *
 * Derived from Particle class                                                *
 *                                                                            *
 * Public Member Functions of Neutrino class                                       *
 *    Neutrino()                              -- Default Constructor               *
 *    Neutrino()                              -- Parameterized Constructor         *
 *    Neutrino()                              -- Constructor with Particle         *
 *    ~Neutrino()                             -- Destructor                        *
 *    Clear()                            -- Set contents to default           *
 *    Fill                               -- Fill neutrino from event tree          *
 *    FillFastSim                        -- Fill neutrino from FastSim tree        *
 *    operator+=                         -- Overload +=                       *
 *    operator+                          -- Overload +                        *
 *    operator=                          -- Overload = Particle               *
 *    operator=                          -- Overload = const Neutrino              *
 *    operator=                          -- Overload =                        *
 *                                                                            *
 * Private Data Members of Neutrino class                                          *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst                              *
 *      20 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      20 Mar 2007 - RS: Add filling from fastsim tree                       *
 *      11 Jan 2008 - RS: Fill only if object passes ID cuts                  *
 *                        Update neutrino filling method.                     *
 *****************************************************************************/
#include "SingleTopRootAnalysis/Particles/Recon/Neutrino.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include "SingleTopRootAnalysis/Base/Histograms/utils.hpp"


using namespace std;

// Integrate classes into the Root system
ClassImp(Neutrino)
/******************************************************************************
 * Neutrino::Neutrino()                                                       *
 *                                                                            *
 * Default Constructor                                                        *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
Neutrino::Neutrino() : Particle::Particle()
{

} //Neutrino()

/******************************************************************************
 * Neutrino::~Neutrino()                                                              *
 *                                                                            *
 * Destructor                                                                 *  
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
Neutrino::~Neutrino()
{

} //~Neutrino()


/******************************************************************************
 * Neutrino::Neutrino(const Neutrino& other)                                              *
 *                                                                            *
 * Copy constructor Neutrino const                                                *
 *                                                                            *
 * Input:  Neutrino                                                               *
 * Output: None                                                               *
 ******************************************************************************/
Neutrino::Neutrino(const Neutrino& other): Particle(other)
{

} //Neutrino()

/******************************************************************************
 * Neutrino::Neutrino(const Neutrino& other)                                                 *
 *                                                                            *
 * Copy constructor Particle const                                            *
 *                                                                            *
 * Input:  Neutrino                                                                *
 * Output: None                                                               *
 ******************************************************************************/
Neutrino::Neutrino(const Particle& other): Particle(other)
{

} //Neutrino() particle

/******************************************************************************
 * Neutrino& Neutrino::operator+=(const Neutrino& other)                                     *
 *                                                                            *
 * Overload +=                                                                *
 *                                                                            *
 * Input:  Neutrino                                                                *
 * Output: Neutrino                                                                *
 ******************************************************************************/
Neutrino& Neutrino::operator+=(const Neutrino& other)
{
  Particle::operator+=(other);
  return *this;
} //+=

/******************************************************************************
 * Neutrino& Neutrino::operator+(const Neutrino& other)                                      *
 *                                                                            *
 * Overload +                                                                 *
 *                                                                            *
 * Input:  Neutrino                                                                *
 * Output: Neutrino                                                                *
 ******************************************************************************/
Neutrino Neutrino::operator+(const Neutrino& other)
{

  Particle particleTemp = Particle::operator+(other);
  Neutrino neutrinoTemp = particleTemp;
  
   // Return the temp particle
  return neutrinoTemp;
} //+

/******************************************************************************
 * Neutrino& Neutrino::operator=(const Particle& other)                                 *
 *                                                                            *
 * Overload = Particle                                                        *
 *                                                                            *
 * Input:  Particle                                                           *
 * Output: Neutrino                                                                *
 ******************************************************************************/
Neutrino& Neutrino::operator=(const Particle& other)
{
  Particle::operator=(other);

  return *this;
} // = Particle

/******************************************************************************
 * Neutrino& Neutrino::operator=(const Neutrino& other)                                      *
 *                                                                            *
 * Overload = const Neutrino                                                       *
 *                                                                            *
 * Input:  Neutrino                                                                *
 * Output: Neutrino                                                                *
 ******************************************************************************/
Neutrino& Neutrino::operator=(const Neutrino& other)
{
  
  Particle::operator=(other);
 
  return *this;
} //= const neutrino

/******************************************************************************
 * Neutrino& Neutrino::operator=(const Neutrino& other)                                      *
 *                                                                            *
 * Overload = non-const Neutrino                                                   *
 *                                                                            *
 * Input:  Neutrino                                                                *
 * Output: Neutrino                                                                *
 ******************************************************************************/
Neutrino& Neutrino::operator=(Neutrino& other)
{
  
  Particle::operator=(other);
 
  return *this;
} //= neutrino non-const

/******************************************************************************         
 * void Neutrino::Fill(std::vector<Muon>& muons,std::vector<Electron>& electrons,
 *                     const Double_t& missingEx, const Double_t& missingEy,) *         
 *                                                                            *         
 * Fill Neutrino vector from tree                                             *         
 *                                                                            *         
 * Input:  electrons and muons and mex and mey                                *         
 * Output: kTRUE                                                              *         
 ******************************************************************************/
Bool_t Neutrino::Fill(const vector<Muon>& tightMuons,const vector<Electron>& tightElectrons,const Double_t& missingEx,const Double_t& missingEy)
{

  Double_t leadlep_px = 0;
  Double_t leadlep_py = 0;
  Double_t leadlep_pz = 0;
  Double_t leadlep_e = 0;
  
  //choose which is the leading high pt lepton
  if((tightMuons.size() > 0)&&(tightElectrons.size() == 0)){
    leadlep_px = tightMuons[0].Px();
    leadlep_py = tightMuons[0].Py();
    leadlep_pz = tightMuons[0].Pz();
    leadlep_e = tightMuons[0].E();
  }
  else if((tightMuons.size()== 0)&&(tightElectrons.size() > 0)){
    leadlep_px = tightElectrons[0].Px();
    leadlep_py = tightElectrons[0].Py();
    leadlep_pz = tightElectrons[0].Pz();
    leadlep_e = tightElectrons[0].E()/1000;
  }
  else if((tightMuons.size() > 0)&&(tightElectrons.size() > 0)){
    if(tightMuons[0].Pt() > tightElectrons[0].Pt()){
      leadlep_px = tightMuons[0].Px();
      leadlep_py = tightMuons[0].Py();
      leadlep_pz = tightMuons[0].Pz();
      leadlep_e = tightMuons[0].E();
    }else{
      leadlep_px = tightElectrons[0].Px();
      leadlep_py = tightElectrons[0].Py();
      leadlep_pz = tightElectrons[0].Pz();
      leadlep_e = tightElectrons[0].E();
    }
  }
    
  // last argument to neutrino_pz is a large negative number
  // to indicate that we don't want to use a scale factor.
  
  //I send in one for the true neutrino pz, at least for now.  not used in calculation
  Double_t scf=-99999.;
  //cout<<"calling pzneutrino"<<endl;
  Double_t nu_pz = -1;
  Double_t nu_e = 0;

  //  cout<<"leadlep "<<leadlep_e<<endl;
  if(( leadlep_px == 0)&&( leadlep_py == 0)&&( leadlep_pz == 0)&&( leadlep_e == 0)){
    //no leading lepton (or is tau)
    nu_pz = 0.; 
  }
  else{
    nu_pz=Neutrino_Pz(missingEx, missingEy,leadlep_px,leadlep_py,leadlep_pz,leadlep_e,scf,1.0);
     //   cout<<"pzpypxE "<<nu_pz<<"  "<<nu_py<<"  "<<nu_px<<"  "<<nu_e<<endl;
  }

  nu_e=sqrt(missingEx*missingEx + missingEy*missingEy + nu_pz*nu_pz);
    
  SetPxPyPzE(missingEx,missingEy,nu_pz,nu_e);

//   finalNeutrino.SetPdgId(-12);
  
  return kTRUE;
} //Fill()

