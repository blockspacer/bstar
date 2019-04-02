//////////////////////////////////////////////////////////////////
//
// file: EventContainer.cpp
// Description: Create TLorentzVectors from particles
//             in the tree. 
//
// Author: Reinhard Schwienhorst
// Created: 10/19/2006
//
// Modifications:
//
//
//////////////////////////////////////////////////////////////////
#include "SingleTopRootAnalysis/Particles/Recon/Particle.hpp"
#include <iostream>
#include <string>
#include <fstream>


using namespace std;

// needed to run with compiled Root:
ClassImp(Particle)

//-------------------------------------------------------------------------
// default constructor 
// 
  Particle::Particle() : TLorentzVector::TLorentzVector(),_charge(0.),_debugLevel(0)
{
}

//-------------------------------------------------------------------------
// default destructor 
// 
Particle::~Particle()
{
  //cout<<"In particle destructor"<<endl;
}

//-------------------------------------------------------------------------
// constructor from other particle 
// 
Particle::Particle(const Particle& other): TLorentzVector(other),_charge(other.Charge()),
					   _debugLevel(other.GetDebugLevel())
{
}

//-------------------------------------------------------------------------
// constructor from other particle, non-const version
// 
Particle::Particle(Particle& other): TLorentzVector(other),_charge(other.Charge()),
					   _debugLevel(other.GetDebugLevel())
{
}

//-------------------------------------------------------------------------
// constructor from other TLorentzVector
// 
Particle::Particle(const TLorentzVector& other): TLorentzVector(other),_charge(0.),
						 _debugLevel(0)
{
}

//-------------------------------------------------------------------------
// operator overloading
// 
// +=
Particle& Particle::operator+=(const Particle& other)
{
  
  TLorentzVector::operator+=(other);
  // also add charges
  SetCharge(Charge()+other.Charge());
  return *this;
}

//-------------------------------------------------------------------------
// 
// addition operator overloading
Particle Particle::operator+(const Particle& other)
{
  //cout<<"beginning of  Particle::operator+"<<endl;
  TLorentzVector tltemp=TLorentzVector::operator+(other);
  Particle ptemp=tltemp;
  
  ptemp.SetCharge(Charge()+other.Charge());

  // done! return the temp particle
  //cout<<"returning from  Particle::operator+"<<endl;
  return ptemp;
}

//-------------------------------------------------------------------------
// assignment operator overloading from a TLorentzVector
// 
Particle& Particle::operator=(TLorentzVector& other)
{
  
  TLorentzVector::operator=(other);
  SetCharge(0.);
  SetDebugLevel(0);
  return *this;
}

//-------------------------------------------------------------------------
// assignment operator overloading from another Particle
// 
Particle& Particle::operator=(const Particle& other)
{
  
  TLorentzVector::operator=(other);
  SetCharge(other.Charge());
  SetDebugLevel(other.GetDebugLevel());
  return *this;
}
