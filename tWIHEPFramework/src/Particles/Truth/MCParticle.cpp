//////////////////////////////////////////////////////////////////
//
// file: MCParticle.cpp
// Description: Create a container for an MCParticle. 
//
// Author: Reinhard Schwienhorst
// Created: 11/15/2006
//
// Modifications:
//    19 Jan 2007 - Modified by R. Schwienhorst: Adapt for change in truth    *
//                  tree used for MC truth information to "Tru_".             *
//
//
//////////////////////////////////////////////////////////////////
#include "SingleTopRootAnalysis/Particles/Truth/MCParticle.hpp"
#include <iostream>
#include <string>
#include <fstream>

 
using namespace std;

// needed to run with compiled Root:
ClassImp(MCParticle)

//-------------------------------------------------------------------------
// default constructor 
// 
  MCParticle::MCParticle() : Particle::Particle(),_PdgId(0),
					       _Status(0),_BarCode(0)
{
}

//-------------------------------------------------------------------------
// default destructor 
// 
MCParticle::~MCParticle()
{
  //cout<<"In MCParticle destructor"<<endl;
}

//-------------------------------------------------------------------------
// constructor from other MCParticle 
// 
MCParticle::MCParticle(const MCParticle& other): Particle(other),
						 _PdgId(other.PdgId()),_BarCode(other.BarCode()),
						 _Status(other.Status())
{
}

//-------------------------------------------------------------------------
// constructor from other particle
// 
MCParticle::MCParticle(const Particle& other): Particle(other),_PdgId(0),
					       _Status(0),_BarCode(0)
{
}

//-------------------------------------------------------------------------
// operator overloading
// 
// +=
MCParticle& MCParticle::operator+=(const MCParticle& other)
{
  
  Particle::operator+=(other);
  // don't change pdgid or status or stable flags.

  return *this;
}

//-------------------------------------------------------------------------
// 
// addition operator overloading
MCParticle MCParticle::operator+(const MCParticle& other)
{
  //cout<<"beginning of  MCParticle::operator+"<<endl;
  Particle tltemp=Particle::operator+(other);
  MCParticle ptemp=tltemp;

  // I don't know how to add PDG IDs or flags, so for now simply don't modify
  // anything from the first one.
  ptemp.SetPdgId( _PdgId );
  ptemp.SetStatus( _Status );
  ptemp.SetBarCode( _BarCode );


  // done! return the temp particle
  //cout<<"returning from  MCParticle::operator+"<<endl;
  return ptemp;
}

//-------------------------------------------------------------------------
// assignment operator overloading from another MCParticle
// 
MCParticle& MCParticle::operator=(const Particle& other)
{
  
  Particle::operator=(other);
  SetPdgId(0);
  SetStatus(0);
  SetBarCode(0);
  return *this;
}

//-------------------------------------------------------------------------
// assignment operator overloading from another MCParticle
// (const version)
MCParticle& MCParticle::operator=(const MCParticle& other)
{
  
  Particle::operator=(other);
  SetPdgId(other.PdgId());
  SetStatus(other.Status());
  SetBarCode(other.BarCode());
  return *this;
}

//-------------------------------------------------------------------------
// assignment operator overloading from another MCParticle
// (non-const version)
MCParticle& MCParticle::operator=(MCParticle& other)
{
  
  Particle::operator=(other);
  SetPdgId(other.PdgId());
  SetStatus(other.Status());
  SetBarCode(other.BarCode());
  return *this;
}

//-------------------------------------------------------------------------
// initialize this MCParticle from the truth tree
// 
void MCParticle::Fill(TruthTree *trtr,int iE)
{

  // Not in current version
  /*
  Double_t pPt  = trtr->Tru_p_T->operator[](iE)/1000.;
  Double_t pEta = trtr->Tru_eta->operator[](iE);
  Double_t pPhi = trtr->Tru_phi->operator[](iE);
  Double_t pM   = trtr->Tru_m->operator[](iE)/1000.;
  Double_t pCharge= trtr->Tru_charge->operator[](iE);
  
  _PdgId = static_cast<Int_t>(trtr->Tru_pdgId->operator[](iE));
  _Status= static_cast<Int_t>(trtr->Tru_status->operator[](iE));
  _BarCode= trtr->Tru_barcode->operator[](iE);
  */

  
  Double_t pPt     = 0.0;
  Double_t pEta    = 0.0;
  Double_t pPhi    = 0.0;
  Double_t pM      = 0.0;
  Double_t pCharge = 0.0;

  _PdgId = 0;
  _Status= 0;
  _BarCode= 0;


    SetPtEtaPhiM(pPt,pEta,pPhi,pM);
    SetCharge(pCharge);

}



