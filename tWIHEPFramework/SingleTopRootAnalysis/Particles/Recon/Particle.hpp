///////////////////////////////////////////////////////////////////////////
// header file: Particle
// 
// Description: 
//  store information about final state objects.
//
// Modifications:
//
//
////////////////////////////////////////////////////////////////////////////
#ifndef particle_h
#define particle_h

#include <TLorentzVector.h>
#include <TEnv.h>
//

class Particle: public TLorentzVector
{
 public:
  //default constructor
  Particle();
  //default constructor
  Particle(const Particle& other);
  Particle(Particle& other);
  // constructor from TLorentzVector
  Particle(const TLorentzVector& other);
  //destructor
  ~Particle();

  // overloaded operators:
  // +=
  Particle& operator+=(const Particle& other);
  // addition:
  Particle operator+(const Particle& other);
  // overloaded assignment operator:
  Particle& operator=(TLorentzVector& other);
  // overloaded assignment operator:
  Particle& operator=(const Particle& other);

  // clear: set all contents to their defaults
  inline void Clear() { TLorentzVector::Clear(); _charge=0.; };

  // Setters
  // Set the charge of this particle
  void SetCharge(int iCh) { _charge = static_cast<Double_t>(iCh); };
  void SetCharge(const Double_t& ch) { _charge = ch; };

  // Getters
  // Get the charge of this particle
  Double_t GetCharge() const { return _charge; };
  Double_t Charge() const { return _charge; };


  //
  // debugging
  inline void SetDebugLevel(int i) { _debugLevel=i; };
  inline int GetDebugLevel() const { return _debugLevel; };
  inline int debugLevel() const { return _debugLevel; };

 private:
  
  Double_t _charge;  // charge of particle

  int _debugLevel; // debugging
  ////////////////////////////////////////
  // needed to run with Root:
  ClassDef(Particle,0)
};


#endif
