/******************************************************************************
 * HistogrammingJetAngular.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingJetAngular()              -- Parameterized Constructor            *
 *    ~HistogrammingJetAngular()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingJetAngular"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - histograms of MET                                                       *
 *                                                                            *
 * History                                                                    *
 *      12 Jul 2016 - Created by D. Leggat for IHEP                           *
 *****************************************************************************/

#ifndef HistogrammingJetAngular_h
#define HistogrammingJetAngular_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Jet.hpp"
#include "SingleTopRootAnalysis/Particles/Recon/Particle.hpp"

class HistogrammingJetAngular : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingJetAngular(EventContainer *obj, bool unisolated = false);
  
  // Destructor
  ~HistogrammingJetAngular();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingJetAngular"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hUntaggedJetDelR;  // Histogram of delta R between the untagged jets
  myTH1F* _hJet1bJetDelR; // delta R between leading untagged jet and b tagged jet
  myTH1F* _hJet2bJetDelR; // delta R between sub-leading untagged jet and b tagged jet
  myTH1F* _hlepJetsDelR; 
  myTH1F* _hlepbJetDelR; // del R between the b jet and the lepton
  myTH1F* _hJet1lepDelR; // delta R between leading untagged jet and lepton
  myTH1F* _hJet2lepDelR; // delta R between sub-leading untagged jet and lepton

  //Now those same histograms but for delta phi, because why not?
  myTH1F* _hUntaggedJetDelPhi;  // Histogram of delta phi between the untagged jets
  myTH1F* _hJet1bJetDelPhi; // delta phi between leading untagged jet and b tagged jet
  myTH1F* _hJet2bJetDelPhi; // delta phi between sub-leading untagged jet and b tagged jet
  myTH1F* _hlepJetsDelPhi; 
  myTH1F* _hlepbJetDelPhi; // del R between the b jet and the lepton
  myTH1F* _hJet1lepDelPhi; // delta phi between leading untagged jet and lepton
  myTH1F* _hJet2lepDelPhi; // delta phi between sub-leading untagged jet and lepton
  myTH1F* _hJet1METDelPhi; // delta phi between the leading untagged jet and the MET
  myTH1F* _hJet2METDelPhi; // delta phi between the sub-leading untagged jet and the MET
  myTH1F* _hbJetMETDelPhi; // delta phi between the b jet and the MET
  myTH1F* _hJetsMETDelPhi; // delta phi between the untagged jets and the MET
  myTH1F* _hLepMETDelPhi; // delta phi between NET and lepton. Not really a jet angle, but whatevs.

  //Also going to make some combination plots here. Because JetAngular is a terrible name for this class anyway.
  myTH1F* _diJetMass; // mass of the untagged jet pair
  myTH1F* _threeJetMass; // mass of the three jet system
  myTH1F* _mTt; //Transverse mass of the reconstructed top
  myTH1F* _ptSys; // p_T of the system

  bool _unisolated; //Whether we are using tight or unisolated leptons

  Jet bJet;
  Jet leadingJet;
  Jet subleadingJet;
  Particle lepton;
  TLorentzVector met;

};


#endif







