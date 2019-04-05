/******************************************************************************
 * HistogrammingMCParticle.hpp                                                *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMCParticle()              -- Parameterized Constructor     *
 *    ~HistogrammingMCParticle()             -- Destructor                    *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMCParticle"  *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of MC particle properties                            *
 *                                                                            *
 * History                                                                    *
 *      15 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *****************************************************************************/

#ifndef HistogrammingMCParticle_h
#define HistogrammingMCParticle_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMCParticle : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMCParticle(EventContainer *obj);
  
  // Destructor
  ~HistogrammingMCParticle();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMCParticle"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of all MC particles

  myTH1F *_hPt;       // Pt of all MC particles
  myTH1F *_hEta;      // Eta of all MC particles
  myTH1F *_hPhi;      // Phi of all MC particles
  myTH1F *_hCharge;   // Charge of all MC particles 
  myTH1F *_hPdgId;    // PDG ID of all MC particles 

  // and some histograms for specific particles:
  myTH1F* _hTopN;     // number of top quarks
  myTH1F* _hTopM;     // mass of the top quark
  myTH1F* _hTopPt;    // pt of the top quark
  myTH1F* _hTopEta;   // eta of the top quark
  myTH1F* _hTopPhi;   // phi of the top quark
  myTH1F* _hTopCharge; // charge of the top quark


};


#endif







