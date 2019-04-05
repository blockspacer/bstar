/******************************************************************************
 * HistogrammingNPvtx.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingNPvtx()              -- Parameterized Constructor            *
 *    ~HistogrammingNPvtx()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingNPvtx"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - histograms of MET                                                       *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      20 Nov 2006 - Modified by Bernard Pope                                *
 *      21 Mar 2007 - RS: Fill from event container, add sumET, mex, mey      *
 *****************************************************************************/

#ifndef HistogrammingNPvtx_h
#define HistogrammingNPvtx_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingNPvtx : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingNPvtx(EventContainer *obj);
  
  // Destructor
  ~HistogrammingNPvtx();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingNPvtx"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNPvtx;       // Histogram of number of Primary Vetex
  myTH2F* _hNPvtx_nJet;  // Histogram of NPvtx vs Jet number
  myTH2F* _hNPvtx_nElectron;  // Histogram of NPvtx vs Electron number
  myTH2F* _hNPvtx_nMuon;  // Histogram of NPvtx vs Muon number

  myTH1F* _hNTrueInteractions; //Number of true interactions in the event
  myTH2F* _hNTrueInteractions_nElectron; // Number of true interactions vs n electrons
  myTH2F* _hNTrueInteractions_nMuon; // Number of true interactions vs n muons
  myTH2F* _hNTrueInteractions_nJet; // Number of true interactions vs n jets

  myTH1F* _hNPUVertices; //Number of pileup vertices
  myTH2F* _hNPUVertices_nElectron; // Number of pileup interactions vs n electrons
  myTH2F* _hNPUVertices_nMuon; // Number of pileup interactions vs n muons
  myTH2F* _hNPUVertices_nJet; // Number of pileup interactions vs n jets
  
};


#endif







