/******************************************************************************
 * HistogrammingTopology.hpp                                                  *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTopology()              -- Parameterized Constructor       *
 *    ~HistogrammingTopology()             -- Destructor                      *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTopology"    *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - lots of histograms of topological variables                             *
 *                                                                            *
 * History                                                                    *
 *      21 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *****************************************************************************/

#ifndef HistogrammingTopology_h
#define HistogrammingTopology_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingTopology : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingTopology(EventContainer *obj);
  
  // Destructor
  ~HistogrammingTopology();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingTopology"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
  // Save Histograms
  void SaveHistogram();

private:
  // Histograms declarations 
  myTH1F* _hNJet;         // Histogram of number of jets
  myTH1F* _hNUnTaggedJet; // Histogram of number of untagged jets
  myTH1F* _hNTaggedJet;   // Histogram of number of tagged jets
  myTH1F* _hNEle;         // Histogram of number of electrons
  myTH1F* _hNMu;          // Histogram of number of muons
  myTH1F* _hNTau;         // Histogram of number of taus


  // in order to write out a skim file:
  TString _skimFileName;   // Name of the output skim file
  TFile* _skimFile;        // Output skim file
  TTree* _skimTree;        // Output skim tree

  //
  // and all of the variables that go in the skim tree
  Double_t _eventWeight;
  Int_t _runNumber;
  Int_t _eventNumber;
  Int_t _sourceNumber;
  
  Int_t _nJet;
  Int_t _nUnTaggedJet;
  Int_t _nTaggedJet;
  Int_t _nEle;
  Int_t _nMu;
  Int_t _nTau;

};


#endif







