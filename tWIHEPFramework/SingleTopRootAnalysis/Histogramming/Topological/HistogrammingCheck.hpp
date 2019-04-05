/******************************************************************************
 * HistogrammingCheck.hpp                                                  *
 *                                                                            *
 * Fills a tree with specific variables.                                      *
 * Should be called after preselection cuts in bin/analysis                   *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingCheck()              -- Parameterized Constructor       *
 *    ~HistogrammingCheck()             -- Destructor                      *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingCheck"    *
 *                                                                            *
 * Private Data Members of this class                                         *
 *        TFile *_skimFile                                                    *
 *        TTree *_skimTree;                                                   *
 *         Lots of variables, histograms                                       *
 *                                                                            *
 * History                                                                    *
 *      21 Nov 2009 - Created by J. Holzbauer for overlap, duplicate event checks*
 *****************************************************************************/

#ifndef HistogrammingCheck_h
#define HistogrammingCheck_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;


class HistogrammingCheck : public HistoCut 
{

public:

  HistogrammingCheck(EventContainer *obj,TFile *skimFilez);
  
  // destructor
  ~HistogrammingCheck();

  // methods:
  void BookHistogram();
  bool Apply();
  
  inline std::string GetCutName (void) { return "HistogrammingCheck"; };

private:
  TFile *_skimFile;
  TTree *_skimTree;
   
  //------------------------------------------------------
  // Event Weight
  //------------------------------------------------------
  Float_t _EventWeight;
  Int_t _RunNumber;
  Int_t _EventNumber;
  
  Float_t _Jet1Pt;
  Float_t _Jet2Pt;
  Float_t _Jet3Pt;
  Float_t _Jet4Pt;
  Float_t _Jet1Eta;
  Float_t _Jet2Eta;
  Float_t _Jet3Eta;
  Float_t _Jet4Eta;
  Float_t _DeltaRJet1Jet2;
  Float_t _DeltaRJet1Electron1;
  Float_t _DeltaRJet1Muon1;
  Float_t _MET;
  

  ////////////////////////////////////////////////////////
  // Histograms declaration 
  ////////////////////////////////////////////////////////

  myTH1F* _hEventWeight;
  myTH1F* _hRunNumber;
  myTH1F* _hEventNumber;
  myTH1F* _hAbsLogEventWeight;
  
  myTH1F* _hJet1Pt;
  myTH1F* _hJet2Pt;
  myTH1F* _hJet3Pt;
  myTH1F* _hJet4Pt;
  myTH1F* _hJet1Eta;
  myTH1F* _hJet2Eta;
  myTH1F* _hJet3Eta;
  myTH1F* _hJet4Eta;
  myTH1F* _hDeltaRJet1Jet2;
  myTH1F* _hDeltaRJet1Electron1;
  myTH1F* _hDeltaRJet1Muon1;
  myTH1F* _hMET;
  

};


#endif


