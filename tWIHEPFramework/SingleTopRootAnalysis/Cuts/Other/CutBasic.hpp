/******************************************************************************
 * CutJetN.hpp                                                               *
 *                                                                            *
 * Apply a default set of cuts to all classes deriving from BaseCuts          *
 *                                                                            *
 * Derived from BaseCut                                                       *
 *                                                                            *
 * Public Member Functions of CutBasic class                                  *
 *    CutBasic()                        -- Parameterized Constructor          *
 *    ~CutBasic()                       -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutBasic"                 *
 *                                                                            *
 * Private Data Members of CutBasic                                           *
 *    Double_t _minElePt                                                      *  
 *    Double_t _maxEleEta                                                     *
 *    Double_t _minMET                                                        *
 *    Double_t _minJetPt                                                      *
 *                                                                            *
 * History                                                                    *
 *      28 Nov 2004 - Created by R. Schwienhorst                              *
 *      10 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#ifndef CutBasic_h
#define CutBasic_h

#include "SingleTopRootAnalysis/Base/Dictionary/BaseCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

class CutBasic : public BaseCut
{

public:

  // Parameterized constructor
  CutBasic(EventContainer *obj, Double_t elept = -1.0, Double_t eleeta = 99.0, Double_t met = -1.0,
	   Double_t jetpt = -1.0, Double_t jeteta = 99.0, 
	   Int_t minjets = 0, Int_t maxjets = 99, 
	   Double_t jet1pt = 0.0, Double_t jet1eta = 99.0);
  
  // Destructor
  ~CutBasic();

  // Book histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Return name of Cut
  inline std::string GetCutName() { return "CutBasic"; };


private:

  // Cut parameters
  Double_t _minElePt;   // min electron pt
  Double_t _maxEleEta;  // min electron eta
  Double_t _minMET;     // min missing ET
  Double_t _minJetPt;   // min Jet pT
  Double_t _maxJetEta;  // max Jet eta

  Int_t _minNjets;      // max number of jets
  Int_t _maxNjets;      // min number of jets

  Double_t _minJet1Pt;  // leading jet pt
  Double_t _maxJet1Eta; // leading jet eta

};


#endif







