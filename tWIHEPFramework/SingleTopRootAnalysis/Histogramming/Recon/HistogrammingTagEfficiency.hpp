/******************************************************************************
 * HistogrammingTagEfficiency.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTagEfficiency()              -- Parameterized Constructor            *
 *    ~HistogrammingTagEfficiency()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTagEfficiency"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - histograms of MET                                                       *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      20 Nov 2006 - Modified by Bernard Pope                                *
 *      21 Mar 2007 - RS: Fill from event container, add sumET, mex, mey      *
 *****************************************************************************/

#ifndef HistogrammingTagEfficiency_h
#define HistogrammingTagEfficiency_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingTagEfficiency : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingTagEfficiency(EventContainer *obj);
  
  // Destructor
  ~HistogrammingTagEfficiency();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingTagEfficiency"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 

  Float_t _bTagCut;

  //For jets originating from a b quark
  myTH2F* _hTagEff_bFlavNum;  // b jets that pass the b-tagging requirememt
  myTH2F* _hTagEff_bFlavDen;  // All b jets
  
  //For jets originating from a c quark
  myTH2F* _hTagEff_cFlavNum;  // c jets that pass the b-tagging requirememt
  myTH2F* _hTagEff_cFlavDen;  // All c jets

  //For jets originating from a light quark
  myTH2F* _hTagEff_lightFlavNum;  // light jets that pass the b-tagging requirememt
  myTH2F* _hTagEff_lightFlavDen;  // All light jets
  
};


#endif







