/**********************************************************************************
 * HistogrammingElectronRegions.hpp                                               *
 *                                                                                *
 * Books and fills histograms                                                     *
 * Books histograms in different pT, eta, and phi regions                         *
 * Only makes hists for Number of Events - used for Efficiency Calculations       *
 *                                                                                *
 * Derived from HistoCut which is in turn derived from BaseCut                    *
 *                                                                                *
 * Public Member Functions of AnalysisMain class                                  *
 *    HistogrammingElectronRegions()  -- Parameterized Constructor                *
 *    ~HistogrammingElectronRegions() -- Destructor                               *
 *    BookHistogramming()             -- Book histograms                          *
 *    Apply()                         -- Fill histograms only (No Cuts)           *
 *    GetCutName()                    -- Returns "HistogrammingElectronRegions"   *
 *                                                                                *
 * Private Data Members of HistogrammingElectronRegions class                     *
 * - lots of histograms of electron propreties                                    *
 *                                                                                *
 * History                                                                        *
 *      23 Jan 2008 - Created by P. Ryan                                          *
 **********************************************************************************/

#ifndef HistogrammingElectronRegions_h
#define HistogrammingElectronRegions_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include <sstream>

class HistogrammingElectronRegions : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingElectronRegions(EventContainer *obj, TString electronTypePassed);
  
  // Destructor
  ~HistogrammingElectronRegions();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingElectronRegions"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:

  //////////////////////////////////////
  // All regions
  //////////////////////////////////////
  myTH1F* _hNObj; 


  //////////////////////////////////////
  // Histogram and Cut Vectors
  //////////////////////////////////////
  std::vector<myTH1F*> ptHistVector;
  std::vector<myTH1F*> ptTurnOnHistVector;
  std::vector<myTH1F*> etaHistVector;
  std::vector<myTH1F*> phiHistVector;
  std::vector<myTH1F*> ptEtaHistVector;
  std::vector<myTH1F*> ptPhiHistVector;
  std::vector<myTH1F*> etaPhiHistVector;
  std::vector< std::vector<myTH1F*> > ptEtaHistVectorVector;
  std::vector< std::vector<myTH1F*> > ptPhiHistVectorVector;
  std::vector< std::vector<myTH1F*> > etaPhiHistVectorVector;

  std::vector<Double_t> ptCutsVector;
  std::vector<Double_t> etaCutsVector;
  std::vector<Double_t> phiCutsVector;

  TString electronType;  // Type of electron being plotted
  
};


#endif






