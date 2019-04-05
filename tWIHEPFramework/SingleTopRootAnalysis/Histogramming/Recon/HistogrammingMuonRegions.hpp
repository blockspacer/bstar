/**********************************************************************************
 * HistogrammingMuonRegions.hpp                                                   *
 *                                                                                *
 * Books and fills histograms                                                     *
 * Books histograms in different pT, eta, and phi regions                         *
 * Only makes hists for Number of Events - used for Efficiency Calculations       *
 *                                                                                *
 * Derived from HistoCut which is in turn derived from BaseCut                    *
 *                                                                                *
 * Public Member Functions of AnalysisMain class                                  *
 *    HistogrammingMuonRegions()    -- Parameterized Constructor                  *
 *    ~HistogrammingMuonRegions()   -- Destructor                                 *
 *    BookHistogramming()                   -- Book histograms                    *
 *    Apply()                               -- Fill histograms only (No Cuts)     *
 *    GetCutName()                 -- Returns "HistogrammingMuonRegions"          *
 *                                                                                *
 * Private Data Members of HistogrammingMuonRegions class                         *
 * - lots of histograms of muon propreties                                        *
 *                                                                                *
 * History                                                                        *
 *      10 Dec 2007 - Created by P. Ryan                                          *
 *      22 Jan 2007 - Updated for new object definitions                          *
 *      2  Dec 2008 - P. Ryan - added 2D rejection histograms                     *
 **********************************************************************************/

#ifndef HistogrammingMuonRegions_h
#define HistogrammingMuonRegions_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include <sstream>

class HistogrammingMuonRegions : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMuonRegions(EventContainer *obj, TString muonTypePassed);
  
  // Destructor
  ~HistogrammingMuonRegions();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMuonRegions"; };

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
  std::vector<myTH1F*> ptEtaRejectionHistVector;
  std::vector<myTH1F*> ptPhiRejectionHistVector;
  std::vector<myTH1F*> etaPhiRejectionHistVector;

  std::vector< std::vector<myTH1F*> > ptEtaHistVectorVector;
  std::vector< std::vector<myTH1F*> > ptPhiHistVectorVector;
  std::vector< std::vector<myTH1F*> > etaPhiHistVectorVector;
  std::vector< std::vector<myTH1F*> > ptEtaHistRejectionVectorVector;
  std::vector< std::vector<myTH1F*> > ptPhiHistRejectionVectorVector;
  std::vector< std::vector<myTH1F*> > etaPhiHistRejectionVectorVector;

  std::vector<Double_t> ptCutsVector;
  std::vector<Double_t> etaCutsVector;
  std::vector<Double_t> phiCutsVector;

  TString muonType;  // Type of muon being plotted
  
};


#endif






