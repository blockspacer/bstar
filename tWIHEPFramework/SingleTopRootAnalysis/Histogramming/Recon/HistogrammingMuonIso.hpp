/**********************************************************************************
 * HistogrammingMuonIso.hpp                                                       *
 *                                                                                *
 * Books and fills histograms                                                     *
 * Books histograms related to isolation                                          *
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
 *      14 Apr 2008 - Created by P. Ryan                                          *
 **********************************************************************************/

#ifndef HistogrammingMuonIso_h
#define HistogrammingMuonIso_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include <sstream>

class HistogrammingMuonIso : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMuonIso(EventContainer *obj, TString muonTypePassed);
  
  // Destructor
  ~HistogrammingMuonIso();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMuonIso"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:

  //////////////////////////////////////
  // All regions
  //////////////////////////////////////
  myTH1F* _hNObj; 
  myTH1F* _hEtCone; 
  myTH1F* _hDeltaR; 
  myTH1F* _hDeltaRFine; 
  myTH1F* _hEta; 
  myTH1F* _hEtaAbs; 


  //////////////////////////////////////
  // Histogram and Cut Vectors
  //////////////////////////////////////
  std::vector<myTH1F*> EtConeHistVector;
  std::vector<myTH1F*> DeltaRHistVector;
  std::vector<myTH1F*> DeltaRFineHistVector;
  std::vector<myTH1F*> EtaAbsHistVector;

  std::vector<Double_t> EtConeCutsVector;
  std::vector<Double_t> DeltaRCutsVector;
  std::vector<Double_t> DeltaRFineCutsVector;
  std::vector<Double_t> EtaAbsCutsVector;

  TString muonType;  // Type of muon being plotted
  
};


#endif






