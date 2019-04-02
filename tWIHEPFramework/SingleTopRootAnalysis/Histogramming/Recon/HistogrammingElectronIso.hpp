/**********************************************************************************
 * HistogrammingElectronIso.hpp                                                   *
 *                                                                                *
 * Books and fills histograms                                                     *
 * Books histograms related to isolation                                          *
 *                                                                                *
 * Derived from HistoCut which is in turn derived from BaseCut                    *
 *                                                                                *
 * Public Member Functions of AnalysisMain class                                  *
 *    HistogrammingElectronRegions()    -- Parameterized Constructor              *
 *    ~HistogrammingElectronRegions()   -- Destructor                             *
 *    BookHistogramming()                   -- Book histograms                    *
 *    Apply()                               -- Fill histograms only (No Cuts)     *
 *    GetCutName()                 -- Returns "HistogrammingElectronRegions"      *
 *                                                                                *
 * Private Data Members of HistogrammingElectronRegions class                     *
 * - lots of histograms of muon propreties                                        *
 *                                                                                *
 * History                                                                        *
 *      14 Apr 2008 - Created by P. Ryan                                          *
 **********************************************************************************/

#ifndef HistogrammingElectronIso_h
#define HistogrammingElectronIso_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include <sstream>

class HistogrammingElectronIso : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingElectronIso(EventContainer *obj, TString electronTypePassed);
  
  // Destructor
  ~HistogrammingElectronIso();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingElectronIso"; };

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

  TString electronType;  // Type of electron being plotted
  
};


#endif






