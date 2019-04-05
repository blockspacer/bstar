/******************************************************************************
 * HistogrammingMCNeutrino.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMCNeutrino()              -- Parameterized Constructor            *
 *    ~HistogrammingMCNeutrino()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMCNeutrino"         *
 *                                                                            *
 * Private Data Members of HistogrammingMCNeutrino class                             *
 * - lots of histograms of w propreties                                       *
 *                                                                            *
 * History                                                                    *
 *      15 May 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/

#ifndef HistogrammingMCNeutrino_h
#define HistogrammingMCNeutrino_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMCNeutrino : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMCNeutrino(EventContainer *obj);
  
  // Destructor
  ~HistogrammingMCNeutrino();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMCNeutrino"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of ws
  myTH1F* _hPtObj1;   // W 1 PT
  myTH1F* _hPtObj2;   // W 2 PT
  myTH1F* _hEtaObj1;  // W 1 Eta
  myTH1F* _hEtaObj2;  // W 2 Eta
  myTH1F* _hPhiObj1;  // W 1 Phi
  myTH1F* _hPhiObj2;  // W 2 Phi

  myTH1F *_hPt;       // Pt of all ws
  myTH1F *_hPz;       // Pt of all ws
  myTH1F *_hEta;      // Eta of all ws
  myTH1F *_hPhi;      // Phi of all ws
  myTH1F *_hCharge;   // Charge of all ws 
  myTH1F *_hNeutrinoPdgId;  //PDG ID of all ws
 myTH1F *_hE;   // Charge of all ws 
  myTH1F *_hM;   // Charge of all ws 
};


#endif


