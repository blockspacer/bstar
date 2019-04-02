/******************************************************************************
 * HistogrammingNeutrino.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingNeutrino()              -- Parameterized Constructor            *
 *    ~HistogrammingNeutrino()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingNeutrino"         *
 *                                                                            *
 * Private Data Members of HistogrammingNeutrino class                             *
 * - lots of histograms of w propreties                                       *
 *                                                                            *
 * History                                                                    *
 *      15 May 2007 - Created by J. Holzbauer                                 *
 *****************************************************************************/

#ifndef HistogrammingNeutrino_h
#define HistogrammingNeutrino_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingNeutrino : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingNeutrino(EventContainer *obj);
  
  // Destructor
  ~HistogrammingNeutrino();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingNeutrino"; };

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
  myTH1F *_hE;   // Charge of all ws 
  myTH1F *_hM;   // Charge of all ws 
};


#endif


