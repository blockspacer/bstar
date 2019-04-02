/********************************************************************************
 * HistogrammingResolutionJet.hpp                                               *
 *                                                                              *
 * Books and fills histograms for jet resolutions (recon - truth)               *
 * Used for events passing cuts applied in other classes                        *
 *                                                                              *
 * Derived from HistoCut which is in turn derived from BaseCut                  *
 *                                                                              *
 * Public Member Functions of AnalysisMain class                                *
 *    HistogrammingResolutionJet()      -- Parameterized Constructor            *
 *    ~HistogrammingResolutionJet()     -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                      *
 *    Apply()                           -- Fill histograms only (No Cuts)       *
 *    GetCutName()                      -- Returns "HistogrammingResolutionJet" *
 *                                                                              *
 * Private Data Members of this class                                           *
 *   - lots of histograms of jet properties                                     *
 *                                                                              *
 * History                                                                      *
 *      26 Mar 2008 - Created by P. Ryan                                        *
 ********************************************************************************/

#ifndef HistogrammingResolutionJet_h
#define HistogrammingResolutionJet_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingResolutionJet : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingResolutionJet(EventContainer *obj);
  
  // Destructor
  ~HistogrammingResolutionJet();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingResolutionJet"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hNObj;     // Histogram of number of jets

  myTH1F* _hPtObj1;   // Jet 1 PT
  myTH1F* _hEtaObj1;  // Jet 1 Eta
  myTH1F* _hPhiObj1;  // Jet 1 Phi

  myTH1F* _hEtaObj2;  // Jet 2 Eta
  myTH1F* _hPtObj2;   // Jet 2 PT
  myTH1F* _hPhiObj2;  // Jet 2 Phi

  myTH1F* _hPtObj3;   // Jet 3 PT
  myTH1F* _hEtaObj3;  // Jet 3 Eta
  myTH1F* _hPhiObj3;  // Jet 3 Phi

  myTH1F* _hPtObj4;   // Jet 4 PT
  myTH1F* _hEtaObj4;  // Jet 4 Eta
  myTH1F* _hPhiObj4;  // Jet 4 Phi

  myTH1F* _hPtObj5;   // Jet 5 PT
  myTH1F* _hEtaObj5;  // Jet 5 Eta
  myTH1F* _hPhiObj5;  // Jet 5 Phi

  myTH1F* _hPtObj6;   // Jet 6 PT
  myTH1F* _hEtaObj6;  // Jet 6 Eta
  myTH1F* _hPhiObj6;  // Jet 6 Phi

  myTH1F* _hPtObj7;   // Jet 7 PT
  myTH1F* _hEtaObj7;  // Jet 7 Eta
  myTH1F* _hPhiObj7;  // Jet 7 Phi

  myTH1F *_hPt;       // Pt of all jets
  myTH1F *_hEta;      // Eta of all jets
  myTH1F *_hPhi;      // Phi of all jets
  myTH1F *_hCharge;   // Charge of all jets 
  myTH1F *_hEtEM0;    // ET in EM layer 0 for all jets



};


#endif







