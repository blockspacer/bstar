/******************************************************************************
 * CutTriangularLepMET.hpp                                                    *
 *                                                                            *
 * Triangular cut related to deltaphi(lepton,MET)                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTriangularLepMET class                       *
 *    CutTriangularLepMET()             -- Parameterized Constructor          *
 *    ~CutTriangularLepMET()            -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutTriangularLepMET"      *
 *                                                                            *
 * Private Data Members of CutTriangularLepMET                                *
 *    myTH1F* _hTriangularLepMETBefore         -- Hist before cuts            *
 *    myTH1F* _hTriangularLepMETAfter          -- Hist after cuts             *
 *    Int_t _TriangularLepMETMin;              -- Minimum requirement         *
 *                                                                            *
 * History                                                                    *
 *      02 Mar 2007 - Created by J. Holzbauer, based on CutJetN                 *
 *****************************************************************************/

#ifndef CutTriangularLepMET_h
#define CutTriangularLepMET_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutTriangularLepMET : public HistoCut 
{

public:

  // Parameterized Constructor
  CutTriangularLepMET(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutTriangularLepMET();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutTriangularLepMET"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hTriangularLepMETBefore;    // Histogram before the cut
  myTH1F* _hTriangularLepMETAfter;     // Histogram after the cut

  myTH2F* _hTriangularLepMET2DBefore;   // 2D Histogram before the cut
  myTH2F* _hTriangularLepMET2DAfter;    // 2D Histogram after the cut

  myTH1F* _hTriangularLepMETCutParamBefore;    // Histogram before the cut
  myTH1F* _hTriangularLepMETCutParamAfter;     // Histogram after the cut

  // Cut parameters
  Double_t _TriangularLepMETMin;           // Minimum condition


};


#endif







