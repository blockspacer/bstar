/******************************************************************************
 * CutBkgdTriangularWtMET.hpp                                                    *
 *                                                                            *
 * BkgdTriangular cut related to deltaphi(Wtton,MET)                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutBkgdTriangularWtMET class                       *
 *    CutBkgdTriangularWtMET()             -- Parameterized Constructor          *
 *    ~CutBkgdTriangularWtMET()            -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutBkgdTriangularWtMET"      *
 *                                                                            *
 * Private Data Members of CutBkgdTriangularWtMET                                *
 *    myTH1F* _hBkgdTriangularWtMETBefore         -- Hist before cuts            *
 *    myTH1F* _hBkgdTriangularWtMETAfter          -- Hist after cuts             *
 *    Int_t _BkgdTriangularWtMETMin;              -- Minimum requirement         *
 *                                                                            *
 * History                                                                    *
 *      02 Mar 2007 - Created by J. Holzbauer, based on CutJetN                 *
 *****************************************************************************/

#ifndef CutBkgdTriangularWtMET_h
#define CutBkgdTriangularWtMET_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutBkgdTriangularWtMET : public HistoCut 
{

public:

  // Parameterized Constructor
  CutBkgdTriangularWtMET(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutBkgdTriangularWtMET();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutBkgdTriangularWtMET"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hBkgdTriangularWtMETBefore;    // Histogram before the cut
  myTH1F* _hBkgdTriangularWtMETAfter;     // Histogram after the cut

  myTH2F* _hBkgdTriangularWtMET2DBefore;   // 2D Histogram before the cut
  myTH2F* _hBkgdTriangularWtMET2DAfter;    // 2D Histogram after the cut

  // Cut parameters
  Double_t _BkgdTriangularWtMETMin;           // Minimum condition
  Double_t _WTransverseMass;
  Double_t _MET;
};


#endif







