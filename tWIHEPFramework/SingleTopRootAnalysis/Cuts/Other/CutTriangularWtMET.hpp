/******************************************************************************
 * CutTriangularWtMET.hpp                                                    *
 *                                                                            *
 * Triangular cut related to deltaphi(Wtton,MET)                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTriangularWtMET class                       *
 *    CutTriangularWtMET()             -- Parameterized Constructor          *
 *    ~CutTriangularWtMET()            -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutTriangularWtMET"      *
 *                                                                            *
 * Private Data Members of CutTriangularWtMET                                *
 *    myTH1F* _hTriangularWtMETBefore         -- Hist before cuts            *
 *    myTH1F* _hTriangularWtMETAfter          -- Hist after cuts             *
 *    Int_t _TriangularWtMETMin;              -- Minimum requirement         *
 *                                                                            *
 * History                                                                    *
 *      02 Mar 2007 - Created by J. Holzbauer, based on CutJetN                 *
 *****************************************************************************/

#ifndef CutTriangularWtMET_h
#define CutTriangularWtMET_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutTriangularWtMET : public HistoCut 
{

public:

  // Parameterized Constructor
  CutTriangularWtMET(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutTriangularWtMET();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutTriangularWtMET"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Histograms
  myTH1F* _hTriangularWtMETBefore;    // Histogram before the cut
  myTH1F* _hTriangularWtMETAfter;     // Histogram after the cut

  myTH2F* _hTriangularWtMET2DBefore;   // 2D Histogram before the cut
  myTH2F* _hTriangularWtMET2DAfter;    // 2D Histogram after the cut

  // Cut parameters
  Double_t _TriangularWtMETMin;           // Minimum condition
  Double_t _WTransverseMass;

};


#endif







