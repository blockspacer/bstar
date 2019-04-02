/******************************************************************************
 * CutTriangularSumDeltaPhiLepMET.hpp                                                    *
 *                                                                            *
 * Triangular cut related to the sum of all deltaphi(lepton,MET)                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutTriangularSumDeltaPhiLepMET class                       *
 *    CutTriangularSumDeltaPhiLepMET()             -- Parameterized Constructor          *
 *    ~CutTriangularSumDeltaPhiLepMET()            -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutTriangularSumDeltaPhiLepMET"      *
 *                                                                            *
 * Private Data Members of CutTriangularSumDeltaPhiLepMET                                *
 *    myTH2F* _hTriangularSumDeltaPhiLepMETBefore         -- Hist before cuts            *
 *    myTH2F* _hTriangularSumDeltaPhiLepMETAfter          -- Hist after cuts             *
 *    Int_t _TriangularSumDeltaPhiLepMETMin;              -- Minimum requirement         *
 *    Int_t _TriangularSumDeltaPhiLepMETMax;              -- Maximum requirement         *
 *                                                                            *
 * History                                                                    *
 *      8 Jun 2011, created by Huaqiao ZHANG @ CERN                           *
 *****************************************************************************/

#ifndef CutTriangularSumDeltaPhiLepMET_h
#define CutTriangularSumDeltaPhiLepMET_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutTriangularSumDeltaPhiLepMET : public HistoCut 
{

public:

  // Parameterized Constructor
  CutTriangularSumDeltaPhiLepMET(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutTriangularSumDeltaPhiLepMET();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutTriangularSumDeltaPhiLepMET"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();

  double xy_to_phi(double x, double y);
  double dPhi(double phi1,double phi2);

private:

  myTH2F* _hTriangularSumDeltaPhiLepMET2DBefore;   // 2D Histogram before the cut
  myTH2F* _hTriangularSumDeltaPhiLepMET2DAfter;    // 2D Histogram after the cut

  // Cut parameters
  Double_t _TriangularSumDeltaPhiLepMETMin;           // Minimum condition
  Double_t _TriangularSumDeltaPhiLepMETMax;           // Maximum condition


};


#endif







