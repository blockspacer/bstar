/******************************************************************************
 * CutPrimaryVertex.hpp                                                         *
 *                                                                            *
 * Apply a default set of cuts to all classes deriving from BaseCuts          *
 *                                                                            *
 * Derived from BaseCut                                                       *
 *                                                                            *
 * Public Member Functions of CutPrimaryVertex class                            *
 *    CutPrimaryVertex()                        -- Parameterized Constructor    *
 *    ~CutPrimaryVertex()                       -- Destructor                   *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutPrimaryVertex"           *
 *                                                                            *

 * History                                                                    *
 *      16 June 2011 - Created by J. Holzbauer                                     *
 ******************************************************************************/

#ifndef CutPrimaryVertex_h
#define CutPrimaryVertex_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>
class CutPrimaryVertex : public HistoCut
{

public:
  
  // Parameterized constructor
  CutPrimaryVertex(EventContainer  *EventContainerObj);
  
  // Destructor
  ~CutPrimaryVertex();

  // Book histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Return name of Cut
  inline std::string GetCutName() { return "CutPrimaryVertex"; };
  

private:

  // Histograms
  myTH1F* _hPrimaryVertexBefore;     // Histogram Before GRL Cut
  myTH1F* _hPrimaryVertexAfter;      // Histogram After GRL Cut
  myTH1F* _hNTrueInteractions;       //Histogram containing the true number of interactions
  myTH1F* _hNTrueIntsUnweighted;       //Histogram containing the true number of interactions without weight applide. this is for PU reweighting purposes.
  

  // Variables to cut PVs on 
  Double_t _ndofCut;
  Double_t _zCut;
  Double_t _dxyCut;
  
  // Tmp
  Int_t tmpCounter;
    
  
  Int_t PrimaryVertexNumber;

};


#endif







