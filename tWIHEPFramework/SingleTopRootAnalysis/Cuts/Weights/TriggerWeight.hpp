/******************************************************************************
 * TriggerWeight.hpp                                                            *
 *                                                                            *
 * Weight events in the histograms                                            *
 *                                                                            *
 * This class is derived from the BaseCut Class and simply                    *
 * cuts on if the event has reconstructed objects                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    TriggerWeight()                     -- Parameterized Constructor          *
 *    ~TriggerWeight()                    -- Destructor                         *
 *    BookHistogram()                   -- Book histograms of weights         *
 *    Apply()                           -- Apply weight and fill histograms   *
 *    GetCutName()                      -- Returns "TriggerWeight"              *
 *                                                                            *
 * Private Data Members of AnalysisMain class                                 *
 *    myTH1F* _hTreeWeight              -- Histogram of Input Tree Weights    *
 *    myTH1F* _hGlobalWeight            -- Histogram of Global Weights        *
 *    myTH1F* _hOutputWeight            -- Histogram of Output Weights        *
 *                                                                            *
 * History                                                                    *
 *     28 Apr 2005 - Created by R. Schwienhorst at D0                         *
 *      7 Nov 2006 - Modified by R. Schwienhorst for ATLAS                    *
 *      9 Nov 2006 - Modified by P. Ryan - cleanup and reorganized            *
 *****************************************************************************/

#ifndef TriggerWeight_h
#define TriggerWeight_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class TriggerWeight : public HistoCut
{

public:

  // Parameterized Constructor
  TriggerWeight(EventContainer *obj);
  
  // Default Destructor
  ~TriggerWeight();

  // Book Histograms
  void BookHistogram(void );
  
  // Apply weight and fill histograms
  Bool_t Apply();

   
  // Get the name describing the cut
  inline std::string GetCutName() { return "TriggerWeight"; };


private:

  // Histograms
  myTH1F* _hTriggerWeight;   // Histogram of input tree weights
  myTH1F* _hOutputWeight;   // Histogram of input tree weights
  TF1 *_hMuCentral;
  TF1 *_hMuICR;
  TF1 *_hMuForward;
  TF1 *_hEMCentral;
  TF1 *_hEMICR;
  TF1 *_hEMForward;
  TFile *_f1;
};


#endif







