/******************************************************************************
 * TriggerWeightConst.hpp                                                            *
 *                                                                            *
 * Weight events in the histograms                                            *
 *                                                                            *
 * This class is derived from the BaseCut Class and simply                    *
 * cuts on if the event has reconstructed objects                             *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    TriggerWeightConst()                     -- Parameterized Constructor          *
 *    ~TriggerWeightConst()                    -- Destructor                         *
 *    BookHistogram()                   -- Book histograms of weights         *
 *    Apply()                           -- Apply weight and fill histograms   *
 *    GetCutName()                      -- Returns "TriggerWeightConst"              *
 *                                                                            *
 * Private Data Members of AnalysisMain class                                 *
 *    myTH1F* _hTreeWeight              -- Histogram of Input Tree Weights    *
 *    myTH1F* _hGlobalWeight            -- Histogram of Global Weights        *
 *    myTH1F* _hOutputWeight            -- Histogram of Output Weights        *
 *                                                                            *
 * History                                                                    *
 *         12 Feb 08 - created by J. Holzbauer, based on EventWeight          *
 *                     This class weights by a constant only                  *
 *****************************************************************************/

#ifndef TriggerWeightConst_h
#define TriggerWeightConst_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class TriggerWeightConst : public HistoCut
{

public:

  // Parameterized Constructor
  TriggerWeightConst(EventContainer *obj);
  
  // Default Destructor
  ~TriggerWeightConst();

  // Book Histograms
  void BookHistogram(void );
  
  // Apply weight and fill histograms
  Bool_t Apply();

   
  // Get the name describing the cut
  inline std::string GetCutName() { return "TriggerWeightConst"; };


private:

  // Histograms
  myTH1F* _hTriggerWeight;   // Histogram of input tree weights
  myTH1F* _hOutputWeight;   // Histogram of input tree weights
  double _ConstTriggerWeight; //the constant weight
};


#endif







