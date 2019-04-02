/******************************************************************************
 * CutNN.hpp                                                               *
 *                                                                            *
 * Cuts on neural network variable                                            *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutNJet class                                   *
 *    CutNN()                        -- Parameterized Constructor          *
 *    ~CutNN()                       -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutNN"                 *
 *                                                                            *
 * Private Data Members of CutNJet class                                      *
 *    myTH1F* _hNJetsBefore              -- Hist of Num of jets before cuts   *
 *    myTH1F* _hNJetsAfter               -- Hist of Num of jets before cuts   *
 *    Int_t _nJetsMin;                   -- Minimum number of jets to require *
 *    Int_t _nJetsMax;                   -- Maximum number of jets allowed    *
 *                                                                            *
 * History                                                                    *
 *       7 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      10 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *      21 Nov 2006 - Modified by N. Triplett - changed to NNCut              *
 *****************************************************************************/
#ifndef CutNN_h
#define CutNN_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include "TMultiLayerPerceptron.h"

#include <string>
#include <vector>
#include <sstream>

using namespace std;


class CutNN : public HistoCut 
{

public:

  // Parameterized Constructor
  CutNN(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutNN();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutNN"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  TTree * simu;
  TMultiLayerPerceptron * mlp;

  // Histograms
  myTH1F* NNOutput;    // Histogram NN output                        
  myTH1F* NNOutputPost;  // Histogram NN output after the cut                       

  // Cut parameters
  Double_t _NNValueMin;       // Min value of the NN output required


};


#endif







