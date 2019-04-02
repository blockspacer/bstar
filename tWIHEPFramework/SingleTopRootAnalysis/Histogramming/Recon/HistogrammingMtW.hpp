/******************************************************************************
 * HistogrammingMtW.hpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMtW()              -- Parameterized Constructor            *
 *    ~HistogrammingMtW()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingMtW"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *  - histograms of MwT                                                       *
 *  - Will probably also include other histograms to be filled at all stages  *
 *    here                                                                    *
 *                                                                            *
 * History                                                                    *
 *      08 Jul 2016 - Created by D. Leggat for CMS IHEP                       *
 *****************************************************************************/

#ifndef HistogrammingMtW_h
#define HistogrammingMtW_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMtW : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMtW(EventContainer *obj, bool unisolated = false);
  
  // Destructor
  ~HistogrammingMtW();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMtW"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  myTH1F* _hMtW;  // Histogram of MtW

  Particle lepton; // Assign the TLorentzVector of the primary lepton to this lepton in order to correctly calculate the mtw in case of unisolated leptons.

  bool _unisolated; //Whether we are using unisolated leptons

};


#endif







