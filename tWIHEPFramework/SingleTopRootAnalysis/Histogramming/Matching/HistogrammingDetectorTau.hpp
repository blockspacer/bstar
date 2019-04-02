/******************************************************************************
 * HistogrammingDetectorTau.hpp                                               *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingDetectorTau()    -- Parameterized Constructor              *
 *    ~HistogrammingDetectorTau()   -- Destructor                             *
 *    BookHistogram()                         -- Book histograms              *
 *    Apply()                          -- Fill histograms only (No Cuts)      *
 *    GetCutName()           -- Returns "HistogrammingDetectorTau"            *
 *                                                                            *
 * Private Data Members of HistogrammingDetectorTau class                     *
 * - lots of histograms for detector for taus                                 *
 *                                                                            *
 * History                                                                    *
 *           Apr 10 2007 - Created by J. Holzbauer                            *
 *****************************************************************************/

#ifndef HistogrammingDetectorTau_h
#define HistogrammingDetectorTau_h


#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingDetectorTau : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingDetectorTau(EventContainer *obj);
  
  // Destructor
  ~HistogrammingDetectorTau();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingDetectorTau"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  //here, particle is a tau
  //miss is where there is a MCTau but no Tau- the detector missed the particle
  //false is where there is a Tau but no MCTau- the detector falsely ID'ed the particle
 
  myTH1F* _hmissMCetaTau;     //eta for particles missed by detector in events where there are no reco particles
  myTH1F* _hfalseMCetaTau;    //eta for particles falsely identified by detector in events where there are no mc particles
  myTH1F* _hmissMCetaTauAll;   //eta for all MC Taus in events with at least one MCTau missed by detector
  myTH1F* _hfalseMCetaTauAll;  //eta for all Reco Taus in events with at least one Tau falsely identified by detector
  
  myTH1F* _hmissMCptTau;  //pt for particles missed by detector in events where there are no reco particles
  myTH1F* _hfalseMCptTau;   //pt for particles falsely identified by detector in events where there are no mc particles
  myTH1F* _hmissMCptTauAll; //pt for all MC Taus in events with at least one MCTau missed by detector
  myTH1F* _hfalseMCptTauAll;  //pt for all Reco Taus in events with at least one Tau falsely identified by detector
  
  myTH1F* _hdeltaRTau;     //deltaR for MCTau[i] and Tau[i], as long as both exist
  myTH1F* _hdeltaPtTau;     //deltaPt for MCTau[i] and Tau[i], as long as both exist
  
  TH1F* _hsummaryTau; //total number of events that contain at least one: reco tau, false tau, MCTau, missed MCTau

};

#endif
