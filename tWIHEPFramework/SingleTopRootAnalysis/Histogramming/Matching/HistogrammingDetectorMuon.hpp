/******************************************************************************
 * HistogrammingDetectorMuon.hpp                                              *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingDetectorMuon()    -- Parameterized Constructor             *
 *    ~HistogrammingDetectorMuon()   -- Destructor                            *
 *    BookHistogram()                         -- Book histograms              *
 *    Apply()                          -- Fill histograms only (No Cuts)      *
 *    GetCutName()           -- Returns "HistogrammingDetectorMuon"           *
 *                                                                            *
 * Private Data Members of HistogrammingDetectorMuon class                    *
 * - lots of histograms for detector for muons                                *
 *                                                                            *
 * History                                                                    *
 *           Apr 10 2007 - Created by J. Holzbauer                            *
 *****************************************************************************/

#ifndef HistogrammingDetectorMuon_h
#define HistogrammingDetectorMuon_h


#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingDetectorMuon : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingDetectorMuon(EventContainer *obj);
  
  // Destructor
  ~HistogrammingDetectorMuon();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingDetectorMuon"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
  // Histograms declarations 
  //here, particle is a muon
  //miss is where there is a MCMuon but no Muon- the detector missed the particle
  //false is where there is a Muon but no MCMuon- the detector falsely ID'ed the particle
 
  myTH1F* _hmissMCetaMu;     //eta for particles missed by detector in events where there are no reco particles
  myTH1F* _hfalseMCetaMu;    //eta for particles falsely identified by detector in events where there are no mc particles
  myTH1F* _hmissMCetaMuAll;   //eta for all MC Muons in events with at least one MCMuon missed by detector
  myTH1F* _hfalseMCetaMuAll;  //eta for all Reco Muons in events with at least one Muon falsely identified by detector
  
  myTH1F* _hmissMCptMu;  //pt for particles missed by detector in events where there are no reco particles
  myTH1F* _hfalseMCptMu;   //pt for particles falsely identified by detector in events where there are no mc particles
  myTH1F* _hmissMCptMuAll; //pt for all MC Muons in events with at least one MCMuon missed by detector
  myTH1F* _hfalseMCptMuAll;  //pt for all Reco Muons in events with at least one Muon falsely identified by detector
  
  myTH1F* _hdeltaRMu;     //deltaR for MCMuon[i] and Muon[i], as long as both exist
  myTH1F* _hdeltaPtMu;     //deltaPt for MCMuon[i] and Muon[i], as long as both exist
  
  TH1F* _hsummaryMu; //total number of events that contain at least one: reco muon, false muon, MCMuon, missed MCMuon

};

#endif
