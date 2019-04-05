/******************************************************************************
 * HistogrammingDetectorElectron.hpp                                          *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingDetectorElectron()    -- Parameterized Constructor         *
 *    ~HistogrammingDetectorElectron()   -- Destructor                        *
 *    BookHistogram()                         -- Book histograms              *
 *    Apply()                          -- Fill histograms only (No Cuts)      *
 *    GetCutName()           -- Returns "HistogrammingDetectorElectron"       *
 *                                                                            *
 * Private Data Members of HistogrammingDetectorElectron class                *
 * - lots of histograms for detector for electrons                            *
 *                                                                            *
 * History                                                                    *
 *           Apr 10 2007 - Created by J. Holzbauer                            *
 *****************************************************************************/

#ifndef HistogrammingDetectorElectron_h
#define HistogrammingDetectorElectron_h


#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingDetectorElectron : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingDetectorElectron(EventContainer *obj);
  
  // Destructor
  ~HistogrammingDetectorElectron();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingDetectorElectron"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();
  
private:
// Histograms declarations 
  //here, particle is a electron
  //miss is where there is a MCElectron but no Electron- the detector missed the particle
  //false is where there is a Electron but no MCElectron- the detector falsely ID'ed the particle
 
  myTH1F* _hmissMCetaEle;     //eta for particles missed by detector in events where there are no reco particles
  myTH1F* _hfalseMCetaEle;    //eta for particles falsely identified by detector in events where there are no mc particles
  myTH1F* _hmissMCetaEleAll;   //eta for all MC Electrons in events with at least one MCElectron missed by detector
  myTH1F* _hfalseMCetaEleAll;  //eta for all Reco Electrons in events with at least one Electron falsely identified by detector
  
  myTH1F* _hmissMCptEle;  //pt for particles missed by detector in events where there are no reco particles
  myTH1F* _hfalseMCptEle;   //pt for particles falsely identified by detector in events where there are no mc particles
  myTH1F* _hmissMCptEleAll; //pt for all MC Electrons in events with at least one MCElectron missed by detector
  myTH1F* _hfalseMCptEleAll;  //pt for all Reco Electrons in events with at least one Electron falsely identified by detector
  
  myTH1F* _hdeltaREle;     //deltaR for MCElectron[i] and Electron[i], as long as both exist
  myTH1F* _hdeltaPtEle;     //deltaPt for MCElectron[i] and Electron[i], as long as both exist
  
  TH1F* _hsummaryEle; //total number of events that contain at least one: reco electron, false electron, MCElectron, missed MCElectron

};

#endif
