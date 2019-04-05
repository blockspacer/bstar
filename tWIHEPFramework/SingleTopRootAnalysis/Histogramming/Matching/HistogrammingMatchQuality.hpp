/******************************************************************************
 * HistogrammingMatchQuality.hpp                                              *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingMatchQuality()      -- Parameterized Constructor           *
 *    ~HistogrammingMatchQuality()     -- Destructor                          *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms                    *
 *    GetCutName()            -- Returns "HistogrammingMatchQuality"          *
 *                                                                            *
 * Private Data Members of this class                                         *
 * - lots of histograms of particle properties                                *
 *                                                                            *
 * Purpose: makes plots of deltaR, deltaPt, DeltaEta, and DeltaPhi.  Also     *
 *   produces these for various Pt ranges from 10-75 GeV, intervals of 5 GeV. *
 *   Makes plots for all possible combinations of the Reco and MC particles   *
 *   of a particular type (Mu or Ele for now) in an event.                    *
 *                                                                            *
 * History                                                                    *
 *         23 May 2007-- created by J. Holzbaer                               *
 *****************************************************************************/

#ifndef HistogrammingMatchQuality_h
#define HistogrammingMatchQuality_h


#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class HistogrammingMatchQuality : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingMatchQuality(EventContainer *obj);
  
  // Destructor
  ~HistogrammingMatchQuality();

  // Get the name describing the cut
  inline std::string GetCutName() { return "HistogrammingMatchQuality"; };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();

  

  //the following are used mostly in booking the histogram but also in drawing it
  //so they are declared here.

  const static int binMudetapt = 150;
  const static double xminMudetapt = -0.005, xmaxMudetapt = 0.005;
  const static int binEledetapt = 150;
  const static double xminEledetapt = -0.005, xmaxEledetapt = 0.005;
 
   const static int   binMudphipt = 300;
  const static double xminMudphipt = -0.005, xmaxMudphipt = 0.005;
  const static int   binEledphipt = 300;
  const static double xminEledphipt = -0.005, xmaxEledphipt = 0.005;
 
   //TDR parameters- -.5 to .5 at 100 bins
  const static int binMudptpt = 250;
  const static double xminMudptpt = -0.25, xmaxMudptpt = 0.25;
  const static int binEledptpt = 250;
  const static double xminEledptpt = -0.25, xmaxEledptpt = 0.25;

  const static int binMudpteta = 250;
  const static double xminMudpteta = -0.25, xmaxMudpteta = 0.25;
  const static int binEledpteta = 250;
  const static double xminEledpteta = -0.25, xmaxEledpteta = 0.25;

  const static int binMudetaeta = 150;
  const static double xminMudetaeta = -0.005, xmaxMudetaeta = 0.005;
  const static int binEledetaeta = 150;
  const static double xminEledetaeta = -0.005, xmaxEledetaeta = 0.005;
  
  const static int binMudphiphi = 300;
  const static double xminMudphiphi = -0.005, xmaxMudphiphi = 0.005;
  const static int binEledphiphi = 300;
  const static double xminEledphiphi = -0.005, xmaxEledphiphi = 0.005;

  const static int binMudrpt = 200;
  const static double xminMudrpt = 0.0, xmaxMudrpt = 0.01;
  const static int binEledrpt = 200;
  const static double xminEledrpt = 0.0, xmaxEledrpt = 0.01;
 
 
  const static int  binbkgdMudetapt = 100;
  const static double xminbkgdMudetapt = -5, xmaxbkgdMudetapt = 5;
  const static int  binbkgdEledetapt = 100;
  const static double xminbkgdEledetapt = -5, xmaxbkgdEledetapt = 5;

  const static int  binbkgdMudphipt = 80;
  const static double xminbkgdMudphipt = -6, xmaxbkgdMudphipt = 6;
  const static int  binbkgdEledphipt = 80;
  const static double xminbkgdEledphipt = -6, xmaxbkgdEledphipt = 6;
  
  const static int  binbkgdMudptpt = 50;
  const static double xminbkgdMudptpt = -1.5, xmaxbkgdMudptpt = 6.5;
  const static int   binbkgdEledptpt = 50;
  const static double xminbkgdEledptpt = -1.5, xmaxbkgdEledptpt = 6.5;

  const static int binbkgdMudpteta = 50;
  const static double xminbkgdMudpteta = -1.5, xmaxbkgdMudpteta = 6.5;
  const static int binbkgdEledpteta = 50;
  const static double xminbkgdEledpteta = -1.5, xmaxbkgdEledpteta = 6.5;

  const static int binbkgdMudetaeta = 100;
  const static double xminbkgdMudetaeta = -5, xmaxbkgdMudetaeta = 5;
  const static int binbkgdEledetaeta = 100;
  const static double xminbkgdEledetaeta = -5, xmaxbkgdEledetaeta = 5;

  const static int binbkgdMudphiphi = 80;
  const static double xminbkgdMudphiphi = -6, xmaxbkgdMudphiphi = 6;
  const static int binbkgdEledphiphi = 80;
  const static double xminbkgdEledphiphi = -6, xmaxbkgdEledphiphi = 6;

  const static double xminbkgdMudrpt = 0.01, xmaxbkgdMudrpt = 0.1;
  const static int  binbkgdEledrpt = 200;
  const static double xminbkgdEledrpt = 0.01, xmaxbkgdEledrpt = 0.1;
  const static int  binbkgdMudrpt = 200;
 
 /**
  
  const static int  binbkgdMudetapt = 18000;
  const static double xminbkgdMudetapt = -6, xmaxbkgdMudetapt = 6;
  const static int  binbkgdEledetapt = 18000;
  const static double xminbkgdEledetapt = -6, xmaxbkgdEledetapt = 6;

  const static int  binbkgdMudphipt = 24000;
  const static double xminbkgdMudphipt = -4, xmaxbkgdMudphipt = 4;
  const static int  binbkgdEledphipt = 24000;
  const static double xminbkgdEledphipt = -4, xmaxbkgdEledphipt = 4;
  
  const static int  binbkgdMudptpt = 2000;
  const static double xminbkgdMudptpt = -1.5, xmaxbkgdMudptpt = 6.5;
  const static int   binbkgdEledptpt = 2000;
  const static double xminbkgdEledptpt = -1.5, xmaxbkgdEledptpt = 6.5;

  const static int binbkgdMudpteta = 4000;
  const static double xminbkgdMudpteta = -1.5, xmaxbkgdMudpteta = 14.5;
  const static int binbkgdEledpteta = 4000;
  const static double xminbkgdEledpteta = -1.5, xmaxbkgdEledpteta = 14.5;

  const static int binbkgdMudetaeta = 18000;
  const static double xminbkgdMudetaeta = -6, xmaxbkgdMudetaeta = 6;
  const static int binbkgdEledetaeta = 18000;
  const static double xminbkgdEledetaeta = -6, xmaxbkgdEledetaeta = 6;

  const static int binbkgdMudphiphi = 24000;
  const static double xminbkgdMudphiphi = -4, xmaxbkgdMudphiphi = 4;
  const static int binbkgdEledphiphi = 24000;
  const static double xminbkgdEledphiphi = -4, xmaxbkgdEledphiphi = 4;

  const static double xminbkgdMudrpt = 0.01, xmaxbkgdMudrpt = 0.1;
  const static int  binbkgdEledrpt = 200;
  const static double xminbkgdEledrpt = 0.01, xmaxbkgdEledrpt = 0.1;
  const static int  binbkgdMudrpt = 200;
  **/
 /**
  const static int  binbkgdMudetapt = 100;
  const static double xminbkgdMudetapt = -6, xmaxbkgdMudetapt = 6;
  const static int  binbkgdEledetapt = 100;
  const static double xminbkgdEledetapt = -6, xmaxbkgdEledetapt = 6;

  const static int  binbkgdMudphipt = 70;
  const static double xminbkgdMudphipt = -7, xmaxbkgdMudphipt = 7;
  const static int  binbkgdEledphipt = 70;
  const static double xminbkgdEledphipt = -7, xmaxbkgdEledphipt = 7;
  
  const static int  binbkgdMudptpt = 100;
  const static double xminbkgdMudptpt = -1.5, xmaxbkgdMudptpt = 7;
  const static int   binbkgdEledptpt = 100;
  const static double xminbkgdEledptpt = -1.5, xmaxbkgdEledptpt = 7;

  const static int binbkgdMudpteta = 200;
  const static double xminbkgdMudpteta = -1.5, xmaxbkgdMudpteta = 30;
  const static int binbkgdEledpteta = 200;
  const static double xminbkgdEledpteta = -1.5, xmaxbkgdEledpteta = 25;

  const static int binbkgdMudetaeta = 100;
  const static double xminbkgdMudetaeta = -6, xmaxbkgdMudetaeta = 6;
  const static int binbkgdEledetaeta = 100;
  const static double xminbkgdEledetaeta = -6, xmaxbkgdEledetaeta = 6;

  const static int binbkgdMudphiphi = 70;
  const static double xminbkgdMudphiphi = -7, xmaxbkgdMudphiphi = 7;
  const static int binbkgdEledphiphi = 70;
  const static double xminbkgdEledphiphi = -7, xmaxbkgdEledphiphi = 7;

  const static double xminbkgdMudrpt = 0.01, xmaxbkgdMudrpt = 0.1;
  const static int  binbkgdEledrpt = 200;
  const static double xminbkgdEledrpt = 0.01, xmaxbkgdEledrpt = 0.1;
  const static int  binbkgdMudrpt = 200;
 **/
private:

  // Histograms declarations 
 
  myTH1F* _hdeltaRMuall;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMuall;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREleall;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEleall;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaPtMuodd;     //deltaPt for odd posible MC Reco matches
  myTH1F* _hdeltaPtEleodd;     //deltaPt for odd posible MC Reco matches
  myTH1F* _hdeltaEtaMuodd;     //deltaPt for odd posible MC Reco matches
  myTH1F* _hdeltaEtaEleodd;     //deltaPt for odd posible MC Reco matches
  myTH1F* _hdeltaPhiMuodd;     //deltaPt for odd posible MC Reco matches
  myTH1F* _hdeltaPhiEleodd;     //deltaPt for odd posible MC Reco matches

  myTH1F* _hdeltaEtaMuall;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaEtaEleall;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaPhiMuall;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaPhiEleall;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaRMu10to15;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu10to15;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle10to15;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle10to15;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaRMu15to20;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu15to20;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle15to20;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle15to20;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaRMu20to25;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu20to25;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle20to25;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle20to25;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaRMu25to30;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu25to30;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle25to30;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle25to30;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaRMu30to35;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu30to35;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle30to35;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle30to35;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaRMu35to40;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu35to40;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle35to40;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle35to40;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaRMu40to45;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu40to45;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle40to45;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle40to45;     //deltaPt for all posible MC Reco matches
  
  myTH1F* _hdeltaRMu45to50;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu45to50;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle45to50;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle45to50;     //deltaPt for all posible MC Reco matches
  
  myTH1F* _hdeltaRMu50to55;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu50to55;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle50to55;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle50to55;     //deltaPt for all posible MC Reco matches
 
  myTH1F* _hdeltaRMu55to60;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu55to60;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle55to60;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle55to60;     //deltaPt for all posible MC Reco matches
 
  myTH1F* _hdeltaRMu60to65;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu60to65;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle60to65;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle60to65;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaRMu65to70;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu65to70;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle65to70;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle65to70;     //deltaPt for all posible MC Reco matches

  myTH1F* _hdeltaRMu70to75;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtMu70to75;     //deltaPt for all posible MC Reco matches
  myTH1F* _hdeltaREle70to75;     //deltaR for all posible MC Reco matches
  myTH1F* _hdeltaPtEle70to75;     //deltaPt for all posible MC Reco matches

 myTH1F* _hbkgddeltaRMuall;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMuall;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREleall;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEleall;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMuall;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEleall;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMuall;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEleall;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaRMu10to15;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu10to15;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle10to15;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle10to15;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaRMu15to20;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu15to20;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle15to20;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle15to20;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaRMu20to25;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu20to25;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle20to25;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle20to25;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaRMu25to30;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu25to30;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle25to30;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle25to30;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaRMu30to35;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu30to35;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle30to35;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle30to35;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaRMu35to40;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu35to40;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle35to40;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle35to40;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaRMu40to45;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu40to45;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle40to45;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle40to45;     //bkgddeltaPt for all posible MC Reco matches
  
  myTH1F* _hbkgddeltaRMu45to50;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu45to50;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle45to50;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle45to50;     //bkgddeltaPt for all posible MC Reco matches
  
  myTH1F* _hbkgddeltaRMu50to55;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu50to55;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle50to55;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle50to55;     //bkgddeltaPt for all posible MC Reco matches
 
  myTH1F* _hbkgddeltaRMu55to60;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu55to60;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle55to60;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle55to60;     //bkgddeltaPt for all posible MC Reco matches
 
  myTH1F* _hbkgddeltaRMu60to65;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu60to65;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle60to65;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle60to65;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaRMu65to70;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu65to70;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle65to70;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle65to70;     //bkgddeltaPt for all posible MC Reco matches

  myTH1F* _hbkgddeltaRMu70to75;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtMu70to75;     //bkgddeltaPt for all posible MC Reco matches
  myTH1F* _hbkgddeltaREle70to75;     //bkgddeltaR for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEle70to75;     //bkgddeltaPt for all posible MC Reco matches
 
  myTH1F* _hdeltaEtaMu10to15;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu10to15;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle10to15;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle10to15;     //deltaPhi for all posible MC Reco matches

  myTH1F* _hdeltaEtaMu15to20;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu15to20;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle15to20;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle15to20;     //deltaPhi for all posible MC Reco matches

  myTH1F* _hdeltaEtaMu20to25;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu20to25;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle20to25;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle20to25;     //deltaPhi for all posible MC Reco matches

  myTH1F* _hdeltaEtaMu25to30;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu25to30;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle25to30;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle25to30;     //deltaPhi for all posible MC Reco matches

  myTH1F* _hdeltaEtaMu30to35;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu30to35;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle30to35;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle30to35;     //deltaPhi for all posible MC Reco matches

  myTH1F* _hdeltaEtaMu35to40;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu35to40;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle35to40;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle35to40;     //deltaPhi for all posible MC Reco matches

  myTH1F* _hdeltaEtaMu40to45;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu40to45;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle40to45;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle40to45;     //deltaPhi for all posible MC Reco matches
  
  myTH1F* _hdeltaEtaMu45to50;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu45to50;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle45to50;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle45to50;     //deltaPhi for all posible MC Reco matches
  
  myTH1F* _hdeltaEtaMu50to55;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu50to55;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle50to55;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle50to55;     //deltaPhi for all posible MC Reco matches
 
  myTH1F* _hdeltaEtaMu55to60;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu55to60;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle55to60;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle55to60;     //deltaPhi for all posible MC Reco matches
 
  myTH1F* _hdeltaEtaMu60to65;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu60to65;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle60to65;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle60to65;     //deltaPhi for all posible MC Reco matches

  myTH1F* _hdeltaEtaMu65to70;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu65to70;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle65to70;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle65to70;     //deltaPhi for all posible MC Reco matches

  myTH1F* _hdeltaEtaMu70to75;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiMu70to75;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaEtaEle70to75;     //deltaEta for all posible MC Reco matches
  myTH1F* _hdeltaPhiEle70to75;     //deltaPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMu10to15;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu10to15;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle10to15;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle10to15;     //bkgddeltaPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMu15to20;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu15to20;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle15to20;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle15to20;     //bkgddeltaPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMu20to25;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu20to25;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle20to25;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle20to25;     //bkgddeltaPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMu25to30;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu25to30;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle25to30;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle25to30;     //bkgddeltaPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMu30to35;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu30to35;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle30to35;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle30to35;     //bkgddeltaPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMu35to40;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu35to40;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle35to40;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle35to40;     //bkgddeltaPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMu40to45;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu40to45;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle40to45;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle40to45;     //bkgddeltaPhi for all posible MC Reco matches
  
  myTH1F* _hbkgddeltaEtaMu45to50;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu45to50;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle45to50;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle45to50;     //bkgddeltaPhi for all posible MC Reco matches
  
  myTH1F* _hbkgddeltaEtaMu50to55;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu50to55;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle50to55;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle50to55;     //bkgddeltaPhi for all posible MC Reco matches
 
  myTH1F* _hbkgddeltaEtaMu55to60;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu55to60;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle55to60;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle55to60;     //bkgddeltaPhi for all posible MC Reco matches
 
  myTH1F* _hbkgddeltaEtaMu60to65;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu60to65;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle60to65;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle60to65;     //bkgddeltaPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMu65to70;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu65to70;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle65to70;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle65to70;     //bkgddeltaPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaMu70to75;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiMu70to75;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEle70to75;     //bkgddeltaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiEle70to75;     //bkgddeltaPhi for all posible MC Reco matches

   myTH1F* _hdeltaEtaEtaMuN4toN3_5;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEleN4toN3_5;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMuN3_5toN3;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEleN3_5toN3;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMuN3toN2_5;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEleN3toN2_5;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMuN2_5toN2;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEleN2_5toN2;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMuN2toN1_5;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEleN2toN1_5;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMuN1_5toN1;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEleN1_5toN1;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMuN1toN0_5;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEleN1toN0_5;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMuN0_5to0;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEleN0_5to0;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMu0to0_5;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEle0to0_5;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMu0_5to1;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEle0_5to1;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMu1to1_5;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEle1to1_5;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMu1_5to2;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEle1_5to2;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMu2to2_5;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEle2to2_5;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMu2_5to3;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEle2_5to3;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMu3to3_5;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEle3to3_5;     //deltaEtaEta for all posible MC Reco matches
 myTH1F* _hdeltaEtaEtaMu3_5to4;     //deltaEtaEta for all posible MC Reco matches
  myTH1F* _hdeltaEtaEtaEle3_5to4;     //deltaEtaEta for all posible MC Reco matches

  myTH1F* _hdeltaPtEtaMuN4toN3_5;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEleN4toN3_5;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMuN3_5toN3;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEleN3_5toN3;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMuN3toN2_5;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEleN3toN2_5;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMuN2_5toN2;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEleN2_5toN2;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMuN2toN1_5;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEleN2toN1_5;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMuN1_5toN1;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEleN1_5toN1;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMuN1toN0_5;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEleN1toN0_5;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMuN0_5to0;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEleN0_5to0;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMu0to0_5;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEle0to0_5;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMu0_5to1;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEle0_5to1;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMu1to1_5;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEle1to1_5;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMu1_5to2;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEle1_5to2;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMu2to2_5;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEle2to2_5;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMu2_5to3;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEle2_5to3;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMu3to3_5;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEle3to3_5;     //deltaPtEta for all posible MC Reco matches
 myTH1F* _hdeltaPtEtaMu3_5to4;     //deltaPtEta for all posible MC Reco matches
  myTH1F* _hdeltaPtEtaEle3_5to4;     //deltaPtEta for all posible MC Reco matches

 myTH1F* _hdeltaPhiPhiMuN3_5toN3;     //deltaPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEleN3_5toN3;     //deltaPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMuN3toN2_5;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEleN3toN2_5;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMuN2_5toN2;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEleN2_5toN2;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMuN2toN1_5;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEleN2toN1_5;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMuN1_5toN1;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEleN1_5toN1;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMuN1toN0_5;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEleN1toN0_5;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMuN0_5to0;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEleN0_5to0;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMu0to0_5;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEle0to0_5;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMu0_5to1;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEle0_5to1;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMu1to1_5;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEle1to1_5;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMu1_5to2;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEle1_5to2;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMu2to2_5;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEle2to2_5;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMu2_5to3;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEle2_5to3;     //deltaPhiPhi for all posible MC Reco matches
 myTH1F* _hdeltaPhiPhiMu3to3_5;     //deltaPhiPhi for all posible MC Reco matches
  myTH1F* _hdeltaPhiPhiEle3to3_5;     //deltaPhiPhi for all posible MC Reco matches

  myTH1F* _hbkgddeltaEtaEtaMuN4toN3_5;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEleN4toN3_5;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMuN3_5toN3;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEleN3_5toN3;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMuN3toN2_5;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEleN3toN2_5;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMuN2_5toN2;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEleN2_5toN2;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMuN2toN1_5;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEleN2toN1_5;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMuN1_5toN1;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEleN1_5toN1;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMuN1toN0_5;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEleN1toN0_5;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMuN0_5to0;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEleN0_5to0;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMu0to0_5;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEle0to0_5;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMu0_5to1;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEle0_5to1;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMu1to1_5;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEle1to1_5;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMu1_5to2;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEle1_5to2;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMu2to2_5;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEle2to2_5;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMu2_5to3;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEle2_5to3;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMu3to3_5;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEle3to3_5;     //bkgddeltaEtaEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaEtaEtaMu3_5to4;     //bkgddeltaEtaEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaEtaEtaEle3_5to4;     //bkgddeltaEtaEta for all posible MC Reco matches

  myTH1F* _hbkgddeltaPtEtaMuN4toN3_5;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEleN4toN3_5;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMuN3_5toN3;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEleN3_5toN3;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMuN3toN2_5;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEleN3toN2_5;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMuN2_5toN2;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEleN2_5toN2;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMuN2toN1_5;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEleN2toN1_5;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMuN1_5toN1;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEleN1_5toN1;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMuN1toN0_5;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEleN1toN0_5;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMuN0_5to0;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEleN0_5to0;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMu0to0_5;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEle0to0_5;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMu0_5to1;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEle0_5to1;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMu1to1_5;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEle1to1_5;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMu1_5to2;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEle1_5to2;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMu2to2_5;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEle2to2_5;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMu2_5to3;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEle2_5to3;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMu3to3_5;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEle3to3_5;     //bkgddeltaPtEta for all posible MC Reco matches
 myTH1F* _hbkgddeltaPtEtaMu3_5to4;     //bkgddeltaPtEta for all posible MC Reco matches
  myTH1F* _hbkgddeltaPtEtaEle3_5to4;     //bkgddeltaPtEta for all posible MC Reco matches

 myTH1F* _hbkgddeltaPhiPhiMuN3_5toN3;     //bkgddeltaPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEleN3_5toN3;     //bkgddeltaPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMuN3toN2_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEleN3toN2_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMuN2_5toN2;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEleN2_5toN2;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMuN2toN1_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEleN2toN1_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMuN1_5toN1;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEleN1_5toN1;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMuN1toN0_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEleN1toN0_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMuN0_5to0;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEleN0_5to0;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMu0to0_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEle0to0_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMu0_5to1;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEle0_5to1;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMu1to1_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEle1to1_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMu1_5to2;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEle1_5to2;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMu2to2_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEle2to2_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMu2_5to3;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEle2_5to3;     //bkgddeltaPhiPhi for all posible MC Reco matches
 myTH1F* _hbkgddeltaPhiPhiMu3to3_5;     //bkgddeltaPhiPhi for all posible MC Reco matches
  myTH1F* _hbkgddeltaPhiPhiEle3to3_5;     //bkgddeltaPhiPhi for all posible MC Reco matches


};

#endif
