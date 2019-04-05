/******************************************************************************
 * HistogrammingTagEfficiency.cpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTagEfficiency()              -- Parameterized Constructor            *
 *    ~HistogrammingTagEfficiency()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTagEfficiency"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *    myTH1F* _hNPvtx                   -- Hist of MET                          *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *      20 Nov 2006 - Modified by Bernard Pope                                *
 *      21 Mar 2007 - RS: Fill from event container, add sumET, mex, mey      *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingTagEfficiency.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingTagEfficiency::HistogrammingTagEfficiency(EventContainer *obj)                     *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Particle class                                                     *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingTagEfficiency::HistogrammingTagEfficiency(EventContainer *obj)
{
  SetEventContainer(obj);
  _bTagCut = 0.8484;
} //HistogrammingTagEfficiency()

/******************************************************************************
 * HistogrammingTagEfficiency::~HistogrammingTagEfficiency()                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingTagEfficiency::~HistogrammingTagEfficiency()
{
  
} //HistogrammingTagEfficiency

/******************************************************************************
 * void HistogrammingTagEfficiency::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingTagEfficiency::BookHistogram(){

  //Declare by 2D histograms
  _hTagEff_bFlavNum = DeclareTH2F("hTagEff_bFlavNum","Number of b flavour jets passing b-taggin requirements",100,0.,400.,20,0.,2.4);
  _hTagEff_bFlavNum->SetXAxisTitle("JetPt");
  _hTagEff_bFlavNum->SetYAxisTitle("Eta");

  _hTagEff_bFlavDen = DeclareTH2F("hTagEff_bFlavDwn","Number of b flavour jets",100,0.,400.,20,0.,2.4);
  _hTagEff_bFlavDen->SetXAxisTitle("JetPt");
  _hTagEff_bFlavDen->SetYAxisTitle("Eta");

  //Declare by 2D histograms
  _hTagEff_cFlavNum = DeclareTH2F("hTagEff_cFlavNum","Number of c flavour jets passing b-taggin requirements",100,0.,400.,20,0.,2.4);
  _hTagEff_cFlavNum->SetXAxisTitle("JetPt");
  _hTagEff_cFlavNum->SetYAxisTitle("Eta");

  _hTagEff_cFlavDen = DeclareTH2F("hTagEff_cFlavDwn","Number of c flavour jets",100,0.,400.,20,0.,2.4);
  _hTagEff_cFlavDen->SetXAxisTitle("JetPt");
  _hTagEff_cFlavDen->SetYAxisTitle("Eta");

  //Declare by 2D histograms
  _hTagEff_lightFlavNum = DeclareTH2F("hTagEff_lightFlavNum","Number of light flavour jets passing b-taggin requirements",100,0.,400.,20,0.,2.4);
  _hTagEff_lightFlavNum->SetXAxisTitle("JetPt");
  _hTagEff_lightFlavNum->SetYAxisTitle("Eta");

  _hTagEff_lightFlavDen = DeclareTH2F("hTagEff_lightFlavDwn","Number of light flavour jets",100,0.,400.,20,0.,2.4);
  _hTagEff_lightFlavDen->SetXAxisTitle("JetPt");
  _hTagEff_lightFlavDen->SetYAxisTitle("Eta");

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingTagEfficiency::Apply()                                         *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingTagEfficiency::Apply()
{
  //cout<<"Begin of HistogrammingTagEfficiency::Apply()"<<endl;
  
  EventContainer *evc = GetEventContainer();

  //std::cout << evc->alljets.size() << " " << evc->taggedJets.size() << std::endl;
  
  //Loop over the jet collection and fill in the efficiency histograms as required
  for (auto const & jet : evc->alljets){
    if (fabs(jet.Eta()) > 2.4) continue;
    //    std::cout << jet.hadronFlavour() << " " << jet.GetbDiscriminator() << " " << jet.IsTagged() << " " << jet.tagged() << std::endl;
    switch (jet.hadronFlavour()){
    case 5: //b flavour quark
      _hTagEff_bFlavDen->FillWithoutWeight(jet.Pt(),jet.Eta());
      if (jet.GetbDiscriminator() > _bTagCut) _hTagEff_bFlavNum->FillWithoutWeight(jet.Pt(),jet.Eta());
      break;
    case 4:
      _hTagEff_cFlavDen->FillWithoutWeight(jet.Pt(),jet.Eta());
      if (jet.GetbDiscriminator() > _bTagCut) _hTagEff_cFlavNum->FillWithoutWeight(jet.Pt(),jet.Eta());
      break;
    default:
      _hTagEff_lightFlavDen->FillWithoutWeight(jet.Pt(),jet.Eta());
      if (jet.GetbDiscriminator() > _bTagCut) _hTagEff_lightFlavNum->FillWithoutWeight(jet.Pt(),jet.Eta());
      
    }
  }

  //cout<<"End of HistogrammingTagEfficiency::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












