/******************************************************************************
 * HistogrammingJetAngular.cpp                                                       *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingJetAngular()              -- Parameterized Constructor            *
 *    ~HistogrammingJetAngular()             -- Destructor                           *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingJetAngular"         *
 *                                                                            *
 * Private Data Members of this class                                         *
 *    myTH1F* _hMtW                   -- Hist of MtW                          *
 *                                                                            *
 * History                                                                    *
 *      08 Jul 2016 - Created by D. Leggat for CMS IHEP                       *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Histogramming/Recon/HistogrammingJetAngular.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * HistogrammingJetAngular::HistogrammingJetAngular(EventContainer *obj)                     *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Particle class                                                     *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingJetAngular::HistogrammingJetAngular(EventContainer *obj, bool unisolated)
{
  SetEventContainer(obj);
  _unisolated = unisolated;
} //HistogrammingJetAngular()

/******************************************************************************
 * HistogrammingJetAngular::~HistogrammingJetAngular()                                  *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
HistogrammingJetAngular::~HistogrammingJetAngular()
{
  
} //HistogrammingJetAngular

/******************************************************************************
 * void HistogrammingJetAngular::BookHistogram()                                          *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void HistogrammingJetAngular::BookHistogram(){
  
  //Histogram of delta R between untagged jets
  _hUntaggedJetDelR = DeclareTH1F("untaggedJetsDelR","#Delta R untagged jets",200,0.,5.);
  _hUntaggedJetDelR->SetXAxisTitle("#Delta R_{j,j}");
  _hUntaggedJetDelR->SetYAxisTitle("Events");
  
  //Histogram of delta R between leading untagged jet and b jet
  _hJet1bJetDelR = DeclareTH1F("bJet1DelR","#Delta R leading untagged jet b-jet",200,0.,5.);
  _hJet1bJetDelR->SetXAxisTitle("#Delta R_{b,j_{1}}");
  _hJet1bJetDelR->SetYAxisTitle("Events");

  //Histogram of delta R between subleading untagged jet and b jet
  _hJet2bJetDelR = DeclareTH1F("bJet2DelR","#Delta R subleading untagged jet b-jet",200,0.,5.);
  _hJet2bJetDelR->SetXAxisTitle("#Delta R_{b,j_{2}}");
  _hJet2bJetDelR->SetYAxisTitle("Events");

  //Histogram of delta R between lepton and b jet
  _hlepJetsDelR = DeclareTH1F("lepJetsDelR","lepton jets delta R",200,0.,5.);
  _hlepJetsDelR->SetXAxisTitle("#Delta R_{l,jj}");
  _hlepJetsDelR->SetYAxisTitle("Events");

  //Histogram of delta R between lepton and b jet
  _hlepbJetDelR = DeclareTH1F("lepbJetDelR","lepton b jet delta R",200,0.,5.);
  _hlepbJetDelR->SetXAxisTitle("#Delta R_{l,b}");
  _hlepbJetDelR->SetYAxisTitle("Events");

  //Histogram of delta R between lepton and leading jet
  _hJet1lepDelR = DeclareTH1F("Jet1lepDelR","lepton leading jet delta R",200,0.,5.);
  _hJet1lepDelR->SetXAxisTitle("#Delta R_{l,j_{1}}");
  _hJet1lepDelR->SetYAxisTitle("Events");

  //Histogram of delta R between lepton and sub-leading jet
  _hJet2lepDelR = DeclareTH1F("Jet2lepDelR","lepton sub-leading jet delta R",200,0.,5.);
  _hJet2lepDelR->SetXAxisTitle("#Delta R_{l,j_{2}}");
  _hJet2lepDelR->SetYAxisTitle("Events");

  //Histogram of delta phi between untagged jets
  _hUntaggedJetDelPhi = DeclareTH1F("untaggedJetsDelPhi","#Delta #Phi untagged jets",200,-3.2,3.2);
  _hUntaggedJetDelPhi->SetXAxisTitle("#Delta #Phi_{j,j}");
  _hUntaggedJetDelPhi->SetYAxisTitle("Events");
  
  //Histogram of delta phi between leading untagged jet and b jet
  _hJet1bJetDelPhi = DeclareTH1F("bJet1DelPhi","#Delta #Phi leading untagged jet b-jet",200,-3.2,3.2);
  _hJet1bJetDelPhi->SetXAxisTitle("#Delta #Phi_{b,j_{1}}");
  _hJet1bJetDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between subleading untagged jet and b jet
  _hJet2bJetDelPhi = DeclareTH1F("bJet2DelPhi","#Delta #Phi subleading untagged jet b-jet",200,-3.2,3.2);
  _hJet2bJetDelPhi->SetXAxisTitle("#Delta #Phi_{b,j_{2}}");
  _hJet2bJetDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between lepton and b jet
  _hlepJetsDelPhi = DeclareTH1F("lepJetsDelPhi","lepton jets delta phi",200,-3.2,3.2);
  _hlepJetsDelPhi->SetXAxisTitle("#Delta #Phi_{l,jj}");
  _hlepJetsDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between lepton and b jet
  _hlepbJetDelPhi = DeclareTH1F("lepbJetDelPhi","lepton b jet delta phi",200,-3.2,3.2);
  _hlepbJetDelPhi->SetXAxisTitle("#Delta #Phi_{l,b}");
  _hlepbJetDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between lepton and leading jet
  _hJet1lepDelPhi = DeclareTH1F("Jet1lepDelPhi","lepton leading jet delta phi",200,-3.2,3.2);
  _hJet1lepDelPhi->SetXAxisTitle("#Delta #Phi_{l,j_{1}}");
  _hJet1lepDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between lepton and sub-leading jet
  _hJet2lepDelPhi = DeclareTH1F("Jet2lepDelPhi","lepton sub-leading jet delta phi",200,-3.2,3.2);
  _hJet2lepDelPhi->SetXAxisTitle("#Delta #Phi_{l,j_{2}}");
  _hJet2lepDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between MET and leading untagged jet
  _hJet1METDelPhi = DeclareTH1F("Jet1METDelPhi","MET leading untagged jet delta phi",200,-3.2,3.2);
  _hJet1METDelPhi->SetXAxisTitle("#Delta #Phi_{MET,j_{1}}");
  _hJet1METDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between MET and sub-leading untagged jet
  _hJet2METDelPhi = DeclareTH1F("Jet1METDelPhi","MET sub-leading untagged jet delta phi",200,-3.2,3.2);
  _hJet2METDelPhi->SetXAxisTitle("#Delta #Phi_{MET,j_{2}}");
  _hJet2METDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between MET and sub-leading untagged jet
  _hbJetMETDelPhi = DeclareTH1F("Jet1METDelPhi","MET b jet delta phi",200,-3.2,3.2);
  _hbJetMETDelPhi->SetXAxisTitle("#Delta #Phi_{MET,b}");
  _hbJetMETDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between MET and tagged jets
  _hJetsMETDelPhi = DeclareTH1F("JetsMETDelPhi","MET jets delta phi",200,-3.2,3.2);
  _hJetsMETDelPhi->SetXAxisTitle("#Delta #Phi_{MET,jj}");
  _hJetsMETDelPhi->SetYAxisTitle("Events");

  //Histogram of delta phi between MET and lepton
  _hLepMETDelPhi = DeclareTH1F("LepMETDelPhi","MET lep delta phi",200,-3.2,3.2);
  _hLepMETDelPhi->SetXAxisTitle("#Delta #Phi_{MET,l}");
  _hLepMETDelPhi->SetYAxisTitle("Events");

  //Mass of the untagged jet pair
  _diJetMass = DeclareTH1F("diJetMass", "Dijet mass", 200, 0., 200.);
  _diJetMass->SetXAxisTitle("m_{dijet}");
  _diJetMass->SetYAxisTitle("Events");

  //Mass of the three jet system
  _threeJetMass = DeclareTH1F("threeJetMass", "Three jet mass", 200,0.,500.);
  _threeJetMass->SetXAxisTitle("m_{jjj}");
  _threeJetMass->SetYAxisTitle("Events");

  // transverse mass of the reconstructed top quark
  _mTt = DeclareTH1F("mTt", "Transverse top mass", 200,0.,200.);
  _mTt->SetXAxisTitle("m_{T,t}");
  _mTt->SetYAxisTitle("Events");

  // pt of the system
  _ptSys = DeclareTH1F("ptSys","p_{T} of the system", 200, 0., 200.);
  _ptSys->SetXAxisTitle("p_{T,sys}");
  _ptSys->SetYAxisTitle("Events");

  //cout<<"end of HistogrammingJetAngular::BookHistogram"<<endl;

} //BookHistogram

/******************************************************************************
 * Bool_t HistogrammingJetAngular::Apply()                                         *
 *                                                                            *
 * Fill histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t HistogrammingJetAngular::Apply()
{
  //cout<<"Begin of HistogrammingJetAngular::Apply()"<<endl;
  
  EventContainer *evc = GetEventContainer();
  
  bJet = evc->taggedJets[0];
  int jetCount = 0;
  leadingJet = evc->unTaggedJets[0];
  subleadingJet = evc->unTaggedJets[1];
  /*  for (auto const jet : evc->jets){
    std::cout << jet.Pt();
    if (jet.IsTagged()) continue;
    std::cout << jet.Pt() << " ";
    if (jetCount == 0) {
      leadingJet = jet;
      jetCount++;
    }
    else subleadingJet = jet;
    std::cout << leadingJet.Pt() << " " << subleadingJet.Pt() << std::endl;
    }*/

  if (!_unisolated){
    if (evc->tightMuons.size() > 0){
      lepton = evc->tightMuons[0];
    }
    else {
      lepton = evc->tightElectrons[0];
    }
  }
  else {
    if (evc->unIsolatedMuons.size() > 0){
      lepton = evc->unIsolatedMuons[0];
    }
    else {
      lepton = evc->unIsolatedElectrons[0];
    }
  }
    
  met.SetPtEtaPhiE(evc->missingEt_xy,0,evc->missingPhi_xy,evc->missingEt_xy);
  /*
  std::cout << "Lead " <<leadingJet.Pt() << " " << lepton.DeltaR(leadingJet) << std::endl;
  std::cout << "Second " << subleadingJet.Pt() << " " << lepton.DeltaR(subleadingJet) << std::endl;
  std::cout << "bjet "  <<bJet.Pt() << " " << lepton.DeltaR(bJet) << std::endl;
  */
  
  _hUntaggedJetDelR   	-> Fill(leadingJet.DeltaR(subleadingJet));
  _hJet1bJetDelR	-> Fill(leadingJet.DeltaR(bJet));
  _hJet2bJetDelR	-> Fill(subleadingJet.DeltaR(bJet));
  _hlepJetsDelR		-> Fill(lepton.DeltaR(leadingJet+subleadingJet));
  _hlepbJetDelR		-> Fill(lepton.DeltaR(bJet));
  _hJet1lepDelR		-> Fill(lepton.DeltaR(leadingJet));
  _hJet2lepDelR		-> Fill(lepton.DeltaR(subleadingJet));

  _hUntaggedJetDelPhi   -> Fill(leadingJet.DeltaPhi(subleadingJet));       
  _hJet1bJetDelPhi      -> Fill(leadingJet.DeltaPhi(bJet));                
  _hJet2bJetDelPhi      -> Fill(subleadingJet.DeltaPhi(bJet));             
  _hlepJetsDelPhi       -> Fill(lepton.DeltaPhi(leadingJet+subleadingJet));
  _hlepbJetDelPhi       -> Fill(lepton.DeltaPhi(bJet));                    
  _hJet1lepDelPhi       -> Fill(lepton.DeltaPhi(leadingJet));              
  _hJet2lepDelPhi       -> Fill(lepton.DeltaPhi(subleadingJet));           
  _hJet1METDelPhi	-> Fill(met.DeltaPhi(leadingJet));
  _hJet2METDelPhi	-> Fill(met.DeltaPhi(subleadingJet));
  _hbJetMETDelPhi	-> Fill(met.DeltaPhi(bJet));
  _hJetsMETDelPhi	-> Fill(met.DeltaPhi(leadingJet + subleadingJet));
  _hLepMETDelPhi	-> Fill(met.DeltaPhi(lepton));

  _diJetMass		-> Fill((subleadingJet + leadingJet).M());
  _threeJetMass		-> Fill((subleadingJet + leadingJet + bJet).M());
  _mTt			-> Fill(std::sqrt(2*met.Pt()*(lepton + bJet).Pt()*(1-cos(met.Phi() - (lepton+bJet).Phi()))));
  
  float totPx = lepton.Px() + met.Px() + leadingJet.Py() + subleadingJet.Px() + bJet.Px();
  float totPy = lepton.Py() + met.Py() + leadingJet.Py() + subleadingJet.Py() + bJet.Py();
  
  _ptSys		-> Fill(std::sqrt(totPx * totPx + totPy * totPy));
  
  
  //  }
  //cout<<"End of HistogrammingJetAngular::Apply()"<<endl;
  return kTRUE;  
  
} //Apply












