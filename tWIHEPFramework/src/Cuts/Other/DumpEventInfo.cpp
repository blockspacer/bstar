/******************************************************************************
 * DumpEventInfo.cpp                                                       *
 *                                                                            *
 * Cuts on the requirement of a good primary vertex                           *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of DumpEventInfo class                              *
 *    DumpEventInfo()                     -- Parameterized Constructor         *
 *    ~DumpEventInfo()                    -- Destructor                        *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "DumpEventInfo"                *
 *                                                                            *
 * Private Data Members of DumpEventInfo class                                 *
 *    myTH1F* _hPrimaryVertexBefore;    -- Hist of PV before cut        *
 *    myTH1F* _hPrimaryVertexAfter;     -- Hist of PV after cut         *
 *                                                                            *
 * History                                                                    *
 *      5th July 2016 - Created by Duncan Leggat                              *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Cuts/Other/DumpEventInfo.hpp"
#include<iostream>

using namespace std;

/******************************************************************************
 * DumpEventInfo::DumpEventInfo(EventContainer *EventContainerObj, TString electronType)    *
 *                                                                            *
 * Parameterized Constructor                                                  *
 *                                                                            *
 * Input:  Event Object class                                                 *
 * Output: None                                                               *
 ******************************************************************************/
DumpEventInfo::DumpEventInfo(EventContainer *EventContainerObj)
{
  // Set Event Container
  SetEventContainer(EventContainerObj);
} // DumpEventInfo


/******************************************************************************
 * DumpEventInfo::~DumpEventInfo()                                              *
 *                                                                            *
 * Destructor                                                                 *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
DumpEventInfo::~DumpEventInfo()
{
  
}//~DumpEventInfo

/******************************************************************************
 * void DumpEventInfo::BookHistogram()                                         *
 *                                                                            *
 * Book Histograms                                                            *
 *                                                                            *
 * Input:  None                                                               *
 * Output: None                                                               *
 ******************************************************************************/
void DumpEventInfo::BookHistogram(){
  


}//BookHistograms()

/******************************************************************************
 * Bool_t DumpEventInfo::Apply()                                               *
 *                                                                            *
 * Apply cuts and fill histograms                                             *
 *                                                                            *
 * Cut on the selected trigger                                                *
 *                                                                            *
 * Input:  None                                                               *
 * Output: kTrue if successful                                                *
 ******************************************************************************/
Bool_t DumpEventInfo::Apply()
{

  //Get event container
  EventContainer *EventContainerObj = GetEventContainer();

  EventTree* evtTree = EventContainerObj->GetEventTree();
  Bool_t passMet = (evtTree->Flag_goodVertices && evtTree->Flag_HBHENoiseFilter && evtTree->Flag_HBHENoiseIsoFilter && evtTree->Flag_EcalDeadCellTriggerPrimitiveFilter);

  //Doesit pass the MET filters?
  std::cout << "MET filters: " << EventContainerObj->passesMETFilters << " MET filter no halo: " << passMet << std::endl;
  //print leptons
  std::cout << "Tight electrons:" << std::endl;
  Particle leadingLep;
  leadingLep.SetPxPyPzE(0.,0.,0.,0.);
  Particle subLeadingLep;
  subLeadingLep.SetPxPyPzE(0.,0.,0.,0.);

  for (auto ele : EventContainerObj->tightElectrons){
    std::cout << "Pt: " << ele.Pt() << " eta: " << ele.Eta() << " phi: " << ele.Phi() << " Id: " << (ele.passTightId()?"Pass":"Fail") << std::endl;
    if (ele.Pt() > leadingLep.Pt()) leadingLep = ele;
    else if (ele.Pt() > subLeadingLep.Pt()) subLeadingLep = ele;
  }
  std::cout << "Tight muons:" << std::endl;
  for (auto mu : EventContainerObj->tightMuons){
    std::cout << "Pt: " << mu.Pt() << " eta: " << mu.Eta() << " phi: " << mu.Phi() << std::endl;
    if (mu.Pt() > leadingLep.Pt()) leadingLep = mu;
    else if (mu.Pt() > subLeadingLep.Pt()) subLeadingLep = mu;
  }
  std::cout << "leading Pt : " << leadingLep.Pt() << " sub: " << subLeadingLep.Pt() << " mll: " << (leadingLep+subLeadingLep).M() << std::endl; 
  std::cout << "MET type 1: " << EventContainerObj->missingEt << " (xy corrected: " << EventContainerObj->missingEt_xy << ")" << std::endl;
  std::cout << "Selected jets:" << std::endl;
  for (auto jet : EventContainerObj->jets){
    std::cout << "Pt: " << jet.Pt() << " eta: " << jet.Eta() << " phi: " << jet.Phi() << " CSVv2: " << jet.bDiscriminator() << " closest lepton: " << jet.closestLep() << std::endl;
  }

  std::cout << "Electrons to use size: " << EventContainerObj->electronsToUsePtr->size() << std::endl;
  std::cout << "Muons to use size: " <<EventContainerObj->muonsToUsePtr->size() << std::endl;
  /*
  std::vector<Particle> particleVec;
  std::map<int,std::vector<std::string> > additionalInfo;
  if (_particle == "electron"){
    if (_particleType == "Tight"){
      std::cout << "Tight electrons:" <<std::endl;
      for (int it = 0; it < EventContainerObj->tightElectrons.size(); it++) {
	particleVec.push_back(EventContainerObj->tightElectrons[it]);
	std::string tempString = EventContainerObj->tightElectrons[it].passTightId()?"Pass":"Fail";
	additionalInfo[it].push_back(" Id: " + tempString);
      }
      //      particleVec = EventContainerObj->tightElectrons;
    }
    if (_particleType == "All"){
      std::cout << "All electrons:" << std::endl;
      for (int it = 0; it < EventContainerObj->electrons.size(); it++) {
	particleVec.push_back(EventContainerObj->electrons[it]);
	std::string tempString = EventContainerObj->electrons[it].passTightId()?"Pass":"Fail";
	additionalInfo[it].push_back(" Id: " + tempString);
      }
    }
  }
  if (_particle == "muon"){
    if (_particleType == "Tight"){
      std::cout << "Tight muons:" <<std::endl;
      for (auto ele : EventContainerObj->tightMuons) particleVec.push_back(ele);
      //  particleVec = EventContainerObj->tightMuons;
    }
  }
  if (_particle == "jet"){
    if (_particleType == "Tight"){
      std::cout << "Tight jets:" << std::endl;
      for (auto jet : EventContainerObj->jets) particleVec.push_back(jet);
    }
  }
  
  //  for (auto const particle : particleVec){
  for (int i= 0 ; i < particleVec.size(); i++){
    auto particle = particleVec[i];
    std::cout << "Pt: " << particle.Pt() << " eta: " << particle.Eta();
    for (auto stringStr : additionalInfo[i]){
      std::cout << stringStr;
    }
    std::cout << std::endl;
    }*/
  
  return true;

} //Apply












