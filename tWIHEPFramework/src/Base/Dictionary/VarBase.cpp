/******************************************************************************
 * VarBase.cpp
 *
 *
 * Abstract class
 *
 * The definition of a class that will add variables to the skim tree
 * Works in a similar manner to BaseCut, i.e. doesn't really do anything
 * itself, but allows variables to be booked and saved later on.
 * 
 * Something about public member functions:
 *
 *
 * History
 *      10 Nov 2016 - Created by D Leggat
 *****************************************************************************/

#include "SingleTopRootAnalysis/Base/Dictionary/VarBase.hpp"

#include "TROOT.h"

// Integrate classes into the Root system
ClassImp(VarBase);


/******************************************************************************  
 * void VarBase::VarBase ()		                                      *  
 *                                                                            *  
 * Default constructor                                                        *  
 *                                                                            *  
 * Input:  None                                                               *  
 * Output: None                                                               *  
 ******************************************************************************/ 
VarBase::VarBase(bool addHists){
  //Set up whether we want to make histograms of the BDT variables
  _makeHists = addHists;
//Nothing immediately happens in here.
}



/****************************************************************************** 
 * void VarBase::SetName (const char *name)                                   * 
 *                                                                            * 
 * Set Name of Object                                                         * 
 *                                                                            * 
 * Input:  name to set                                                        * 
 * Output: None                                                               * 
 ******************************************************************************/
void VarBase::SetName(const char *name){
  //Set name
  TNamed::SetName(name);
  
} // end SetName

/****************************************************************************** 
 * void VarBase::BookBranches (TTree * skimTree)                              * 
 *                                                                            * 
 * Book the branches into skim tree                                           * 
 *                                                                            * 
 * Input:  The skimTree is the TTree to associate the variables with          *
 * Output: None                                                               * 
 ******************************************************************************/
void VarBase::BookBranches(TTree * skimTree){

  for (auto intVar : _intVars){
    string tempString = intVar.first;
    _branchVec[tempString.c_str()] = skimTree->Branch(tempString.c_str(),&(_intVars[tempString.c_str()]),(tempString+"/I").c_str());
    if (DoHists()) {
      _histograms[tempString] = BookTH1FHistogram(tempString,tempString,10,0.,intVar.second);
    }
  }

  for (auto floatVar: _floatVars){
    string tempString = floatVar.first;
    _branchVec[tempString.c_str()] = skimTree->Branch(tempString.c_str(),&(_floatVars[tempString.c_str()]),(tempString+"/F").c_str());
    if (DoHists()) {
      float startVar = 0.;
      if (floatVar.second < 0.) startVar = floatVar.second;
      _histograms[tempString] = BookTH1FHistogram(tempString,tempString,100,startVar,std::fabs(floatVar.second));
    }
  }


}

/****************************************************************************** 
 * void VarBase::ResetBranches ()                                             * 
 *                                                                            * 
 * Reset the branches in skim tree                                            * 
 *                                                                            * 
 * Input:  name to set                                                        * 
 * Output: None                                                               * 
 ******************************************************************************/
void VarBase::ResetBranches(){
  for (auto intVar : _intVars){
    intVar.second = -999;
  }
  for (auto floatVar : _floatVars){
    floatVar.second = -999.;
  }
}

/****************************************************************************** 
 * void VarBase::OutputBranches ()                                             * 
 *                                                                            * 
 * Output the value of the variables
 *                                                                            * 
 * Input:  None                                                               * 
 * Output: None                                                               * 
 ******************************************************************************/
void VarBase::OutputBranches(){
  for (auto floatVar : _floatVars){
    std::cout << floatVar.first << " " << floatVar.second << std::endl;
  }

  for (auto intVar : _intVars){
    std::cout << intVar.first << " " << intVar.second << std::endl;
  }
}

/****************************************************************************** 
 * void VarBase::TouchBranches ()                                             * 
 *                                                                            * 
 * Touch every filled variable because that helped for some reason
 *                                                                            * 
 * Input:  None                                                               * 
 * Output: None                                                               * 
 ******************************************************************************/
void VarBase::TouchBranches(){
  for (auto floatVar : _floatVars){
    float temp = floatVar.second;
  }
  for (auto intVar : _intVars){
    int temp = intVar.second;
  }
}

/****************************************************************************** 
 * void VarBase::FillHistograms ()                                             * 
 *                                                                            * 
 * Fill the histograms we're making here *
 *                                                                            * 
 * Input:  None                                                               * 
 * Output: None                                                               * 
 ******************************************************************************/
void VarBase::FillHistograms(double weight){
  double fillWeight = weight;
  if (weight == -999.) weight = 1.;
  for (auto intVar : _intVars){
    //    std::cout << intVar.first << " " << intVar.second << std::endl;
    _histograms[intVar.first]->Fill(intVar.second,weight);
  }
  
  for (auto floatVar : _floatVars){
    //std::cout << floatVar.first << " " << floatVar.second << std::endl;
    _histograms[floatVar.first]->Fill(floatVar.second,weight);
  }
}
