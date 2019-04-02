/******************************************************************************
 * bstarVars.hpp                                                                *
 *                                                                            *
 * Add the variables we will use in the BDT to the skim tree                  *
 * 
 * History
 *      25 Nov 2016 - Created by D. Leggat
 ******************************************************************************/

#ifndef bstarVars_h
#define bstarVars_h

#include "SingleTopRootAnalysis/Base/Dictionary/VarBase.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp" 
#include "TLorentzVector.h"

class bstarVars: public VarBase {
  
public:
  bstarVars(bool makeHistos = false);
  
  //void BookBranches(TTree * skimTree);
  void FillBranches(EventContainer * evtObj);
  
private:
  std::vector<TLorentzVector> Jet;

  std::vector<TLorentzVector> BJet;


  TLorentzVector sumJet;
  TLorentzVector sumBJet;
  TLorentzVector totalJets;

};

#endif
