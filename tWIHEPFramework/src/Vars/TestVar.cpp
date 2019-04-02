/******************************************************************************
 * testVar.cpp                                                                *
 *
 * Testing out some variables being added to the skim tree.
 *
 * History
 *      10th Nov 2016 - Created by D. Leggat. Just a test, so probably not that important
 ******************************************************************************/

#include "SingleTopRootAnalysis/Vars/TestVar.hpp"

#include <iostream>

//Test out a couple of variables, one int and one float I guess
TestVar::TestVar(){
  _floatVars["testFloat"] = 0.;
  _intVars["testInt"] = 0;
}

void TestVar::FillBranches(EventContainer * evtObj){
  _floatVars["testFloat"] = evtObj->taggedJets[0].Pt();
  _intVars["testInt"] = evtObj->nPvtx;
}
