/******************************************************************************
 * testVar.cpp                                                                *
 *
 * Testing out some variables being added to the skim tree.
 *
 * History
 *      10th Nov 2016 - Created by D. Leggat. Just a test, so probably not that important
 ******************************************************************************/

#include "SingleTopRootAnalysis/Vars/ChannelFlag.hpp"

#include <iostream>

//Test out a couple of variables, one int and one float I guess
ChannelFlag::ChannelFlag(){

  _intVars["channel"] = 0;
  
}

void ChannelFlag::FillBranches(EventContainer * evtObj){

  //Get the number of jets and tagged jets to determine value of channel.
  int nJets = evtObj->jets.size();
  int nbJets = evtObj->taggedJets.size();

  if (nJets == 2 && nbJets == 1) _intVars["channel"] = ChannelFlag::j2t1;
  if (nJets == 3 && nbJets == 1) _intVars["channel"] = ChannelFlag::j3t1;
  if (nJets == 4 && nbJets == 1) _intVars["channel"] = ChannelFlag::j4t1;
  if (nJets == 3 && nbJets == 2) _intVars["channel"] = ChannelFlag::j3t2;
  if (nJets == 4 && nbJets == 2) _intVars["channel"] = ChannelFlag::j4t2;
    

}
