/******************************************************************************
 * ChannelFlag.hpp                                                                *
 *                                                                            *
 * Add various weight variables to the skim tree
 * 
 * History
 *      10 Nov 2016 - Created by D. Leggat
 ******************************************************************************/

#ifndef ChannelFlag_h
#define ChannelFlag_h

#include "SingleTopRootAnalysis/Base/Dictionary/VarBase.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp" 

class ChannelFlag: public VarBase {
  
public:
  
  enum Channel {
    j3t1=0,
    j3t2=1,
    j2t1=2,
    j4t1=3,
    j4t2=4,
  };

  ChannelFlag();
  
  //void BookBranches(TTree * skimTree);
  void FillBranches(EventContainer * evtObj);

private:
  
};

#endif
