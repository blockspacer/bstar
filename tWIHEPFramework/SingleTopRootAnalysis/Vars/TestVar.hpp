/******************************************************************************
 * testVar.hpp                                                                *
 *                                                                            *
 * Puts some testing variables into the skim tree
 * 
 * History
 *      10 Nov 2016 - Created by D. Leggat
 ******************************************************************************/

#ifndef TestVar_h
#define TestVar_h

#include "SingleTopRootAnalysis/Base/Dictionary/VarBase.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp" 

class TestVar: public VarBase {
  
public:
  TestVar();
  
  //void BookBranches(TTree * skimTree);
  void FillBranches(EventContainer * evtObj);
  
};

#endif
