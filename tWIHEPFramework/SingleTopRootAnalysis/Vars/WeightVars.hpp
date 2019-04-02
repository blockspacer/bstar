/******************************************************************************
 * WeightVars.hpp                                                                *
 *                                                                            *
 * Add various weight variables to the skim tree
 * 
 * History
 *      10 Nov 2016 - Created by D. Leggat
 ******************************************************************************/

#ifndef WeightVars_h
#define WeightVars_h

#include "SingleTopRootAnalysis/Base/Dictionary/VarBase.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp" 

class WeightVars: public VarBase {
  
public:
  WeightVars();
  
  //void BookBranches(TTree * skimTree);
  void FillBranches(EventContainer * evtObj);

private:
  std::vector<std::string> _bTagSystNames;
  
};

#endif
