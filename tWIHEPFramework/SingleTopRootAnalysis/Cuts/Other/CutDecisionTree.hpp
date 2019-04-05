/******************************************************************************
 * CutDecisionTree.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutDecisionTree class                              *
 *    CutDecisionTree()                    -- Parameterized Constructor          *
 *    ~CutDecisionTree()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutDecisionTree"             *
 *                                                                            *
 * Private Data Members of CutDecisionTree                                       *
 *    myTH1F* _hDecisionTreeBefore         -- Hist of DecisionTree before cuts      *
 *    myTH1F* _hDecisionTreeAfter          -- Hist of DecisionTree after cuts       *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *    Int_t _missingEtMin;              -- Minimum Missing Et                 *
 *                                                                            *
 * History                                                                    *
 *      15 Dec 2006 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef CutDecisionTree_h
#define CutDecisionTree_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutDecisionTree : public HistoCut 
{

public:

  // Parameterized Constructor
  CutDecisionTree(EventContainer *EventContainerObj, TString branchname = "None", Bool_t cutTrigger = kFALSE);
  
  // Destructor
  ~CutDecisionTree();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutDecisionTree"; };

  void SetDecisionBranch(TString str) { _DecisionBranch=str; };
  TString DecisionBranch() const { return _DecisionBranch; };
  TString GetDecisionBranch() const { return _DecisionBranch; };

  void SetCutTrigger(Bool_t ct) { _CutTrigger=ct; };
  Bool_t CutTrigger() const { return _CutTrigger; };
  Bool_t GetCutTrigger() const { return _CutTrigger; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hDecisionTreeBefore;    // Histogram Missing Et
  myTH1F* _hDecisionTreeAfter;     // Histogram Missing Et

  TString _DecisionBranch;
  Bool_t _CutTrigger;

};


#endif







