/******************************************************************************
 * CutFlowTable.cpp                                                           *
 *                                                                            *
 * Keep track of events passing/failing cuts.                                 *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of this class                                      *
 *    CutFlowTable()      -- Parameterized Constructor                        *
 *    ~CutFlowTable()     -- Destructor                                       *
 *    Clear()            -- Reset internal counters                           *
 *    PrintCutFlow()     -- Print out cut statistics                          *
 *    AddCutToFlow()     -- Add a new cut to the table                        *
 *    PassCut()          -- Tell this class about a cut that passed           *
 *    FailCut()          -- Tell this class about a cut that failed           *
 *    FindCutName()      -- Returns true if this cut already exists           *
 *    SetEventContainer()-- Set the internal eventContainer pointer.          *
 *                                                                            *
 * Private Data Members of this class                                         *
 * vector<pair<TString,CutFlowLine> >: list of cuts                           *
 * EventContainer*: pointer to event container                                *
 *                                                                            *
 * History                                                                    *
 *      21 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *       5 Oct 2007 - Changed PassCut and FailCut for MCatNLO weight          *
 *****************************************************************************/
#include "SingleTopRootAnalysis/Base/Dictionary/CutFlowTable.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <iostream>
#include <iomanip>

using namespace std;

//
// default constructor
//  First one is most likely used only by root during readback.
//
CutFlowTable::CutFlowTable(void)
  : TNamed ("NONE",""), _EventObj(NULL)
{
  _cutLines.clear();
  _cutJetLines.clear();
}

CutFlowTable::CutFlowTable(const TString& name, const TString& title)
  : TNamed (name.Data(), title.Data()), _EventObj(NULL)
{
  _cutLines.clear();
  _cutJetLines.clear();
}

//
// Print
//
void CutFlowTable::PrintCutFlow(Option_t *option) const
{
   // don't print if there are no cuts 
  if(_cutLines.size()==0) return;
  // print
  cout << "<CutFlowJetTable> " << GetName() << " result:"<<endl;
  cout << "-----------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << "|    Cut Name                     | 0 jet | 1 jet | 2 jet | 3 jet | 4 jet | 5 jet | 6 jet | 7 jet | 8+jet | total |" << endl; 
  for(int iCut=0; iCut<_cutLines.size(); iCut++) {
    _cutJetLines[iCut].second.Print(option);
  }
  cout << "-----------------------------------------------------------------------------------------------------------------------------" << endl;

  // print
  cout << "<CutFlowTable> " << GetName() << " result:"<<endl;
  cout << "-----------------------------------------------------------------------------------------------------------------------------" << endl;
  cout << "|    Cut Name                    |        exposed         |        pass            |      cut" << endl; 
  cout << "|                                |  count  |     yield    |  count  |    yield     |" << endl; 
  for(int iCut=0; iCut<_cutLines.size(); iCut++) {
    _cutLines[iCut].second.Print(option);
  }
  cout << "-----------------------------------------------------------------------------------------------------------------------------" << endl;
   
}

//
// AddCut
//
//  Add a new cut to the list of cuts
//
void CutFlowTable::AddCutToFlow(const TString& name, const TString& title)
{

  // first check if this cut is present
  if(FindCutName(name)>=0) {
    cout<<"<CutFlowTable::AddCutToFlow> "<<"cut "<<name<<" is already in the list!"<<endl;    
  }
  else {
    CutFlowLine newline(name,title);
    CutFlowJetLine newjetline(name,title);
    _cutLines.push_back(make_pair(name,newline));
    _cutJetLines.push_back(make_pair(name,newjetline));
  }
}

//
// FailCut
//
//  This cut failed, and we are not given a weight, use the global weight.
//
void CutFlowTable::FailCut(const TString& name)
{
  Double_t wgt=1.;
  if(_EventObj) wgt = _EventObj->GetOutputEventWeight(); // takes into account MC@NLO
  FailCut(name,wgt);
}

//
// FailCut
//
//  This cut failed, and we are given a weight. 
//
void CutFlowTable::FailCut(const TString& name, const Double_t& weight)
{

  // first check if this cut is present
  if(FindCutName(name)<0) {
    cout<<"<CutFlowTable::FailCut> "<<"cut "<<name<<" not found! Adding it to the list"<<endl;    
    CutFlowLine newline(name);
    CutFlowJetLine newjetline(name);
    _cutLines.push_back(make_pair(name,newline));
    _cutJetLines.push_back(make_pair(name,newjetline));
  }
  _cutLines[FindCutName(name)].second.FailCut(weight);
  _cutJetLines[FindCutName(name)].second.FailCut(weight);
}

//
// PassCut
//
//  This cut passed, and we are not given a weight, use the global weight.
//
void CutFlowTable::PassCut(const TString& name)
{
  Double_t wgt=1.;
  if(_EventObj) wgt = _EventObj->GetOutputEventWeight(); // takes into account MC@NLO
  PassCut(name,wgt);
}

//
// PassCut
//
//  This cut passed, and we are given a weight. Inc the guy!
//
void CutFlowTable::PassCut(const TString& name, const Double_t& weight)
{

  // first check if this cut is present
  if(FindCutName(name)<0) {
    cout<<"<CutFlowTable::PassCut> "<<"cut "<<name<<" not found! Adding it to the list"<<endl;    
    CutFlowLine newline(name);
    pair<TString,CutFlowLine> newEntry(name,newline);
    CutFlowJetLine newjetline(name);
    pair<TString,CutFlowJetLine> newJetEntry(name,newjetline);
    _cutLines.push_back(newEntry);
    _cutJetLines.push_back(newJetEntry);
  }
  _cutLines[FindCutName(name)].second.PassCut(weight);
  _cutJetLines[FindCutName(name)].second.PassCut(weight,_EventObj->jets.size());
}

//
// FindCutName
//
// Find the index of the cut with this name in the vector.
// Returns -1 if this cut name is not found.
//
Int_t CutFlowTable::FindCutName(const TString& name) const
{
  Int_t res=-1;
  for(int i=0;i<_cutLines.size();i++) {
    if(_cutLines[i].first==name) {
      return i;
    }
  }
  return -1;
}
