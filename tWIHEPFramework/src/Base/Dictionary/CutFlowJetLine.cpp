/******************************************************************************
 * CutFlowJetLine.hpp                                                         *
 *                                                                            *
 * Holds the pass/fail info for a single cut.                                 *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of this class                                      *
 *    CutFlowJetLine()      -- Parameterized Constructor                      *
 *    ~CutFlowJetLine()     -- Destructor                                     *
 *    Clear()            -- Reset internal counters                           *
 *    Print()            -- Print out cut statistics                          *
 *    PassCut()          -- Tell this class about a passed event              *
 *    FailCut()          -- Tell this class about a failed event              *
 *    GetExposedCount()  -- Returns fail+pass events                          *
 *    GetExposedYield()  -- Returns sum of weights for fail+pass events       *
 *    GetPassCount()  -- Returns passed events                                *
 *    GetPassYield()  -- Returns sum of weights for passed events             *
 *                                                                            *
 * Private Data Members of this class                                         *
 * Counts and yields for exposed and passed events                            *
 *                                                                            *
 * History                                                                    *
 *      14 Sept 2011 - Created by James Koll based off CutFlowLine            *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Base/Dictionary/CutFlowJetLine.hpp"

#include <iostream>
#include <iomanip>

using namespace std;

//
// CutFlowJetLine
//
//  Default constructor
//
CutFlowJetLine::CutFlowJetLine(void)
  : _exposedCount(0), _exposedYield(0.)
  , TNamed ("NONE","")
{
    for (int i=0;i!=10;++i){
      _passCount[i]=0; 
      _passYield[i]=0.0;
    }
}

// Normal constructor
CutFlowJetLine::CutFlowJetLine(const TString& name, const TString& title)
  : _exposedCount(0), _exposedYield(0.)
  , TNamed (name.Data(), title.Data())
{
    for (int i=0;i!=10;++i){
      _passCount[i]=0; 
      _passYield[i]=0.0;
    }
}
void CutFlowJetLine::PassCut(const Double_t& weight, const Int_t NJets) { 
  _exposedCount++;
  _exposedYield+=weight; 
  for (int jets=0;jets!=8;++jets){
    if(NJets==jets){
      _passCount[jets]++; 
      _passYield[jets]+=weight; 
    }
  }  
  if (NJets>=8){
  _passCount[8]++; 
  _passYield[8]+=weight;     
  }
  _passCount[9]++; 
  _passYield[9]+=weight; 
};
 void CutFlowJetLine::FailCut(const Double_t& weight) { _exposedCount++; _exposedYield+=weight; };
//
// Print
//
void CutFlowJetLine::Print (Option_t *option) const
{
  cout<<"|  " << setw(30) << GetName() << " |";
  for (int jets=0;jets!=10;++jets){
    cout << setw(7) << setprecision(6)<<_passCount[jets]<<"|";
  }
  cout << endl;
}

