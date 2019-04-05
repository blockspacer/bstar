/******************************************************************************
 * CutFlowLine.hpp                                                            *
 *                                                                            *
 * Holds the pass/fail info for a single cut.                                 *
 *                                                                            *
 *                                                                            *
 * Public Member Functions of this class                                      *
 *    CutFlowLine()      -- Parameterized Constructor                         *
 *    ~CutFlowLine()     -- Destructor                                        *
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
 *      21 Nov 2006 - Created by R. Schwienhorst for ATLAS                    *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Base/Dictionary/CutFlowLine.hpp"

#include <iostream>
#include <iomanip>

using namespace std;

//
// CutFlowLine
//
//  Default constructor
//
CutFlowLine::CutFlowLine(void)
  : _exposedCount(0), _exposedYield(0.), _passCount(0), _passYield(0.)
  , TNamed ("NONE","")
{
}

// Normal constructor
CutFlowLine::CutFlowLine(const TString& name, const TString& title)
  : _exposedCount(0), _exposedYield(0.), _passCount(0), _passYield(0.)
  , TNamed (name.Data(), title.Data())
{
}

//
// Print
//
void CutFlowLine::Print (Option_t *option) const
{
  cout<<"|  " << setw(29) << GetName() << " | ";
	cout << setw(7) <<_exposedCount;
	cout << " | "<<setw(12) <<_exposedYield;
	cout << " | ";
	cout << setw(7) <<_passCount;
	cout << " | " <<setw(12) <<_passYield;
	cout << " | " <<GetTitle()<<endl;
	//cout<<" |"<<endl;
}

