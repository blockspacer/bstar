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
 *      17 June 2015 - Huaqiao ZHANG Port from STAR                           *
 *****************************************************************************/

#ifndef __CutFlowLine__
#define __CutFlowLine__

#include "TNamed.h"
#include <iostream>
#include <string>

class CutFlowLine : public TNamed {
public:
  CutFlowLine(void);
  CutFlowLine(const TString &name, const TString &title="");

  // reset internal counters
  inline void Clear() { _exposedCount=0; _exposedYield=0.; _passCount=0; _passYield=0.; };

  // print out this line
  void Print(Option_t *option = "") const; 
  
  // call this to tell the line that this cut passed or failed.
  inline void PassCut(const Double_t& weight) { _exposedCount++; _exposedYield+=weight; _passCount++; _passYield+=weight;
  //std::cout<<"weight is :"<<weight<<std::endl;
  };
  inline void FailCut(const Double_t& weight) { _exposedCount++; _exposedYield+=weight; };

  // get the pass count and yield.
  inline const Int_t GetExposedCount(void) const { return _exposedCount; };
  inline const Double_t GetExposedYield(void) const { return _exposedYield; };
  inline const Int_t GetPassCount(void) const { return _passCount; };
  inline const Double_t GetPassYield(void) const { return _passYield; };
  inline const Int_t GetFailCount(void) const { return(_exposedCount - _passCount); };
  inline const Double_t GetFailYield(void) const { return(_exposedYield - _passYield); };

private:
  Int_t _exposedCount;  // number of events exposed to this cut
  Double_t _exposedYield;  // sum of weights for events exposed to this cut
  Int_t _passCount;  // number of events passing this cut
  Double_t _passYield;  // sum of weights for events passing this cut

  ClassDef (CutFlowLine, 0);
};




#endif
