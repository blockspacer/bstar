/******************************************************************************
 * CutRunNumber.hpp                                                              *
 *                                                                            *
 * Cuts on charge and flavor of lepton                                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutRunNumber class                                 *
 *    CutRunNumber()                       -- Parameterized Constructor          *
 *    ~CutRunNumber()                      -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutRunNumber"                *
 *                                                                            *
 * Private Data Members of CutJetN                                           *
 *    Int_t _leptonCharge                -- lepton charge cut value           *
 *    Int_t _leptonFlavor                -- lepton flavor cut value           *
 *    myTH1F* _hLeptonCharge             -- Hist of lepton charge             *
 *    myTH1F* _hLeptonFlavor             -- Hist of lepton flavor             *
 *                                                                            *
 * History                                                                    *
 *      09 Jun 2010 - Created by J. Holzbauer                                *
 *****************************************************************************/

#ifndef CutRunNumber_h
#define CutRunNumber_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include <fstream>

class CutRunNumber : public HistoCut
{

public:

  // Parameterized Constructor
  CutRunNumber(EventContainer *obj);
  
  // Destructor
  ~CutRunNumber();

  // Book Histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Get the name describing the cut
  inline std::string GetCutName() { return "CutRunNumber"; }


private:

  // Cut parameters
   ifstream * _numberFile;
  int _RNint;
  // Histograms
  myTH1F*   _hRunNumberBefore;  // Lepton charge before cut
  myTH1F*   _hRunNumberAfter;   // Lepton charge after cut

};


#endif







