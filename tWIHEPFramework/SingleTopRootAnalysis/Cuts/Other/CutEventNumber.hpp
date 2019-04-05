/******************************************************************************
 * CutEventNumber.hpp                                                              *
 *                                                                            *
 * Cuts on charge and flavor of lepton                                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutEventNumber class                                 *
 *    CutEventNumber()                       -- Parameterized Constructor          *
 *    ~CutEventNumber()                      -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutEventNumber"                *
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

#ifndef CutEventNumber_h
#define CutEventNumber_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include <fstream>

class CutEventNumber : public HistoCut
{

public:

  // Parameterized Constructor
  CutEventNumber(EventContainer *obj);
  
  // Destructor
  ~CutEventNumber();

  // Book Histograms
  void BookHistogram();

  // Make cuts and fill histograms
  Bool_t Apply();
  
  // Get the name describing the cut
  inline std::string GetCutName() { return "CutEventNumber"; }


private:

  // Cut parameters
   ifstream * _numberFile;
  int _ENint;
  // Histograms
  myTH1F*   _hEventNumberBefore;  // Lepton charge before cut
  myTH1F*   _hEventNumberAfter;   // Lepton charge after cut

};


#endif







