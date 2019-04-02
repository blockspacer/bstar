/******************************************************************************
 * CutBkgdQCDShape.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutBkgdQCDShape class                              *
 *    CutBkgdQCDShape()                    -- Parameterized Constructor          *
 *    ~CutBkgdQCDShape()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutBkgdQCDShape"             *
 *                                                                            *
 * Private Data Members of CutBkgdQCDShape                                       *
 *    myTH1F* _hBkgdQCDShapeBefore         -- Hist of BkgdQCDShape before cuts      *
 *    myTH1F* _hBkgdQCDShapeAfter          -- Hist of BkgdQCDShape after cuts       *
 *                                                                            *
 * History                                                                    *
 *      10 Dec 10: Created J. Holzbauer                                       *
 *****************************************************************************/

#ifndef CutBkgdQCDShape_h
#define CutBkgdQCDShape_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutBkgdQCDShape : public HistoCut 
{

public:

  // Parameterized Constructor
  CutBkgdQCDShape(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutBkgdQCDShape();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutBkgdQCDShape"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hBkgdQCDShapeBefore;
  myTH1F* _hBkgdQCDShapeAfter; 

};


#endif







