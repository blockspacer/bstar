/******************************************************************************
 * CutQCDShape.hpp                                                           *
 *                                                                            *
 * Cuts on Missing Et                                                         *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutQCDShape class                              *
 *    CutQCDShape()                    -- Parameterized Constructor          *
 *    ~CutQCDShape()                   -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Apply cuts and fill histograms     *
 *    GetCutName()                      -- Returns "CutQCDShape"             *
 *                                                                            *
 * Private Data Members of CutQCDShape                                       *
 *    myTH1F* _hQCDShapeBefore         -- Hist of QCDShape before cuts      *
 *    myTH1F* _hQCDShapeAfter          -- Hist of QCDShape after cuts       *
 *                                                                            *
 * History                                                                    *
 *      10 Dec 10: Created J. Holzbauer                                       *
 *****************************************************************************/

#ifndef CutQCDShape_h
#define CutQCDShape_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"


class CutQCDShape : public HistoCut 
{

public:

  // Parameterized Constructor
  CutQCDShape(EventContainer *EventContainerObj);
  
  // Destructor
  ~CutQCDShape();

  // Get the name of describing the cut
  inline std::string GetCutName() { return "CutQCDShape"; };

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:


  // Histograms
  myTH1F* _hQCDShapeBefore;
  myTH1F* _hQCDShapeAfter; 

};


#endif







