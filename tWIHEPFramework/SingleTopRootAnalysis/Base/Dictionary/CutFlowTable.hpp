/******************************************************************************
 * CutFlowTable.hpp                                                           *
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
 *      17 June 2015 - Huaqiao ZHANG port from STAR                           *
 *****************************************************************************/

#ifndef __CutFlowTable__
#define __CutFlowTable__

#include "TNamed.h"
#include "TClonesArray.h"

#include <vector>
// pair is defined in utility
#include <utility>

#include "SingleTopRootAnalysis/Base/Dictionary/CutFlowLine.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/CutFlowJetLine.hpp"
class EventContainer;


class CutFlowTable : public TNamed {
public:
  // ctor -- both for readback from root and for internal creation.
  CutFlowTable (void);
  CutFlowTable (const TString& name, const TString& title="");
  
  // Reset the list of cuts
  inline void Clear() { _cutLines.clear();_cutJetLines.clear(); };
  
  // Add a cut to the list of cuts
  void AddCutToFlow(const TString& name,const TString& title);

  // Nice print out for this object only
  void PrintCutFlow(Option_t *option = "") const; // *MENU*
  
  // tell this table that this cut passed or failed, using global event weight
  void PassCut(const TString& name);
  void FailCut(const TString& name);
  // tell this table that this cut passed or failed, with custom event weight
  void PassCut(const TString& name, const Double_t& weight);
  void FailCut(const TString& name, const Double_t& weight);
  
  // find the index at which the cut with this name is
  Int_t FindCutName(const TString& name) const;
  
  // Set Event Container Object
  inline void SetEventContainer(EventContainer *obj) { _EventObj = obj; };

  inline std::vector<std::pair<TString,CutFlowLine> > GetCutLines(){return _cutLines;};
  inline std::vector<std::pair<TString,CutFlowJetLine> > GetCutJetLines(){return _cutJetLines;};
private:
  ClassDef (CutFlowTable, 0);
  
  EventContainer *_EventObj; // event object in order to get access to the event weight

  std::vector<std::pair<TString,CutFlowLine> > _cutLines; //|| CutFlowLine. I would love to use a map here, but that screws up the order
  std::vector<std::pair<TString,CutFlowJetLine> > _cutJetLines; //|| CutFlowLine. I would also love to use a map here, but that screws up the order
};

#endif
