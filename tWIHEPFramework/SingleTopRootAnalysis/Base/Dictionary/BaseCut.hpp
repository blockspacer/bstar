/******************************************************************************
 * BaseCut.hpp                                                                *
 *                                                                            *
 * Abstract class                                                             *
 *                                                                            *
 * Derived from BaseCut                                                       *
 * Parent to HistCut class which is in turn a parent to classes such as       *
 * CutJetN, CutElectron, etc                                                 *
 *                                                                            *
 * Public Member Functions of BaseCut class                                   *
 *    BaseCut()                         -- Constructor                        *
 *    ~BaseCut()                        -- Destructor                         *
 *    SetName()                         -- Set Name of object                 *
 *    GetBaseCut()                      -- Return a BaseCut Object            *
 *    BookHistogram()                   -- Book Histograms                    *
 *    Apply()                           -- Make cuts and fill histograms      *
 *    SaveHistogram()                   -- Post-processing after event loop   *
 *    GetCutName()                      -- Return name of cut as string       *
 *    GetDirectory()                    -- Return directory assoc with task   *
 *    SetDirectory()                    -- Set directory assoc with task      *
 *    PrintCut()                        -- Print Cuts                         *
 *    DoneProcessingLoop()              -- Called after running on full chain *
 *    GetEventContainer()            -- Return Event Container Obj            *
 *    SetEventContainer()            -- Assign Event Container Obj            *
 *                                                                            *
 * Protected Data Members of BaseCut class                                    *
 *    TDirectory *_thisDirecotry         -- Directory in which to put hists   *
 *    EventContainer *_EventContainerObj -- Structure for events              *
 *                                                                            *
 * History                                                                    *
 *      17 June 2015 -  Port from STAR to CMS by Huaqiao ZHANG                *
 *****************************************************************************/

#ifndef BaseCut_h
#define BaseCut_h

#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>

#include <TDirectory.h>

class CutFlowTable;

class BaseCut : public TNamed
{

public:

  // Constructor
  // Pass object or NULL depending upon existance of object
  inline BaseCut(EventContainer *obj = NULL): _EventContainerObj(obj),_cutFlow(NULL) {};

  // Destructor
  inline virtual ~BaseCut(){};

  // Set Name (used if you want to RESET name) of object
  virtual void SetName(const char *name);

  // Return a pointer to a named base cut
  static BaseCut *GetBaseCut (const std::string &name);

  // Makes Cuts and Fills Histograms -- Purely Virtual
  // Matches the list of cuts with available predifined cuts and applies them 
  // if it is implemented in a class derived from AnalysisCut while fills 
  // histograms if it is implemented in a class derived from HistoCut class. 
  virtual bool Apply() = 0;

  // Book Histograms in derived classes -- Purely Virtual
  virtual void BookHistogram() = 0;

  // Save Histograms in derived classes -- Purely virtual
  virtual void SaveHistogram() = 0;

  // Return the cut list file name -- Purely Virtual
  virtual std::string GetCutName() = 0;

  // Return Directory associated with histogram task
  inline virtual TDirectory *GetDirectory() {return _thisDirectory;}

  // Set Directory associated with histogram task
  inline virtual void SetDirectory( TDirectory* directory) {_thisDirectory = directory;}

  // Print Cuts
  inline virtual void PrintCut(const std::string &indent = ""){};

  // Called when a cut processor has completed running on a full chain.
  inline virtual void DoneProcessingLoop() {};

  // Return pointer to Event Container Object
  inline EventContainer *GetEventContainer() {return _EventContainerObj; };

  // Set Event Container Object
  inline void SetEventContainer(EventContainer *obj) { _EventContainerObj = obj; };

  // Set and get the cut flow table object
  inline void SetCutFlowTable(CutFlowTable *cft) { _cutFlow=cft; };
  inline CutFlowTable* GetCutFlowTable(void) { return _cutFlow; };

protected:

  TDirectory     *_thisDirectory;      // Directory in which to put histograms and graphs
  EventContainer *_EventContainerObj;  // Structure for storing event info
  CutFlowTable   *_cutFlow;            // Provides access to the global cutFlow table with detailed cut info 

  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef (BaseCut,1);
  ////////////////////// Do Not Put Anything Below this Line /////////////////////

};

#endif





