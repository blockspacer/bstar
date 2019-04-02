/******************************************************************************
 * CutListProcessor.hpp                                                       *
 *                                                                            *
 * Process the list of cuts and histogramming classes                         *
 *                                                                            *
 * Derived from BaseCut                                                       *
 * Parent to AnalysisMain                                                     *
 *                                                                            *
 * Public Member Functions of CutNJet class                                   *
 *    CutListProcessor()                -- Parameterized Constructor          *
 *    ~CutListProcessor()               -- Destructor                         *
 *    AddCut()                          -- Add a cut to the list to process   *
 *    BookHistogram()                   -- Book Histograms                    *
 *    Apply()                           -- Make cuts and fill histograms      *
 *    GetCutName()                      -- Get name of cut                    *
 *    PrintCut()                        -- Print name of cut                  *
 *    SetDebug()                        -- Set Debug on or off                *
 *    DeclareTH2F()                     -- Wrapper for TH2F                   *
 *    AddDefaultBaseCut                 -- Add template to every CutListProc  *
 *    template<class T>                 -- Template to be added to every list *
 *    GetCurrentPackage()               -- Returns current package            *
 *    GetCurrentEventContainer()        -- Returns event read                 *
 *                                                                            *
 * Private Data Members of CutNJet class                                      *
 *    string name                         -- Name of CutListProcessor         *
 *    vector<BaseCut *> _vectorCutAndHist -- Vector of cuts and histograms    *
 *    TH1F *_hCutAndHist                  -- Number of times each cut runs    *
 *    Bool_t _alwaysPass                  -- Always returns kFALSE            *
 *    CleanUpRun()                        -- Clean up after loop              *
 *    Bool_t _gDegub                      -- If debug is on                   *
 *                                                                            *
 *    template <class T> class AutoInsert_t                                   *
 *    template <class T> static AutoInsert_t<T> *AutoInsert_t_maker           *
 *    static std::vector<AutoInsert*> &defaultCuts(void)                      *
 *    static BaseCut *_currentBaseCut;                                        *
 *    static EventContainer *_currentEventContainer;                          *
 *    CutFlowTable _cutFlowDetails        -- store cut flow statistics.       *
 *                                                                            *
 * History                                                                    *
 *      17 June 2006 - Huaqiao ZHANG : Port from STAR                         *
 *****************************************************************************/

#ifndef __CutListProcessor__
#define __CutListProcessor__

#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/CutFlowTable.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/BaseCut.hpp"
#include "TH1F.h"

#include <vector>

class CutListProcessor : public BaseCut {
public:
  
  // Parameterized Constructor
  CutListProcessor(const std::string &name, const TString& condition="AND",EventContainer *obj=NULL);

  // Destructor
  inline virtual ~CutListProcessor() {};
  
  // Add a new cut to the list to processes
  void AddCut (BaseCut *cut);
  
  // Book Histograms - before looping over events
  virtual void BookHistogram(void);

  Bool_t Apply() { return Apply(NULL); };
  // Make cuts and fill histograms - once per event
  virtual Bool_t Apply (EventContainer *obj);
  
  // Save Histograms - after looping over events
  virtual void SaveHistogram(void);

  // Get name describing the cut
  virtual std::string GetCutName (void);
  
  // Print info from cut classes
  virtual void PrintCut(const std::string &indent = "");
  
  // Print cut details for each one
  void PrintCutDetails();
  
 
  // Set debug
  static void SetDebug(Bool_t debugFlag) {_gDebug = debugFlag;}
  
  // Template to be added as default BaseCut
  template<class T>

  // Add a BaseCut to the start of every CutListProcessor
  // Added class must have at least one method T::operator() (CutListProcessor *p)
  // Called after construction of CutListProcessor
  static inline void AddDefaultBaseCut(T &creator){
    defaultCuts().push_back(AutoInsert_t_maker (creator));
  } //AddDefaultBaseCut
  
  // Return package currently running
  static inline BaseCut *GetGurrentPackage() {
    return _currentBaseCut;
  } // GetCurrentPackage

  // Get current event container
  static inline EventContainer *GetCurrentEventContainer () {
    return _currentEventContainer;
  } // EventContainer
  
  // Return the number of cuts in the cut lsit
  inline int GetNumberOfCuts() {
    return _vectorCutAndHist.size();
  }

private:
  int NumCutLines;
  std::vector<std::pair<TString,CutFlowLine> >  _cutLines;

  // Name of CutListProcessor
  const std::string _name;

  // Vector of tasks (makeing cuts and filling histograms) to be performed
  std::vector<BaseCut *> _vectorCutAndHist;

  // a cut flow table that cuts can use to fill in more details for each cut
  CutFlowTable _cutFlowDetails;

  // Histogram of times Cut and Histogram classes are called
  TH1F *_hCutAndHist;
  TH1F *_hCutEff;
  TH1F *_hCutEffNoEffect;
  
  
  Bool_t _alwaysPass;   // Always return true?
  Bool_t _doOr;         // Do an OR instead of the default AND of the conditions?
  
  // Loop over Vector and perform tasks on each element
  void CleanUpRun(std::vector<BaseCut *>::iterator start,
		  std::vector<BaseCut *>::iterator end);
  
  // Debug Flag
  static Bool_t _gDebug;
  
  //////////////////////////////////////////////////////////////
  // Template classes, etc needed to support Cut ListProcessor
  //////////////////////////////////////////////////////////////
  class AutoInsert {
  public:
    // Destructor
    virtual inline ~AutoInsert(){};
    // Purely Virtual function
    virtual void AddCuts (CutListProcessor *p) = 0;
  }; //AutoInsert
  
  // Template class derived from AutoInsert
  template <class T> class AutoInsert_t : public AutoInsert {
  public:
    inline AutoInsert_t (T &Create_functor) : _creator (Create_functor){};
    inline void AddCuts (CutListProcessor *p) {_creator (p);}
  private:
    T _creator;
  }; // AutoInsert_t
  
  // Return AutoInsert_t
  template<class T> static AutoInsert_t<T> *AutoInsert_t_maker (T &creator)
  {
    return new AutoInsert_t<T> (creator);
  } // AutoInsert_t_maker

  // Add default cuts
  static std::vector<AutoInsert*> &defaultCuts();
  
  // Current Base Cut and Event Container - used for tracking purposes
  static BaseCut *_currentBaseCut;
  static EventContainer *_currentEventContainer;
};


#endif
