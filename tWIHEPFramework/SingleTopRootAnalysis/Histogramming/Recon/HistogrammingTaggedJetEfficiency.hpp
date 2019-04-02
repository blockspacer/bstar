/******************************************************************************
 * HistogrammingTaggedJetEfficiency.hpp                                                 *
 *                                                                            *
 * Books and fills histograms                                                 *
 * Used for events passing cuts applied in other classes                      *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of AnalysisMain class                              *
 *    HistogrammingTaggedJetEfficiency()          -- Parameterized Constructor          *
 *    ~HistogrammingTaggedJetEfficiency()         -- Destructor                         *
 *    BookHistogram()                   -- Book histograms                    *
 *    Apply()                           -- Fill histograms only (No Cuts)     *
 *    GetCutName()                      -- Returns "HistogrammingTaggedJetEfficiency"   *
 *                                                                            *
 * Private Data Members of this class                                         *
 *   - lots of histograms of jet properties                                   *
 *                                                                            *
 * History                                                                    *
 *      21 Feb 2007 - Created by P. Ryan                                      *
 *****************************************************************************/

#ifndef HistogrammingTaggedJetEfficiency_h
#define HistogrammingTaggedJetEfficiency_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include <sstream>

class HistogrammingTaggedJetEfficiency : public HistoCut 
{

public:

  // Parameterized Constructor
  HistogrammingTaggedJetEfficiency(EventContainer *obj, Int_t jetindex);
  
  // Destructor
  ~HistogrammingTaggedJetEfficiency();

  inline std::string cutname(Int_t jetindex){
    std::string s;
    std::stringstream out;
    out << jetindex;
    s = out.str();
    std::string name = "HistogrammingTaggedJetEfficiencyW" + s;
    return name;
  }

  // Get the name describing the cut
  inline std::string GetCutName() { return cutname(_jetindex); };

  // Book Histograms
  void BookHistogram();

  // Fill Histograms
  Bool_t Apply();

  inline Int_t JetIndex() const { return _jetindex; };
  inline void SetJetIndex(Int_t jetindex) { _jetindex=jetindex; };
  inline Int_t GetJetIndex() const {return JetIndex();};

private:
  // Histograms declarations 

  myTH1F *_hBLabelW;
  myTH1F *_hLQLabelW;
  myTH1F *_hLQPurLabelW;
  myTH1F *_hCLabelW;
  myTH1F *_hTauLabelW;

  myTH1F *_hAllWtag;

  myTH1F *_hAllW1;
  myTH1F *_hTwoW1;
  myTH1F *_hThreeW1;
  myTH1F *_hFourW1;

  myTH1F *_hAllW2;
  myTH1F *_hTwoW2;
  myTH1F *_hThreeW2;
  myTH1F *_hFourW2;
   
  std::vector<Double_t> _bweightvec;
  std::vector<Double_t> _bweightvectag;
  Int_t _treeflavorlabel;
  Int_t _jetindex;
 
  Double_t _xmin;
  Double_t _xmax;

};


#endif






