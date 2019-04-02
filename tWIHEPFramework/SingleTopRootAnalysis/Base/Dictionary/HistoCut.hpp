/******************************************************************************
 * HistoCut.hpp                                                               *
 *                                                                            *
 * Abstract class which defines wrappers for TH1F, TH2F, and TGraph.          *
 * Wrappers are of type myTH1F, myTH2F, and TGraph.  The booking and filling  *
 * of the histograms is done in the derived classes.                          *
 *                                                                            *
 * Derived from BaseCut                                                       *
 * Parent to classes such as CutJetN, CutElectron, etc                       *
 *                                                                            *
 * Public Member Functions of HistoCut class                                  *
 *    HistoCut()                        -- Constructor                        *
 *    ~HistoCut()                       -- Destructor                         *
 *    DeclareTH1F()                     -- Wrapper for TH1F                   *
 *    DeclareTH2F()                     -- Wrapper for TH2F                   *
 *    DeclareTGraph()                   -- Wrapper for TGraph                 *
 *                                                                            *
 * Private Data Members of HistoCut class                                     *
 *    TObjArray _histograms              -- Collection of histograms in array *
 *    TObjArray _graphs                  -- Collection of graphs in array     *
 *                                                                            *
 * History                                                                    *
 *      17 June 2015 - Ported from STAR to CMS by Huaqiao ZHANG               *
 *****************************************************************************/

#ifndef HistoCut_h
#define HistoCut_h


#include "SingleTopRootAnalysis/Base/Dictionary/BaseCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <TObjArray.h>
#include "SingleTopRootAnalysis/Base/Histograms/myTH1F.hpp"
#include "SingleTopRootAnalysis/Base/Histograms/myTH2F.hpp"
//#include <TGraph.h>
#include <TDirectory.h>

class HistoCut : public BaseCut
{

public:
 
  // Constructor
  // Pass obj or NULL depending upon existance of object
  HistoCut(EventContainer* obj = NULL);
  
  // Destructor
  virtual ~HistoCut();

  // Book Histograms - Empty for this class
  virtual void BookHistogram();

  // Apply cut and fill histograms - Empty for this class
  virtual Bool_t Apply();
  
  // Save Histograms. Empty for this class
  virtual void SaveHistogram() {};

  // Wrapper for 1-D Histograms
  inline myTH1F* DeclareTH1F(const std::string& name, const std::string& title,
			     Int_t nXBins, Double_t xMin, Double_t xMax) {
    myTH1F* h = new myTH1F(GetDirectory(),GetEventContainer(),name, title, nXBins, xMin, xMax);
    _histograms.Add(h);
    return h; 
  } //DeclareTH1F

  // Wrapper for 2-D Histograms
  inline myTH2F* DeclareTH2F(const std::string& name, const std::string& title,
  			     Int_t nXBins, Double_t xMin, Double_t xMax,
  			     Int_t nYBins, Double_t yMin, Double_t yMax){
    myTH2F *h=new myTH2F(GetDirectory(),GetEventContainer(),name, title, nXBins, xMin, xMax, nYBins, yMin, yMax);
    _histograms.Add(h);
    return h;
  } //DeclareTH2F

  // Wrapper for Graphs
  //inline TGraph* DeclareTGraph(const Int_t& nBins);


private:

  // Collection of histograms in an array
  TObjArray  _histograms;
  // Collection of graphs in an array
  TObjArray  _graphs;
 
  ////////////////////////////////////////////////////////////////////////////////
  // Integrate classes into the Root system
  // Must come at end of class definition
  ClassDef (HistoCut,1);
  ////////////////////// Do Not Put Anything Below this Line /////////////////////

};

#endif
