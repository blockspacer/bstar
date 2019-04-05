/******************************************************************************
 * myTH2F.hpp                                                                 *
 *                                                                            *
 * Wrapper for TH2F                                                           *
 * Stores histograms in a directory structure                                 *
 *                                                                            *
 * Public Member Functions and data of myTH1F class                           *
 *    myTH2F()                          -- Parameterized Constructor          *
 *    ~myTH2F()                         -- Destructor                         *
 *    Fill()                            -- Fill histograms with weight        *
 *    FillWithoutWeight()               -- Fill histograms without weight     *
 *    Write()                           -- Write to file                      *
 *    Copy()                            -- Copy contents to another hist      *
 *    Divide()                          -- Divide by another histogram        *
 *    SetXAxisTitle()                   -- Set X Axis title                   *
 *    SetYAxisTitle()                   -- Set Y Axis title                   *
 *    SetXAxisBinLabel()                -- Set X Axis individual bin labels   *
 *    SetYAxisBinLabel()                -- Set Y Axis individual bin labels   *
 *    GetNbinsX()                       -- Get number of bins in X            *
 *    TH2F* _histogram                  -- Histogram whihc is operated on     *
 *                                                                            *
 * Private Data Members of myTH1F                                             *
 *    EventContainer *_eventContainerObj -- Event container                   *
 *                                                                            *
 * History                                                                    *
 *      17 June 2015 - Huaqiao ZHANG port from STAR                           *
 *****************************************************************************/
#ifndef myTH2F_h
#define myTH2F_h

#include <TH2.h>
#include <TGraph.h>
#include <string>
#include <TDirectory.h>
#include <TObject.h>

class EventContainer;

class myTH2F: public TObject
{
 public:
  //  constructor
  myTH2F(TDirectory *baseDir,EventContainer *obj,const std::string& name, const std::string& title,
	 Int_t nXBins, Double_t xMin, Double_t xMax,Int_t nYBins, Double_t yMin, Double_t yMax);
  // destructor
  ~myTH2F() {};
  // fill method, fill only some of the histograms depending on npart
  // If the weight is <=0, use the global event weight (default).
  Int_t Fill(Double_t x,Double_t y, Double_t myWeight=-999.0);
  // Fill method that fills a histogram without using any weight
  Int_t FillWithoutWeight(Double_t x,Double_t y);
  //
  // write output to the file
  Int_t Write(TDirectory *baseDir);
  //
  // copy the content of this histogram into another.
  Int_t Copy(myTH2F other);
  //
  // Divide this histogram by the other.
  Int_t Divide(myTH2F other);
  //
  // Set the label for the X axis.
  void SetXAxisTitle(const TString label);
  //
  // Set the label for the Y axis.
  void SetYAxisTitle(const TString label);
  //
  // Set the label for individual bins on the Y axis.
  void SetXAxisBinLabel(Int_t iBin,TString label);
  //
  // Set the label for individual bins on the Y axis.
  void SetYAxisBinLabel(Int_t iBin,TString label);
  //
  // Get the number of bins.
  Int_t GetNbinsX() {return _histogram->GetNbinsX(); };

  // declare the histograms and make them public:
  TH2F* _histogram;           // make every histogram three times:

 private:
  // a pointer to the class storing objects and other stuff
  EventContainer *_eventContainerObj;

};
#endif
