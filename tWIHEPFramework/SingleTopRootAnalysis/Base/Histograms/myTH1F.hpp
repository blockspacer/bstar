/******************************************************************************
 * myTH1F.hpp                                                                 *
 *                                                                            *
 * Wrapper for TH1F                                                           *
 * Stores histograms in a directory structure                                 *
 *                                                                            *
 * Public Member Functions and data of myTH1F class                           *
 *    myTH1F()                          -- Parameterized Constructor          *
 *    ~myTH1F()                         -- Destructor                         *
 *    Fill()                            -- Fill histograms with weight        *
 *    FillWithoutWeight()               -- Fill histograms without weight     *
 *    Write()                           -- Write to file                      *
 *    Copy()                            -- Copy contents to another hist      *
 *    Divide()                          -- Divide by another histogram        *
 *    SetXAxisTitle()                   -- Set X Axis title                   *
 *    SetYAxisTitle()                   -- Set Y Axis title                   *
 *    SetXAxisBinLabel()                -- Set X Axis individual bin labels   *
 *    GetNbinsX()                       -- Get number of bins in X            *
 *    TH1F* _histogram                  -- Histogram whihc is operated on     *
 *                                                                            *
 * Private Data Members of myTH1F                                             *
 *    EventContainer *_eventContainerObj -- Event container                   *
 *                                                                            *
 * History                                                                    *
 *      17 June 2015 - Huaqiao ZHANG port from STAR                           *
 *****************************************************************************/
#ifndef myTH1F_h
#define myTH1F_h

#include <TH1.h>
#include <TGraph.h>
#include <string>
#include <TDirectory.h>
#include <TObject.h>

class EventContainer;

class myTH1F: public TObject
{
public:
  //  Parameterized constructor
  myTH1F(TDirectory *baseDir,EventContainer *obj,const std::string& name, const std::string& title,
	 Int_t nXBins, Double_t xMin, Double_t xMax);

  // Destructor
  ~myTH1F() {};

  // Fill
  // Fill only some of the histograms depending on npart
  // If the weight is <=0, use the global event weight (default).
  Int_t Fill(Double_t x, Double_t myWeight = -999.0);

  // Fill without weight
  Int_t FillWithoutWeight(Double_t x);

  // Write output to file
  Int_t Write(TDirectory *baseDir);

  // Copy the content of this histogram into another.
  Int_t Copy(myTH1F other);

  // Divide this histogram by the other.
  Int_t Divide(myTH1F other);

  // Set the title for the X axis.
  void SetXAxisTitle(const TString label);

  // Set the title for the Yaxis.
  void SetYAxisTitle(const TString label);

  // Set the label for individual bins on the X axis.
  void SetXAxisBinLabel(Int_t iBin,TString label);

  // Get the number of bins.
  Int_t GetNbinsX() {return _histogram->GetNbinsX(); };

  // Histogram which is operaton upon
  TH1F* _histogram;

 private:
  // Pointer to EventContainer
  EventContainer *_eventContainerObj;

};

#endif
