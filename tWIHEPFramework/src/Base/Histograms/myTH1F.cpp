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
 *      19 Oct 2006 - Created by R. Schwienhorst                              *
 *      10 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *       5 Oct 2007 - Changed PassCut and FailCut for MCatNLO weight          *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Base/Histograms/myTH1F.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
#include <iostream>

using namespace std;

/******************************************************************************
 * myTH1F::myTH1F()                                                           *
 *                                                                            *
 * Parameterized Constructor                                                  *
 * Stores histogram in Base Directory                                         *
 *                                                                            *
 * Input:  Name, title, binning info, eventContainerObj                       *
 * Output: None                                                               *
 ******************************************************************************/

// Create the same histogram in several different directories
myTH1F::myTH1F(TDirectory *baseDir,EventContainer *obj,const std::string& name, const std::string& title,
	       Int_t nXBins, Double_t xMin, Double_t xMax): _eventContainerObj(obj)
{
  // Go to Base Directory
  baseDir->cd();

  // If Base Directory is not inclusive make an inclusive directory
  //if(! baseDir->GetKey("inclusive")) baseDir->mkdir("inclusive");  
  // Go to Inclusive directory
  //baseDir->cd("inclusive");

  // Create histogram
  _histogram = new TH1F(name.c_str(), title.c_str(), nXBins, xMin, xMax);  
  _histogram -> Sumw2();
  
} // myTH1F

/******************************************************************************
 * myTH1F::Fill()                                                             *
 *                                                                            *
 * Fill Histogram with weight                                                 *
 * If myWeight < 0 use GlobalWeight                                           *
 *                                                                            *
 * Input:  Value to fill, weight                                              *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH1F::Fill(Double_t x, Double_t myWeight)
{

  // Get Global Weight from EventContainer class.
  Double_t weight = _eventContainerObj -> GetEventWeight();

  // If the input weight is greater than 0, use that instead of global weight
  if( myWeight != -999.0) weight = myWeight;

  // Fill the histogram
  _histogram -> Fill(x, weight);
  return 0;
} //Fill


/******************************************************************************
 * myTH1F::FillWithoutWeight()                                                *
 *                                                                            *
 * Fill Histogram with without weight                                         *
 *                                                                            *
 * Input:  Variable to fill                                                   *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH1F::FillWithoutWeight(Double_t x)
{
  // Fill the histogram
  _histogram -> Fill(x, 1.0);
  return 0;
} //FillWithoutWeight

/******************************************************************************
 * myTH1F::Write()                                                            *
 *                                                                            *
 * Write histogram to output file                                             *
 *                                                                            *
 * Input:  Directory                                                          *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH1F::Write(TDirectory *baseDir)
{
  //baseDir->cd("inclusive");

  _histogram -> Write();  
  
  return 0;
} //Write

/******************************************************************************
 * myTH1F::Copy(myTH1F other)                                                 *
 *                                                                            *
 * Copy histogram                                                             *
 *                                                                            *
 * Input:  Histogram to receive copy                                          *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH1F::Copy(myTH1F other)
{
  // Call the copy method for histogram
  _histogram->Copy(*(other._histogram));  
  
  return 0;
} //Copy

/******************************************************************************
 * myTH1F::Divide(myTH1F other)                                               *
 *                                                                            *
 * Divide histogram                                                           *
 *                                                                            *
 * Input:  Histogram to be divided                                            *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH1F::Divide(myTH1F other)
{
  // Call the divide method for histogram
  _histogram->Divide(other._histogram);  
  
  return 0;
} //Divide

/******************************************************************************
 * myTH1F::SetXAxisTitle(myTH1F other)                                        *
 *                                                                            *
 * Set X Axis Title                                                           *
 *                                                                            *
 * Input:  title string                                                       *
 * Output: None                                                               *
 ******************************************************************************/
void myTH1F::SetXAxisTitle(TString label)
{
  // call the method for histogram
  _histogram->GetXaxis()->SetTitle(label);  
  
  return;
} // SetXAxisTitle

/******************************************************************************
 * myTH1F::SetYAxisTitle(myTH1F other)                                        *
 *                                                                            *
 * Set Y Axis Title                                                           *
 *                                                                            *
 * Input:  title string                                                       *
 * Output: None                                                               *
 ******************************************************************************/
void myTH1F::SetYAxisTitle(TString label)
{
  // call the method for histogram
  _histogram -> GetYaxis() -> SetTitle(label);  
  
  return;
} //SetYAxisTitle

/******************************************************************************
 * myTH1F::SetXAxisBinLabel(Int_t iBin,TString label)                         *
 *                                                                            *
 * Set labels for individual X bins                                           *
 *                                                                            *
 * Input:  bin number, label                                                  *
 * Output: None                                                               *
 ******************************************************************************/
void myTH1F::SetXAxisBinLabel(Int_t iBin,TString label)
{
  // Call the method for each histogram
  _histogram -> GetXaxis() -> SetBinLabel(iBin,label);  
  
  return;
} //SetXAxisBinLabel


