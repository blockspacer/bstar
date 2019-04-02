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
 *      19 Oct 2006 - Created by R. Schwienhorst                              *
 *      10 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Base/Histograms/myTH2F.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

using namespace std;

/******************************************************************************
 * myTH2F::myTH2F()                                                           *
 *                                                                            *
 * Parameterized Constructor                                                  *
 * Stores histogram in Base Directory                                         *
 *                                                                            *
 * Input:  Name, title, binning info, eventContainerObj                       *
 * Output: None                                                               *
 ******************************************************************************/

myTH2F::myTH2F(TDirectory *baseDir,EventContainer *obj,const std::string& name, const std::string& title,
	       Int_t nXBins, Double_t xMin, Double_t xMax,
	       Int_t nYBins, Double_t yMin, Double_t yMax): _eventContainerObj(obj)
{
  // Go to base directory
  baseDir->cd();

  // If Base Directory is not inclusive make an inclusive directory
  // if(! baseDir->GetKey("inclusive")) baseDir->mkdir("inclusive");
  // Go to Inclusive directory
  //baseDir->cd("inclusive");

  // Create histogram
  _histogram=new TH2F(name.c_str(),title.c_str(),nXBins,xMin,xMax,nYBins,yMin,yMax);  
  _histogram -> Sumw2();

} //myTH2F()

/******************************************************************************
 * myTH2F::Fill()                                                             *
 *                                                                            *
 * Fill Histogram with weight                                                 *
 * If myWeight < 0 use GlobalWeight                                           *
 *                                                                            *
 * Input:  Value to fill (x and y), weight                                    *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH2F::Fill(Double_t x, Double_t y, Double_t myWeight)
{
  // Get Global WEight from Event Container class
  Double_t weight=_eventContainerObj->GetEventWeight();

  // If the input weight is greater than 0, use that instead of global weight
  if(myWeight != -999.0) weight = myWeight;

  // Fill the histogram
  _histogram->Fill(x,y,weight);
    
  return 0;
} //Fill()

/******************************************************************************
 * myTH2F::FillWithoutWeight()                                                *
 *                                                                            *
 * Fill Histogram with without weight                                         *
 *                                                                            *
 * Input:  Variable to fill (x and y)                                         *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH2F::FillWithoutWeight(Double_t x, Double_t y)
{
  // Fill the histogram
  _histogram->Fill(x, y, 1.0);

  return 0;
} //myTH2F()

/******************************************************************************
 * myTH2F::Write()                                                            *
 *                                                                            *
 * Write histogram to output file                                             *
 *                                                                            *
 * Input:  Directory                                                          *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH2F::Write(TDirectory *baseDir)
{
  //baseDir->cd("inclusive");

  _histogram->Write();  
  
  return 0;
} //Write()

/******************************************************************************
 * myTH2F::Copy(myTH2F other)                                                 *
 *                                                                            *
 * Copy histogram                                                             *
 *                                                                            *
 * Input:  Histogram to receive copy                                          *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH2F::Copy(myTH2F other)
{
  // Call the copy method for histogram
  _histogram -> Copy(*(other._histogram));  
  
  return 0;
} //Copy()

/******************************************************************************
 * myTH2F::Divide(myTH2F other)                                               *
 *                                                                            *
 * Divide histogram                                                           *
 *                                                                            *
 * Input:  Histogram to be divided                                            *
 * Output: Return 0                                                           *
 ******************************************************************************/
Int_t myTH2F::Divide(myTH2F other)
{
  // Call the divide method for histogram
  _histogram->Divide(other._histogram);  
  
  return 0;
} //Divide()

/******************************************************************************
 * myTH2F::SetXAxisTitle(myTH2F other)                                        *
 *                                                                            *
 * Set X Axis Title                                                           *
 *                                                                            *
 * Input:  title string                                                       *
 * Output: None                                                               *
 ******************************************************************************/
void myTH2F::SetXAxisTitle(TString label)
{
  // Call the method for histogram
  _histogram->GetXaxis()->SetTitle(label);  
  
  return;
} //myTH2F()

/******************************************************************************
 * myTH2F::SetYAxisTitle(myTH1F other)                                        *
 *                                                                            *
 * Set Y Axis Title                                                           *
 *                                                                            *
 * Input:  title string                                                       *
 * Output: None                                                               *
 ******************************************************************************/
void myTH2F::SetYAxisTitle(TString label)
{
  // Call the method for histogram
  _histogram->GetYaxis()->SetTitle(label);  
  
  return;
} //myTH2F

/******************************************************************************
 * myTH2F::SetXAxisBinLabel(Int_t iBin,TString label)                         *
 *                                                                            *
 * Set labels for individual X bins                                           *
 *                                                                            *
 * Input:  bin number, label                                                  *
 * Output: None                                                               *
 ******************************************************************************/
void myTH2F::SetXAxisBinLabel(Int_t iBin,TString label)
{
  // Call the method for histogram
  _histogram->GetXaxis()->SetBinLabel(iBin,label);  
  
  return;
}//SetXAxisBinLabel

/******************************************************************************
 * myTH2F::SetYAxisBinLabel(Int_t iBin,TString label)                         *
 *                                                                            *
 * Set labels for individual X bins                                           *
 *                                                                            *
 * Input:  bin number, label                                                  *
 * Output: None                                                               *
 ******************************************************************************/
void myTH2F::SetYAxisBinLabel(Int_t iBin,TString label)
{
  // Call the method for each histogram
  _histogram->GetYaxis()->SetBinLabel(iBin,label);  
  
  return;
}//SetYAxisBinLabel
