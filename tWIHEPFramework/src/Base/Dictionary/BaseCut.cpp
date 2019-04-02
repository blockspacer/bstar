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
 *    GetCutName()                      -- Return name of cut as string       *
 *    GetDirectory()                    -- Return directory assoc with task   *
 *    SetDirectory()                    -- Set directory assoc with task      *
 *    PrintCut()                        -- Print Cuts                         *
 *    DoneProcessingLoop()              -- Called after running on full chain *
 *    GetEventContainer()            -- Return Event Container Obj         *
 *    SetEventContainer()            -- Assign Event Container Obj         *
 *                                                                            *
 * Protected Data Members of BaseCut class                                    *
 *    TDirectory *_thisDirecotry         -- Directory in which to put hists   *
 *    EventContainer *_EventContainerObj -- Structure for events              *
 *                                                                            *
 * History                                                                    *
 *      10 Oct 2006 - RS: Ported D0/theory original to Atlas.                 *
 *      14 Nov 2006 - Modified by P. Ryan - cleanup and reorganized           *
 *****************************************************************************/

#include "SingleTopRootAnalysis/Base/Dictionary/BaseCut.hpp"

#include "TROOT.h"

#include <iostream>

// Integrate classes into the Root system
ClassImp(BaseCut);


// Direcotry for Base Cut
static TDirectory *gBaseCutDirectory = NULL;

/******************************************************************************
 * void BaseCut::SetName (const char *name)                                   *
 *                                                                            *
 * Set Name of Object                                                         *
 *                                                                            *
 * Input:  name to set                                                        *
 * Output: None                                                               *
 ******************************************************************************/
void BaseCut::SetName(const char *name)
{
  // Set Name
  TNamed::SetName(name);

  ///
  /// Make sure the BaseCuts directory is there and ready to go!
  ///

  // Make BaseCuts Directory if it doesn't already exist
  if (NULL == gBaseCutDirectory) {
    gBaseCutDirectory = gROOT->mkdir ("BaseCuts");
  } //if

  /// Save current object to BaseCuts directory
  gBaseCutDirectory->Append (this);
} //SetName


/******************************************************************************
 * BaseCut *BaseCut::GetBaseCut (const std::string &name)                     *
 *                                                                            *
 * Return pointer to BaseCut object                                           *
 *                                                                            *
 * Input:  Name of object                                                     *
 * Output: Pointer to BaseCut object or NULL if none exists                   *
 ******************************************************************************/
BaseCut *BaseCut::GetBaseCut (const std::string &name)
{
  BaseCut *BaseCutObj = 0;

  if (NULL == gBaseCutDirectory) {
    BaseCutObj = NULL;
  } //if
  else {
    BaseCutObj = static_cast<BaseCut*>(gBaseCutDirectory->FindObject(name.c_str()));
  } //else

  if (NULL == BaseCutObj ) {
    std::cout << "<BaseCut::GetBaseCut> ERROR: Base cut " << name << " does not exist" << std::endl;
  } //if

  return BaseCutObj;
} //GetBaseCut
