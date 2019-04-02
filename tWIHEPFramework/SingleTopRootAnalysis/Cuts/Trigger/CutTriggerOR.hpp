/******************************************************************************
 * CutTriggerOR.hpp                                                           *
 *                                                                            *
 * Cuts on 1 or more CutTriggerORs                                            *
 * Selection perforemed with full trigger (L1, L2, EF)                        *
 *                                                                            *
 * Derived from HistoCut which is in turn derived from BaseCut                *
 *                                                                            *
 * Public Member Functions of CutMuonPt1 class                                *
 *    CutTriggerOROR()                    -- Parameterized Constructor        *
 *    ~CutTriggerOROR()                   -- Destructor                       *
 *    BookHistogram()                -- Book histograms                       *
 *    Apply()                        -- Apply cuts and fill histograms        *
 *    GetCutName()                   -- Returns "CutTriggerOROR"              *
 *                                                                            *
 * Private Data Members of CutTriggerOROR class                               *
 *    myTH1F* _hL1MuonBefore;        -- Hist before trigger cut               *
 *    myTH1F* _hL1MuonAfter;         -- Hist after trigger cut                *
 *                                                                            *
 *    Double_t _L1MuonName;          -- Name of trigger cut                   *
 *                                                                            *
 * History                                                                    *
 *      14 Nov 2007 - Created by P. Ryan                                      *
 *      17 Jun 2008 - P. Ryan - updated for v13.0.30                          *
 *      27 Jan 2009 - P. Ryan - updated for v14.05.006                        *
 *****************************************************************************/

#ifndef CutTriggerOR_h
#define CutTriggerOR_h

#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"

#include <string>
#include <vector>
#include <sstream>

using namespace std;

class CutTriggerOR : public HistoCut 
{

public:

  // Parameterized Constructor
  CutTriggerOR(EventContainer *EventContainerObj, std::string CutTriggerORNamePassed = "NOTHINGPASSED");
  
  // Destructor
  ~CutTriggerOR();

  // Get the name of describing the cut
  ///////////// Change this
  inline std::string GetCutName() { return "CutTriggerOR_" + _TriggerOrName; };
  /**

  // Return Individual CutTriggerORs
  // ******************** Level 1 ********************
  // ********** Muon **********
  inline Bool_t L1_MU15()               { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU15; };
  inline Bool_t L1_MU0()                { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU0; };
  inline Bool_t L1_MU6()                { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU6; };
  inline Bool_t L1_MU10()               { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU10; };
  inline Bool_t L1_MU20()               { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU20; };
  inline Bool_t L1_MU0_COMM()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU0_COMM; };
  inline Bool_t L1_MU6_J5()             { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU6_J5; };
  inline Bool_t L1_MU0_EM3()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU0_EM3; };
  inline Bool_t L1_MU0_UNPAIRED()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU0_UNPAIRED; };
  inline Bool_t L1_MU6_UNPAIRED()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU6_UNPAIRED; };
  inline Bool_t L1_MU0_EMPTY()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU0_EMPTY; };
  inline Bool_t L1_MU6_EMPTY()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU6_EMPTY; };
  inline Bool_t L1_MU10_EMPTY()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU10_EMPTY; };
  inline Bool_t L1_MU0_COMM_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU0_COMM_EMPTY; };
  //  inline Bool_t L1_MU0_COMM_UNPAIRED()  { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU0_COMM_UNPAIRED; };
  inline Bool_t L1_MU6_EM3()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU6_EM3; };
  inline Bool_t L1_MU0_J5()             { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_MU0_J5; };

  // ********** Electron **********
  inline Bool_t L1_EM2()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_EM2; };
  inline Bool_t L1_EM3()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_EM3; };
  inline Bool_t L1_EM3_UNPAIRED()  { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_EM3_UNPAIRED; };
  inline Bool_t L1_EM4()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_EM4; };
  inline Bool_t L1_EM5()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_EM5; };
  inline Bool_t L1_EM5I()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_EM5I; };
  inline Bool_t L1_EM10()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_EM10; };
  inline Bool_t L1_EM10I()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_EM10I; };
  inline Bool_t L1_EM14()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_EM14; };

  // ********** Jet **********
  inline Bool_t L1_J5()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J5; };
  inline Bool_t L1_J5_UNPAIRED()   { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J5_UNPAIRED; };
  inline Bool_t L1_J5_EMPTY()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J5_EMPTY; };
  inline Bool_t L1_J10()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J10; };
  inline Bool_t L1_J15()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J15; };
  inline Bool_t L1_J10_win6()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J10_win6; };
  inline Bool_t L1_J15_MV()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J15_MV; };
  inline Bool_t L1_J15_LV()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J15_LV; };
  inline Bool_t L1_J15_NL()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J15_NL; };
  inline Bool_t L1_J30()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J30; };
  inline Bool_t L1_J75()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J75; };
  inline Bool_t L1_J55()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J55; };
  inline Bool_t L1_J10_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J10_EMPTY; };
  inline Bool_t L1_J30_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_J30_EMPTY; };
  //  inline Bool_t L1_FJ3_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_FJ3_EMPTY; };
  inline Bool_t L1_2J5()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_2J5; };
  inline Bool_t L1_2J10()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_2J10; };
  inline Bool_t L1_2J10_win6()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_2J10_win6; };
  inline Bool_t L1_2J15()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_2J15; };
  inline Bool_t L1_JE120()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_JE120; };
  //  inline Bool_t L1_JE280()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_JE280; };
  //  inline Bool_t L1_FJ3_UNPAIRED()  { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_FJ3_UNPAIRED; };
  //  inline Bool_t L1_FJ18()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_FJ18; };
  //  inline Bool_t L1_FJ3()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_FJ3; };
  //inline Bool_t L1_2FJ18()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_2FJ18; };
  //inline Bool_t L1_2FJ3()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_2FJ3; };
  inline Bool_t L1_3J10()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_3J10; };

  // ********** Tau **********
  inline Bool_t L1_TAU5()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU5; };
  inline Bool_t L1_TAU6()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU6; };
  inline Bool_t L1_TAU6I()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU6I; };
  inline Bool_t L1_TAU8()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU8; };
  inline Bool_t L1_TAU11()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU11; };
  inline Bool_t L1_TAU11I()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU11I; };
  inline Bool_t L1_TAU20()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU20; };
  inline Bool_t L1_TAU30()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU30; };
  inline Bool_t L1_2TAU6I()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_2TAU6I; };
  inline Bool_t L1_2TAU5()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_2TAU5; };
  inline Bool_t L1_2TAU6()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_2TAU6; };
  inline Bool_t L1_TAU5_EMPTY()  { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU5_EMPTY; };
  inline Bool_t L1_TAU6_XE10()   { return GetEventContainer() -> GetTriggerTree() -> Passed_L1_TAU6_XE10; };



  // ********** Missing ET **********
  // Must implement 
  **/
  /**
    // ******************** Level 2 ********************
  // ********** Muon **********
  inline Bool_t L2_mu4()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4; };
  inline Bool_t L2_mu4_MSonly()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_MSonly; };
  inline Bool_t L2_mu4_SiTrk()                     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_SiTrk; };
  inline Bool_t L2_mu4_tile()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_tile; };
  inline Bool_t L2_mu4_trod()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_trod; };
  inline Bool_t L2_mu4_NoIDTrkCut()                { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_NoIDTrkCut; };
  inline Bool_t L2_mu4_MG()                        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_MG; };
  inline Bool_t L2_mu4_tile_SiTrk()                { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_tile_SiTrk; };
  inline Bool_t L2_mu4_trod_SiTrk()                { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_trod_SiTrk; };
  inline Bool_t L2_mu4_passHLT()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_passHLT; };
  inline Bool_t L2_mu6()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6; };
  inline Bool_t L2_mu6_MSonly()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_MSonly; };
  inline Bool_t L2_mu6_SiTrk()                     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_SiTrk; };
  inline Bool_t L2_mu6_MG()                        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_MG; };
  inline Bool_t L2_mu6_passHLT   ()                { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_passHLT   ; };
  inline Bool_t L2_mu10()                          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu10; };
  inline Bool_t L2_mu10i_loose()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu10i_loose; };
  inline Bool_t L2_mu10_MG()                       { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu10_MG; };
  inline Bool_t L2_mu10_MSonly()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu10_MSonly; };
  inline Bool_t L2_mu10_SiTrk()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu10_SiTrk; };
  inline Bool_t L2_mu10_NoIDTrkCut()               { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu10_NoIDTrkCut; };
  inline Bool_t L2_mu10_passHLT()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu10_passHLT; };
  inline Bool_t L2_mu13()                          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu13; };
  inline Bool_t L2_mu15()                          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu15; };
  inline Bool_t L2_mu20()                          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu20; };
  inline Bool_t L2_mu20_passHLT()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu20_passHLT; };
  inline Bool_t L2_mu20_MSonly()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu20_MSonly; };
  inline Bool_t L2_mu4_mu6()                       { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_mu6; };
  inline Bool_t L2_mu4_j20_matched()               { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_j20_matched; };
  inline Bool_t L2_mu6_Jpsie5e3()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_Jpsie5e3; };
  inline Bool_t L2_mu4mu6_DiMu()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4mu6_DiMu; };
  inline Bool_t L2_mu4_j10_matched()               { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_j10_matched; };
  inline Bool_t L2_mu4_j5_matched()                { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_j5_matched; };
  inline Bool_t L2_mu4_DiMu_FS()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_DiMu_FS; };
  inline Bool_t L2_mu6_DiMu_FS()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_DiMu_FS; };
  inline Bool_t L2_mu4_DiMu()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_DiMu; };
  inline Bool_t L2_mu6_DiMu()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_DiMu; };
  inline Bool_t L2_mu4_MSonly_Trk_Jpsi_loose()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_MSonly_Trk_Jpsi_loose; };
  inline Bool_t L2_mu4_MSonly_Trk_Upsi_loose_FS()  { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_MSonly_Trk_Upsi_loose_FS; };
  inline Bool_t L2_mu6_Trk_Jpsi_loose()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_Trk_Jpsi_loose; };
  inline Bool_t L2_mu6_Trk_Upsi_loose_FS()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_Trk_Upsi_loose_FS; };
  inline Bool_t L2_mu4_Trk_Jpsi_loose()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_Trk_Jpsi_loose; };
  inline Bool_t L2_mu4_Trk_Upsi_loose_FS()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_Trk_Upsi_loose_FS; };
  inline Bool_t L2_mu6_Trk_Jpsi_loose_FS()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_Trk_Jpsi_loose_FS; };
  inline Bool_t L2_mu4_DiMu_FS_noOS()              { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_DiMu_FS_noOS; };
  inline Bool_t L2_mu6_DiMu_FS_noOS()              { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_DiMu_FS_noOS; };
  inline Bool_t L2_mu4_L2MSonly_EFFS_passL2()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_L2MSonly_EFFS_passL2; };
  inline Bool_t L2_mu6_L2MSonly_EFFS_passL2()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_L2MSonly_EFFS_passL2; };
  inline Bool_t L2_mu10_L2MSonly_EFFS_passL2()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu10_L2MSonly_EFFS_passL2; };
  inline Bool_t L2_mu4_MSonly_EFFS_passL2()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_MSonly_EFFS_passL2; };
  inline Bool_t L2_mu6_MSonly_EFFS_passL2()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu6_MSonly_EFFS_passL2; };
  inline Bool_t L2_mu10_MSonly_EFFS_passL2()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu10_MSonly_EFFS_passL2; };
  //  inline Bool_t L2_mu4_MB2_noL2_EFFS()             { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_mu4_MB2_noL2_EFFS; };



  // ********** Electron **********
  inline Bool_t L2_e3_NoCut()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e3_NoCut; };
  inline Bool_t L2_e3_NoCut_IdScan()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e3_NoCut_IdScan; };
  inline Bool_t L2_e3_NoCut_TRT()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e3_NoCut_TRT; };
  inline Bool_t L2_e5_medium()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_medium; };
  inline Bool_t L2_e5_NoCut_Ringer()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_NoCut_Ringer; };
  inline Bool_t L2_e5_loose_NoIDTrkCut() { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_loose_NoIDTrkCut; };
  inline Bool_t L2_e5_medium_IdScan()    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_medium_IdScan; };
  inline Bool_t L2_e5_NoCut()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_NoCut; };
  inline Bool_t L2_e5_NoCut_IdScan()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_NoCut_IdScan; };
  inline Bool_t L2_e5_NoCut_TRT()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_NoCut_TRT; };
  inline Bool_t L2_e5_NoCut_SiTrk()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_NoCut_SiTrk; };
  inline Bool_t L2_e5_NoCut_FwdBackTrk() { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_NoCut_FwdBackTrk; };
  inline Bool_t L2_e5_NoCut_cells()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_NoCut_cells; };
  inline Bool_t L2_e5_NoCut_L2SW()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_NoCut_L2SW; };
  inline Bool_t L2_e5_EFfullcalo()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_EFfullcalo; };
  inline Bool_t L2_e5_medium_TRT()       { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_medium_TRT; };
  inline Bool_t L2_e5_loose()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_loose; };
  inline Bool_t L2_e5_loose_IdScan()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_loose_IdScan; };
  inline Bool_t L2_e5_loose_TRT()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e5_loose_TRT; };
  inline Bool_t L2_e10_loose()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e10_loose; };
  inline Bool_t L2_e10i_loose()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e10i_loose; };
  inline Bool_t L2_e10_loose_passL2()    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e10_loose_passL2; };
  inline Bool_t L2_e10_loose_passEF()    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e10_loose_passEF; };
  inline Bool_t L2_e10_medium()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e10_medium; };
  inline Bool_t L2_e15i_loose()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e15i_loose; };
  inline Bool_t L2_e15_loose()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_e15_loose; };

 
  // ********** Jet  **********
  inline Bool_t L2_j7()                              { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j7; };
  inline Bool_t L2_j7_cosmic()                       { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j7_cosmic; };
  inline Bool_t L2_j7_noise()                        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j7_noise; };
  inline Bool_t L2_j15()                             { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j15; };
  inline Bool_t L2_j15_noise()                       { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j15_noise; };
  inline Bool_t L2_j20_Trackless_HV()                { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j20_Trackless_HV; };
  inline Bool_t L2_j30()                             { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j30; };
  inline Bool_t L2_j30_L1TAU_HV()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j30_L1TAU_HV; };
  inline Bool_t L2_j60()                             { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j60; };
  inline Bool_t L2_j90()                             { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j90; };
  inline Bool_t L2_j130()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_j130; };
  inline Bool_t L2_2j7()                             { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_2j7; };
  inline Bool_t L2_2j7_deta3_5()                     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_2j7_deta3_5; };
  inline Bool_t L2_2j7_deta5()                       { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_2j7_deta5; };
  inline Bool_t L2_2j15()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_2j15; };
  inline Bool_t L2_2j15_deta3_5()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_2j15_deta3_5; };
  inline Bool_t L2_2j30()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_2j30; };
  inline Bool_t L2_2j30_deta3_5()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_2j30_deta3_5; };
  inline Bool_t L2_3j15()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_3j15; };
  //  inline Bool_t L2_FJ3()                             { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_FJ3; };
  //  inline Bool_t L2_FJ18()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_FJ18; };
  inline Bool_t L2_2MUL1_j40_HV()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_2MUL1_j40_HV; };
  inline Bool_t L2_SeededStreamerL1CaloJet()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_SeededStreamerL1CaloJet; };
  inline Bool_t L2_L1ItemStreamer_L1_J15_LV()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J15_LV; };
  inline Bool_t L2_L1ItemStreamer_L1_3J10()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_3J10; };
  //  inline Bool_t L2_L1ItemStreamer_L1_JE280()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_JE280; };
  inline Bool_t L2_L1ItemStreamer_L1_JE120()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_JE120; };
  //  inline Bool_t L2_L1ItemStreamer_L1_FJ3()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_FJ3; };
  inline Bool_t L2_L1ItemStreamer_L1_2J10()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_2J10; };
  //  inline Bool_t L2_L1ItemStreamer_L1_FJ3_UNPAIRED()  { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_FJ3_UNPAIRED; };
  inline Bool_t L2_L1ItemStreamer_L1_2J15()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_2J15; };
  inline Bool_t L2_L1ItemStreamer_L1_2J5()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_2J5; };
  inline Bool_t L2_L1ItemStreamer_L1_J30()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J30; };
  //  inline Bool_t L2_L1ItemStreamer_L1_FJ18()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_FJ18; };
  inline Bool_t L2_L1ItemStreamer_L1_MU6_J5()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_MU6_J5; };
  inline Bool_t L2_L1ItemStreamer_L1_J5_UNPAIRED()   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J5_UNPAIRED; };
  inline Bool_t L2_L1ItemStreamer_L1_J15_NL()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J15_NL; };
  inline Bool_t L2_L1ItemStreamer_L1_J10()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J10; };
  inline Bool_t L2_L1ItemStreamer_L1_J15()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J15; };
  //  inline Bool_t L2_L1ItemStreamer_L1_FJ3_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_FJ3_EMPTY; };
  inline Bool_t L2_L1ItemStreamer_L1_J5()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J5; };
  inline Bool_t L2_L1ItemStreamer_L1_J5_EMPTY()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J5_EMPTY; };
  inline Bool_t L2_L1ItemStreamer_L1_2J10_win6()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_2J10_win6; };
  inline Bool_t L2_L1ItemStreamer_L1_J30_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J30_EMPTY; };
  inline Bool_t L2_L1ItemStreamer_L1_J75()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J75; };
  inline Bool_t L2_L1ItemStreamer_L1_J10_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J10_EMPTY; };
  //  inline Bool_t L2_L1ItemStreamer_L1_2FJ18()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_2FJ18; };
  inline Bool_t L2_L1ItemStreamer_L1_J15_MV()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J15_MV; };
  inline Bool_t L2_L1ItemStreamer_L1_J55()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J55; };
  inline Bool_t L2_L1ItemStreamer_L1_J10_win6()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_J10_win6; };
  //  inline Bool_t L2_L1ItemStreamer_L1_2FJ3()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_2FJ3; };
  inline Bool_t L2_L1ItemStreamer_L1_MU0_J5()        { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_MU0_J5; };

  // ********** Tau **********
  inline Bool_t L2_tau12_loose()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tau12_loose; };
  inline Bool_t L2_tau16_loose()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tau16_loose; };
  inline Bool_t L2_tau16i_NoIDTrkCut()                { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tau16i_NoIDTrkCut; };
  inline Bool_t L2_tauNoCut()                         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tauNoCut; };
  inline Bool_t L2_tauNoCut_NoIDTrkCut()              { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tauNoCut_NoIDTrkCut; };
  inline Bool_t L2_tauNoCut_SiTrk()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tauNoCut_SiTrk; };
  inline Bool_t L2_tauNoCut_cosmic()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tauNoCut_cosmic; };
  inline Bool_t L2_tauNoCut_cells()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tauNoCut_cells; };
  inline Bool_t L2_tauNoCut_hasTrk()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tauNoCut_hasTrk; };
  inline Bool_t L2_tauNoCut_hasTrk_SiTrk()            { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_tauNoCut_hasTrk_SiTrk; };
  inline Bool_t L2_trk9_Central_Tau_IDCalib()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_trk9_Central_Tau_IDCalib; };
  inline Bool_t L2_trk9_Fwd_Tau_IDCalib()             { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_trk9_Fwd_Tau_IDCalib; };
  inline Bool_t L2_trk9_Central_Tau_SiTrack_IDCalib() { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_trk9_Central_Tau_SiTrack_IDCalib; };
  inline Bool_t L2_trk9_Fwd_Tau_SiTrack_IDCalib()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_trk9_Fwd_Tau_SiTrack_IDCalib; };
  inline Bool_t L2_SeededStreamerL1CaloTau()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_SeededStreamerL1CaloTau; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU6I()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU6I; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU11()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU11; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU6_XE10()      { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU6_XE10; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU5_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU5_EMPTY; };
  inline Bool_t L2_L1ItemStreamer_L1_2TAU5()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_2TAU5; };
  inline Bool_t L2_L1ItemStreamer_L1_2TAU6()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_2TAU6; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU20()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU20; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU5()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU5; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU6()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU6; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU8()           { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU8; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU11I()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU11I; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU30()          { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU30; };
  inline Bool_t L2_L1ItemStreamer_L1_2TAU6I()         { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_2TAU6I; };
  inline Bool_t L2_L1ItemStreamer_L1_TAU5_UNPAIRED()  { return GetEventContainer() -> GetTriggerTree() -> Passed_L2_L1ItemStreamer_L1_TAU5_UNPAIRED; };




  **/
  /**
  // ******************** Event Filter ********************
  // ********** Muon **********
  inline Bool_t EF_mu4()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4; };
  inline Bool_t EF_mu4_MSonly()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_MSonly; };
  inline Bool_t EF_mu4_SiTrk()                     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_SiTrk; };
  inline Bool_t EF_mu4_tile()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_tile; };
  inline Bool_t EF_mu4_trod()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_trod; };
  inline Bool_t EF_mu4_NoIDTrkCut()                { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_NoIDTrkCut; };
  inline Bool_t EF_mu4_MG()                        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_MG; };
  inline Bool_t EF_mu4_tile_SiTrk()                { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_tile_SiTrk; };
  inline Bool_t EF_mu4_trod_SiTrk()                { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_trod_SiTrk; };
  inline Bool_t EF_mu4_passHLT()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_passHLT; };
  inline Bool_t EF_mu6()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6; };
  inline Bool_t EF_mu6_MSonly()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_MSonly; };
  inline Bool_t EF_mu6_SiTrk()                     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_SiTrk; };
  inline Bool_t EF_mu6_MG()                        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_MG; };
  inline Bool_t EF_mu6_passHLT   ()                { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_passHLT   ; };
  inline Bool_t EF_mu10()                          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu10; };
  inline Bool_t EF_mu10i_loose()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu10i_loose; };
  inline Bool_t EF_mu10_MG()                       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu10_MG; };
  inline Bool_t EF_mu10_MSonly()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu10_MSonly; };
  inline Bool_t EF_mu10_SiTrk()                    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu10_SiTrk; };
  inline Bool_t EF_mu10_NoIDTrkCut()               { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu10_NoIDTrkCut; };
  inline Bool_t EF_mu10_passHLT()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu10_passHLT; };
  inline Bool_t EF_mu13()                          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu13; };
  inline Bool_t EF_mu15()                          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu15; };
  inline Bool_t EF_mu20()                          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu20; };
  inline Bool_t EF_mu20_passHLT()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu20_passHLT; };
  inline Bool_t EF_mu20_MSonly()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu20_MSonly; };
  inline Bool_t EF_mu4_mu6()                       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_mu6; };
  inline Bool_t EF_mu4_j20_matched()               { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_j20_matched; };
  inline Bool_t EF_mu6_Jpsie5e3()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_Jpsie5e3; };
  inline Bool_t EF_mu4mu6_DiMu()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4mu6_DiMu; };
  inline Bool_t EF_mu4_j10_matched()               { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_j10_matched; };
  inline Bool_t EF_mu4_j5_matched()                { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_j5_matched; };
  inline Bool_t EF_mu4_DiMu_FS()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_DiMu_FS; };
  inline Bool_t EF_mu6_DiMu_FS()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_DiMu_FS; };
  inline Bool_t EF_mu4_DiMu()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_DiMu; };
  inline Bool_t EF_mu6_DiMu()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_DiMu; };
  inline Bool_t EF_mu4_MSonly_Trk_Jpsi_loose()     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_MSonly_Trk_Jpsi_loose; };
  inline Bool_t EF_mu4_MSonly_Trk_Upsi_loose_FS()  { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_MSonly_Trk_Upsi_loose_FS; };
  inline Bool_t EF_mu6_Trk_Jpsi_loose()            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_Trk_Jpsi_loose; };
  inline Bool_t EF_mu6_Trk_Upsi_loose_FS()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_Trk_Upsi_loose_FS; };
  inline Bool_t EF_mu4_Trk_Jpsi_loose()            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_Trk_Jpsi_loose; };
  inline Bool_t EF_mu4_Trk_Upsi_loose_FS()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_Trk_Upsi_loose_FS; };
  inline Bool_t EF_mu6_Trk_Jpsi_loose_FS()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_Trk_Jpsi_loose_FS; };
  inline Bool_t EF_mu4_DiMu_FS_noOS()              { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_DiMu_FS_noOS; };
  inline Bool_t EF_mu6_DiMu_FS_noOS()              { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_DiMu_FS_noOS; };
  inline Bool_t EF_mu4_L2MSonly_EFFS_passL2()      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_L2MSonly_EFFS_passL2; };
  inline Bool_t EF_mu6_L2MSonly_EFFS_passL2()      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_L2MSonly_EFFS_passL2; };
  inline Bool_t EF_mu10_L2MSonly_EFFS_passL2()     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu10_L2MSonly_EFFS_passL2; };
  inline Bool_t EF_mu4_MSonly_EFFS_passL2()        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_MSonly_EFFS_passL2; };
  inline Bool_t EF_mu6_MSonly_EFFS_passL2()        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu6_MSonly_EFFS_passL2; };
  inline Bool_t EF_mu10_MSonly_EFFS_passL2()       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu10_MSonly_EFFS_passL2; };
  //  inline Bool_t EF_mu4_MB2_noL2_EFFS()             { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_mu4_MB2_noL2_EFFS; };


  // ********** Electron **********
  inline Bool_t EF_e3_NoCut()            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e3_NoCut; };
  inline Bool_t EF_e3_NoCut_IdScan()     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e3_NoCut_IdScan; };
  inline Bool_t EF_e3_NoCut_TRT()        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e3_NoCut_TRT; };
  inline Bool_t EF_e5_medium()           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_medium; };
  inline Bool_t EF_e5_NoCut_Ringer()     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_NoCut_Ringer; };
  inline Bool_t EF_e5_loose_NoIDTrkCut() { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_loose_NoIDTrkCut; };
  inline Bool_t EF_e5_medium_IdScan()    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_medium_IdScan; };
  inline Bool_t EF_e5_NoCut()            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_NoCut; };
  inline Bool_t EF_e5_NoCut_IdScan()     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_NoCut_IdScan; };
  inline Bool_t EF_e5_NoCut_TRT()        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_NoCut_TRT; };
  inline Bool_t EF_e5_NoCut_SiTrk()      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_NoCut_SiTrk; };
  inline Bool_t EF_e5_NoCut_FwdBackTrk() { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_NoCut_FwdBackTrk; };
  inline Bool_t EF_e5_NoCut_cells()      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_NoCut_cells; };
  inline Bool_t EF_e5_NoCut_L2SW()       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_NoCut_L2SW; };
  inline Bool_t EF_e5_EFfullcalo()       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_EFfullcalo; };
  inline Bool_t EF_e5_medium_TRT()       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_medium_TRT; };
  inline Bool_t EF_e5_loose()            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_loose; };
  inline Bool_t EF_e5_loose_IdScan()     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_loose_IdScan; };
  inline Bool_t EF_e5_loose_TRT()        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e5_loose_TRT; };
  inline Bool_t EF_e10_loose()           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e10_loose; };
  inline Bool_t EF_e10i_loose()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e10i_loose; };
  inline Bool_t EF_e10_loose_passL2()    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e10_loose_passL2; };
  inline Bool_t EF_e10_loose_passEF()    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e10_loose_passEF; };
  inline Bool_t EF_e10_medium()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e10_medium; };
  inline Bool_t EF_e15_loose()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e15_loose; };
  inline Bool_t EF_e15i_loose()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_e15i_loose; };

  // ********** Jet  **********
  inline Bool_t EF_j10()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j10; };
  inline Bool_t EF_j10_larcalib()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j10_larcalib; };
  inline Bool_t EF_j10_noise()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j10_noise; };
  inline Bool_t EF_j10_c4()                         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j10_c4; };
  inline Bool_t EF_j10_cosmic()                     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j10_cosmic; };
  inline Bool_t EF_j20()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j20; };
  inline Bool_t EF_j20_noise()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j20_noise; };
  inline Bool_t EF_j20_c4()                         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j20_c4; };
  inline Bool_t EF_j30_Trackless_HV()               { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j30_Trackless_HV; };
  inline Bool_t EF_j30_L1TAU_HV()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j30_L1TAU_HV; };
  inline Bool_t EF_j40()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j40; };
  inline Bool_t EF_j80()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j80; };
  inline Bool_t EF_j80_larcalib()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j80_larcalib; };
  inline Bool_t EF_j140()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j140; };
  inline Bool_t EF_j200()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_j200; };
  inline Bool_t EF_2j10()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_2j10; };
  inline Bool_t EF_2j10_deta3_5()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_2j10_deta3_5; };
  inline Bool_t EF_2j10_deta5()                     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_2j10_deta5; };
  inline Bool_t EF_2j20()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_2j20; };
  inline Bool_t EF_2j20_deta3_5()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_2j20_deta3_5; };
  inline Bool_t EF_2j40()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_2j40; };
  inline Bool_t EF_2j40_deta3_5()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_2j40_deta3_5; };
  inline Bool_t EF_3j20()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_3j20; };
  //  inline Bool_t EF_FJ3_larcalib()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_FJ3_larcalib; };
  //  inline Bool_t EF_FJ3()                            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_FJ3; };
  //inline Bool_t EF_FJ18()                           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_FJ18; };
  //inline Bool_t EF_FJ18_larcalib()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_FJ18_larcalib; };
  inline Bool_t EF_2MUL1_j40_HV()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_2MUL1_j40_HV; };
  inline Bool_t EF_SeededStreamerL1CaloJet()        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_SeededStreamerL1CaloJet; };
  inline Bool_t EF_L1ItemStreamer_L1_J15_LV()       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J15_LV; };
  inline Bool_t EF_L1ItemStreamer_L1_3J10()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_3J10; };
  //inline Bool_t EF_L1ItemStreamer_L1_JE280()        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_JE280; };
  inline Bool_t EF_L1ItemStreamer_L1_JE120()        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_JE120; };
  //inline Bool_t EF_L1ItemStreamer_L1_FJ3()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_FJ3; };
  inline Bool_t EF_L1ItemStreamer_L1_2J10()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_2J10; };
  //inline Bool_t EF_L1ItemStreamer_L1_FJ3_UNPAIRED() { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_FJ3_UNPAIRED; };
  inline Bool_t EF_L1ItemStreamer_L1_2J15()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_2J15; };
  inline Bool_t EF_L1ItemStreamer_L1_2J5()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_2J5; };
  inline Bool_t EF_L1ItemStreamer_L1_J30()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J30; };
  //inline Bool_t EF_L1ItemStreamer_L1_FJ18()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_FJ18; };
  inline Bool_t EF_L1ItemStreamer_L1_MU6_J5()       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_MU6_J5; };
  inline Bool_t EF_L1ItemStreamer_L1_J5_UNPAIRED()  { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J5_UNPAIRED; };
  inline Bool_t EF_L1ItemStreamer_L1_J15_NL()       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J15_NL; };
  inline Bool_t EF_L1ItemStreamer_L1_J10()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J10; };
  inline Bool_t EF_L1ItemStreamer_L1_J15()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J15; };
  //inline Bool_t EF_L1ItemStreamer_L1_FJ3_EMPTY()    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_FJ3_EMPTY; };
  inline Bool_t EF_L1ItemStreamer_L1_J5()           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J5; };
  inline Bool_t EF_L1ItemStreamer_L1_J5_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J5_EMPTY; };
  inline Bool_t EF_L1ItemStreamer_L1_2J10_win6()    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_2J10_win6; };
  inline Bool_t EF_L1ItemStreamer_L1_J30_EMPTY()    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J30_EMPTY; };
  inline Bool_t EF_L1ItemStreamer_L1_J75()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J75; };
  inline Bool_t EF_L1ItemStreamer_L1_J10_EMPTY()    { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J10_EMPTY; };
  //  inline Bool_t EF_L1ItemStreamer_L1_2FJ18()        { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_2FJ18; };
  inline Bool_t EF_L1ItemStreamer_L1_J15_MV()       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J15_MV; };
  inline Bool_t EF_L1ItemStreamer_L1_J55()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J55; };
  inline Bool_t EF_L1ItemStreamer_L1_J10_win6()     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_J10_win6; };
  //  inline Bool_t EF_L1ItemStreamer_L1_2FJ3()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_2FJ3; };
  inline Bool_t EF_L1ItemStreamer_L1_MU0_J5()       { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_MU0_J5; };


  // ********** Tau **********
  inline Bool_t EF_tau12_loose()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tau12_loose; };
  inline Bool_t EF_tau16_loose()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tau16_loose; };
  inline Bool_t EF_tau16i_NoIDTrkCut()                { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tau16i_NoIDTrkCut; };
  inline Bool_t EF_tauNoCut()                         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tauNoCut; };
  inline Bool_t EF_tauNoCut_NoIDTrkCut()              { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tauNoCut_NoIDTrkCut; };
  inline Bool_t EF_tauNoCut_SiTrk()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tauNoCut_SiTrk; };
  inline Bool_t EF_tauNoCut_cosmic()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tauNoCut_cosmic; };
  inline Bool_t EF_tauNoCut_cells()                   { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tauNoCut_cells; };
  inline Bool_t EF_tauNoCut_hasTrk()                  { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tauNoCut_hasTrk; };
  inline Bool_t EF_tauNoCut_hasTrk_SiTrk()            { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_tauNoCut_hasTrk_SiTrk; };
  inline Bool_t EF_SeededStreamerL1CaloTau()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_SeededStreamerL1CaloTau; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU6I()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU6I; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU11()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU11; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU6_XE10()      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU6_XE10; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU5_EMPTY()     { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU5_EMPTY; };
  inline Bool_t EF_L1ItemStreamer_L1_2TAU5()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_2TAU5; };
  inline Bool_t EF_L1ItemStreamer_L1_2TAU6()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_2TAU6; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU20()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU20; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU5()           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU5; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU6()           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU6; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU8()           { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU8; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU11I()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU11I; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU30()          { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU30; };
  inline Bool_t EF_L1ItemStreamer_L1_2TAU6I()         { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_2TAU6I; };
  inline Bool_t EF_L1ItemStreamer_L1_TAU5_UNPAIRED()  { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_L1ItemStreamer_L1_TAU5_UNPAIRED; };

  // ********** Missing ET **********
  inline Bool_t EF_xe20()                      { return GetEventContainer() -> GetTriggerTree() -> Passed_EF_xe20; };


  // Return Combined L1, L2, EF CutTriggerORs
  // Chains from   http://olya.web.cern.ch/olya/menu/view_menu.php?name=MC_lumi1E31_simpleL1Calib&flavor=15.5.X.Y-VAL&nightly=rel_0&tags=&project=AtlasP1HLT#Chain_L1_MU0_EM3
  // ********** Muon **********
  inline Bool_t Item_mu4()             { return ( L1_MU0()    && EF_mu4()         ); };
  inline Bool_t Item_mu6()             { return ( L1_MU6()         && EF_mu6()         ); };
  inline Bool_t Item_mu10()            { return ( L1_MU10()       && EF_mu10()        ); };
  inline Bool_t Item_mu10i_loose()     { return ( L1_MU10() && EF_mu10i_loose() ); };
  inline Bool_t Item_mu13()            { return ( L1_MU10()        && EF_mu13()        ); };
  inline Bool_t Item_mu15()            { return ( L1_MU15()        && EF_mu15()        ); };
  inline Bool_t Item_mu20()            { return ( L1_MU20()        && EF_mu20()        ); };
  inline Bool_t Item_mu10_MSonly()     { return ( EF_mu10_MSonly() ); };


  // ********** Electron **********
  inline Bool_t Item_e5_medium()           { return (  L1_EM3()          && EF_e5_medium()    ); };  
  inline Bool_t Item_e10_loose()           { return (  L1_EM5()            && EF_e10_loose()    ); };  
  inline Bool_t Item_e10_medium()          { return (  L1_EM5()          && EF_e10_medium()   ); };  
  //  inline Bool_t Item_e10i_loose()          { return  (  L1_EM5()     && L2_e10i_loose()        && EF_e10i_loose()   ); };  
  inline Bool_t Item_e15_loose()           { return  (  L1_EM10()            && EF_e15_loose()    ); };  
  //inline Bool_t Item_e15i_loose()          { return  (  L1_EM5()     && L2_e15i_loose()        && EF_e15i_loose()   ); };  
  //  inline Bool_t Item_e15_medium()          { return (  L1_EM10()    && L2_e15_medium()        && EF_e15_medium()       ); };  


  // ********** Jet **********
  inline Bool_t Item_j40()           { return (  EF_j40()   ); };  
  inline Bool_t Item_j80()           { return (  EF_j80()   ); };  
  inline Bool_t Item_j140()          { return (  EF_j140()  ); };  
  inline Bool_t Item_2j40()          { return (  EF_2j40()  ); };  


  // ********** Tau **********
  // Need to implement
  // ********** B Jet **********
  // Need to implement
  // ********** Missing Et **********
  inline Bool_t Item_xe20()          { return (  EF_xe20()  ); };  

  // ********** Total Et **********
  // Need to implement
  // ********** Total Jet Et **********
  // Need to implement

  // ********** Combinations **********
  inline Bool_t Item_j80_xe20()          { return (  EF_j80()  && EF_xe20()  ); };  
  inline Bool_t Item_2j40_xe20()         { return (  EF_2j40() && EF_xe20()  ); };  
  **/
  // ********** No trigger returns True **********
  inline Bool_t NONE()        { return kTRUE; };
 

  // Book Histograms
  void BookHistogram();

  // Apply cut and fill histogram
  Bool_t Apply();


private:

  // Negation Flat
  Bool_t NegationFlag;

  // Histograms
  myTH1F* _hTriggerOrBefore;    // Histogram Muon Before CutTriggerOR Cut
  myTH1F* _hTriggerOrAfter;          // Histogram Muon After CutTriggerOR Cut

  // Cut parameters
  std::vector<TString> _TriggerOrNameVector;    // Vector of trigger names
  std::string _TriggerOrName;                   // A single trigger name

  std::vector<std::string> ValidTriggerNameVector;  // Vector of trigger names
  std::vector<std::string> TriggerCutVector;        // Vector of trigger names used in cut
  
  
};


#endif
