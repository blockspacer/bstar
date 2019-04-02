#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "TTreePlayer.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iostream>
#include <TLorentzVector.h>
#include <cmath>
#include "TMath.h"
#include "TH2.h"
void sumfile(const char * Input = ""){
 TFile *fileBTagEfficiency = TFile::Open(Input,"UPDATE");
 TH2D *bFlavEfficiency;TH2D *cFlavEfficiency;TH2D *lightFlavEfficiency;
 TH2D *DEN1;TH2D *DEN2;TH2D *DEN3;
 bFlavEfficiency= (TH2D*)fileBTagEfficiency->Get("histM_b_passing_");
 DEN1 =  (TH2D*)fileBTagEfficiency->Get("histM_b_total_");
 cFlavEfficiency = (TH2D*)fileBTagEfficiency->Get("histM_c_passing_");
 DEN2 = (TH2D*)fileBTagEfficiency->Get("histM_c_total_");
 lightFlavEfficiency = (TH2D*)fileBTagEfficiency->Get("histM_udsg_passing_");
 DEN3 = (TH2D*)fileBTagEfficiency->Get("histM_udsg_total_");
 bFlavEfficiency->Divide(DEN1);
 cFlavEfficiency->Divide(DEN2);
 lightFlavEfficiency->Divide(DEN3);
 bFlavEfficiency->Write("bFlavEfficiency");
 cFlavEfficiency->Write("cFlavEfficiency");
 lightFlavEfficiency->Write("lightFlavEfficiency");
 fileBTagEfficiency->Close(); 
}
