void changed(){
TFile *f1= TFile::Open("PileUpMCdistribution.root","UPDATE");
TObject *changed;
changed = f1->Get("PU_MC");
changed->Write("pileup");



}
