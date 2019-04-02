

//CutFraction.hpp  
 
//Based off HistogrammingElectron.hpp

//Brad Schoenrock & Jake Parsons - June 03,2011

//this program cuts the data into thirds for MVA


#include "SingleTopRootAnalysis/Base/Dictionary/HistoCut.hpp"
#include "SingleTopRootAnalysis/Base/Dictionary/EventContainer.hpp"
 
class CutFraction : public HistoCut
{
 
  public:
 
        CutFraction(EventContainer *obj, int WhichThirdPassed);
 
        ~CutFraction();
 
        inline std::string GetCutName()
        {
                return "CutFraction";
        };
 

        void BookHistogram();
 
        Bool_t Apply();

        std::vector<Electron>   temp_electrons;
        std::vector<Electron>   temp_tightElectrons;
        std::vector<Electron>   temp_vetoElectrons;
        std::vector<Electron>   temp_ptetaElectrons;
        std::vector<Electron>   temp_isolatedElectrons;
        std::vector<Electron>   temp_unIsolatedElectrons;
        std::vector<Muon>       temp_muons;
        std::vector<Muon>       temp_tightMuons;
        std::vector<Muon>       temp_vetoMuons;
        std::vector<Muon>       temp_ptetaMuons;
        std::vector<Muon>       temp_isolatedMuons;
        std::vector<Muon>       temp_unIsolatedMuons;
        std::vector<Tau>        temp_taus;
        std::vector<Jet>        temp_jets;
        std::vector<Jet>        temp_jetors;//overlap removal vector for electron/jets
        std::vector<Jet>        temp_jetms;//overlap removal vector for muons
        std::vector<Jet>        temp_taggedJets;
        std::vector<Jet>        temp_unTaggedJets;
        std::vector<Jet>        temp_bLabeledJets;
        std::vector<Jet>        temp_cLabeledJets;
        std::vector<Jet>        temp_tauLabeledJets;
        std::vector<Jet>        temp_lightQuarkLabeledJets;
        std::vector<Neutrino>   temp_neutrinos;
 
  int WhichThird;
 
  private:
 
        myTH1F* _hNObj;
 
};
 
