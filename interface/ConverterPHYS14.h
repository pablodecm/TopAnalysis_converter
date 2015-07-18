
#ifndef ConverterPHYS14_h
#define ConverterPHYS14_h

#include <vector>
#include <string>
#include <utility>
#include <iostream>
#include <memory>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLorentzVector.h>

#include "Math/GenVector/VectorUtil.h"

// mut_dataformats includes
#include "mut_framework/mut_dataformats/interface/EventInfo.h"
#include "mut_framework/mut_dataformats/interface/Lepton.h"
#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_dataformats/interface/MET.h"

#include "EventHighLevel.h"
#include "mt2bisect.h"


using namespace ROOT::Math;

class ConverterPHYS14 : public TSelector {
  public :

    // asociated with a TTree 
    TTreeReader reader;

    // TTreeReader  for event information
    TTreeReaderValue<int> event; 
    TTreeReaderValue<int> lumi;
    TTreeReaderValue<int> run; 
    TTreeReaderValue<int> nPV;
    
    // TTreeReader pointer to Vertex vector
    TTreeReaderValue<std::vector<float>> T_Vertex_z;
    
    // TTreeReader pointer to MET variables
    TTreeReaderValue<float> pf_met_et;
    TTreeReaderValue<float> pf_met_phi;

    // TTreeReader pointer to Muon vectors
    TTreeReaderValue<std::vector<float>> T_Muon_Px;
    TTreeReaderValue<std::vector<float>> T_Muon_Py;
    TTreeReaderValue<std::vector<float>> T_Muon_Pz;
    TTreeReaderValue<std::vector<float>> T_Muon_Energy;
    TTreeReaderValue<std::vector<int>> T_Muon_Charge;
    TTreeReaderValue<std::vector<bool>> T_Muon_IsGlobalMuon;
    TTreeReaderValue<std::vector<bool>> T_Muon_IsTrackerMuonArbitrated;
    TTreeReaderValue<std::vector<bool>> T_Muon_IsGMPTMuons;
    TTreeReaderValue<std::vector<float>> T_Muon_NormChi2GTrk;
    TTreeReaderValue<std::vector<int>> T_Muon_NValidPixelHitsInTrk;
    TTreeReaderValue<std::vector<int>> T_Muon_NValidHitsGTrk;
    TTreeReaderValue<std::vector<float>> T_Muon_IPwrtAveBSInTrack;
    TTreeReaderValue<std::vector<float>> T_Muon_chargedHadronIsoR04;
    TTreeReaderValue<std::vector<float>> T_Muon_neutralHadronIsoR04;
    TTreeReaderValue<std::vector<float>> T_Muon_photonIsoR04;
    TTreeReaderValue<std::vector<float>> T_Muon_sumPUPtR04;
    TTreeReaderValue<std::vector<float>> T_Muon_dxyInTrack;
    TTreeReaderValue<std::vector<float>> T_Muon_dzInTrack;
    TTreeReaderValue<std::vector<int>> T_Muon_NumOfMatchedStations;
    TTreeReaderValue<std::vector<bool>> T_Muon_IsPFMuon;
    TTreeReaderValue<std::vector<int>> T_Muon_NLayers;


    // TTreeReader pointer to Electron vectors
    TTreeReaderValue<std::vector<float>>   T_Elec_Px;
    TTreeReaderValue<std::vector<float>>   T_Elec_Py;
    TTreeReaderValue<std::vector<float>>   T_Elec_Pz;
    TTreeReaderValue<std::vector<float>>   T_Elec_Energy;
    TTreeReaderValue<std::vector<int>>     T_Elec_Charge;
    TTreeReaderValue<std::vector<float>>   T_Elec_IPwrtPV;
    TTreeReaderValue<std::vector<float>>   T_Elec_dzwrtPV;
    TTreeReaderValue<std::vector<float>>   T_Elec_sumChargedHadronPt;
    TTreeReaderValue<std::vector<float>>   T_Elec_sumNeutralHadronEt;
    TTreeReaderValue<std::vector<float>>   T_Elec_sumPhotonEt;
    TTreeReaderValue<std::vector<float>>   T_Elec_sumPUPt;
    TTreeReaderValue<std::vector<float>>   T_Elec_eSuperClusterOverP;
    TTreeReaderValue<std::vector<float>>   T_Elec_ecalEnergy;
    TTreeReaderValue<std::vector<float>>   T_Elec_dr03TkSumPt;
    TTreeReaderValue<std::vector<float>>   T_Elec_dr03EcalSumEt;
    TTreeReaderValue<std::vector<float>>   T_Elec_dr03HcalSumEt;
    TTreeReaderValue<std::vector<bool>>    T_Elec_isEB;
    TTreeReaderValue<std::vector<bool>>    T_Elec_isEE;
    TTreeReaderValue<std::vector<float>>   T_Elec_chargedHadronIso;
    TTreeReaderValue<std::vector<float>>   T_Elec_neutralHadronIso;
    TTreeReaderValue<std::vector<float>>   T_Elec_photonIso;
    TTreeReaderValue<std::vector<bool>>    T_Elec_passConversionVeto;
    TTreeReaderValue<std::vector<float>>   T_Elec_sigmaIetaIeta;
    TTreeReaderValue<std::vector<float>>   T_Elec_sigmaIetaIetaFull5by5;
    TTreeReaderValue<std::vector<float>>   T_Elec_deltaPhiIn;
    TTreeReaderValue<std::vector<float>>   T_Elec_deltaEtaIn;
    TTreeReaderValue<std::vector<float>>   T_Elec_HtoE;
    TTreeReaderValue<std::vector<float>>   T_Elec_vz;
    TTreeReaderValue<std::vector<int>>     T_Elec_nLost;
    TTreeReaderValue<std::vector<int>>     T_Elec_nHits;
    TTreeReaderValue<std::vector<float>>   T_Elec_SC_Et;
    TTreeReaderValue<std::vector<float>>   T_Elec_SC_Eta;
    TTreeReaderValue<std::vector<float>>   T_Elec_PFElecPt;
    TTreeReaderValue<std::vector<bool>>    T_Elec_isPF;

    // TTreeReader pointer to Jet vectors
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_Px;
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_Py;
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_Pz;
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_Energy;
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_Eta;
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_Tag_CombInclusiveSVtxV2;
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_CharHadEnergyFrac;
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_NeutHadEnergyFrac;
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_CharEmEnergyFrac;
    TTreeReaderValue<std::vector<float>>   T_JetAKCHS_NeutEmEnergyFrac;
    TTreeReaderValue<std::vector<int>>     T_JetAKCHS_ChargedMultiplicity;
    TTreeReaderValue<std::vector<int>>     T_JetAKCHS_nDaughters;
 
 
    // output filename
    std::string o_filename;

    // mut objects to save to TTree
    mut::EventInfo * eventInfo = nullptr;
    std::vector<mut::Jet> * pfjets = nullptr;
    std::vector<mut::Lepton> * leptons = nullptr;
    mut::MET * pfmet = nullptr;
    EventHighLevel * ev_high = nullptr;

    // TTree pointer
    TTree * ttree;
    // Number of events
    long int n_ev_tree_ = 0;
    long int n_ev_proc_ = 0;

    // default constructor
    ConverterPHYS14(TTree * /*tree*/ =0) : 
     event(reader,"T_Event_EventNumber"), 
     lumi(reader, "T_Event_LuminosityBlock"),
     run(reader, "T_Event_RunNumber"),
     nPV(reader, "T_Event_nPU"),
     T_Vertex_z(reader,"T_Vertex_z"),
     // met variables
     pf_met_et(reader,"T_MET_ET"),
     pf_met_phi(reader,"T_MET_Phi"),
     // muon variables
     T_Muon_Px(reader,"T_Muon_Px"),
     T_Muon_Py(reader,"T_Muon_Py"),
     T_Muon_Pz(reader,"T_Muon_Pz"),
     T_Muon_Energy(reader,"T_Muon_Energy"),
     T_Muon_Charge(reader,"T_Muon_Charge"),
     T_Muon_IsGlobalMuon(reader,"T_Muon_IsGlobalMuon"),
     T_Muon_IsTrackerMuonArbitrated(reader,"T_Muon_IsTrackerMuonArbitrated"),
     T_Muon_IsGMPTMuons(reader,"T_Muon_IsGMPTMuons"),
     T_Muon_NormChi2GTrk(reader,"T_Muon_NormChi2GTrk"),
     T_Muon_NValidPixelHitsInTrk(reader,"T_Muon_NValidPixelHitsInTrk"),
     T_Muon_NValidHitsGTrk(reader,"T_Muon_NValidHitsGTrk"),
     T_Muon_IPwrtAveBSInTrack(reader,"T_Muon_IPwrtAveBSInTrack"),
     T_Muon_chargedHadronIsoR04(reader,"T_Muon_chargedHadronIsoR04"),
     T_Muon_neutralHadronIsoR04(reader,"T_Muon_neutralHadronIsoR04"),
     T_Muon_photonIsoR04(reader,"T_Muon_photonIsoR04"),
     T_Muon_sumPUPtR04(reader,"T_Muon_sumPUPtR04"),
     T_Muon_dxyInTrack(reader,"T_Muon_dxyInTrack"),
     T_Muon_dzInTrack(reader,"T_Muon_dzInTrack"),
     T_Muon_NumOfMatchedStations(reader,"T_Muon_NumOfMatchedStations"),
     T_Muon_IsPFMuon(reader,"T_Muon_IsPFMuon"),
     T_Muon_NLayers(reader,"T_Muon_NLayers"),
     // electron variables
     T_Elec_Px(reader,"T_Elec_Px"),
     T_Elec_Py(reader,"T_Elec_Py"),
     T_Elec_Pz(reader,"T_Elec_Pz"),
     T_Elec_Energy(reader,"T_Elec_Energy"),
     T_Elec_Charge(reader,"T_Elec_Charge"),
     T_Elec_IPwrtPV(reader,"T_Elec_IPwrtPV"),
     T_Elec_dzwrtPV(reader,"T_Elec_dzwrtPV"),
     T_Elec_sumChargedHadronPt(reader,"T_Elec_sumChargedHadronPt"),
     T_Elec_sumNeutralHadronEt(reader,"T_Elec_sumNeutralHadronEt"),
     T_Elec_sumPhotonEt(reader,"T_Elec_sumPhotonEt"),
     T_Elec_sumPUPt(reader,"T_Elec_sumPUPt"),
     T_Elec_eSuperClusterOverP(reader,"T_Elec_eSuperClusterOverP"),
     T_Elec_ecalEnergy(reader,"T_Elec_ecalEnergy"),
     T_Elec_dr03TkSumPt(reader,"T_Elec_dr03TkSumPt"),
     T_Elec_dr03EcalSumEt(reader,"T_Elec_dr03EcalSumEt"),
     T_Elec_dr03HcalSumEt(reader,"T_Elec_dr03HcalSumEt"),
     T_Elec_isEB(reader,"T_Elec_isEB"),
     T_Elec_isEE(reader,"T_Elec_isEE"),
     T_Elec_chargedHadronIso(reader,"T_Elec_chargedHadronIso"),
     T_Elec_neutralHadronIso(reader,"T_Elec_neutralHadronIso"),
     T_Elec_photonIso(reader,"T_Elec_photonIso"),
     T_Elec_passConversionVeto(reader,"T_Elec_passConversionVeto"),
     T_Elec_sigmaIetaIeta(reader,"T_Elec_sigmaIetaIeta"),
     T_Elec_sigmaIetaIetaFull5by5(reader,"T_Elec_sigmaIetaIetaFull5by5"),
     T_Elec_deltaPhiIn(reader,"T_Elec_deltaPhiIn"),
     T_Elec_deltaEtaIn(reader,"T_Elec_deltaEtaIn"),
     T_Elec_HtoE(reader,"T_Elec_HtoE"),
     T_Elec_vz(reader,"T_Elec_vz"),
     T_Elec_nLost(reader,"T_Elec_nLost"),
     T_Elec_nHits(reader,"T_Elec_nHits"),
     T_Elec_SC_Et(reader,"T_Elec_SC_Et"),
     T_Elec_SC_Eta(reader,"T_Elec_SC_Eta"),
     T_Elec_PFElecPt(reader,"T_Elec_PFElecPt"),
     T_Elec_isPF(reader,"T_Elec_isPF"),
     // jet variables
     T_JetAKCHS_Px(reader,"T_JetAKCHS_Px"),
     T_JetAKCHS_Py(reader,"T_JetAKCHS_Py"),
     T_JetAKCHS_Pz(reader,"T_JetAKCHS_Pz"),
     T_JetAKCHS_Energy(reader,"T_JetAKCHS_Energy"),
     T_JetAKCHS_Eta(reader,"T_JetAKCHS_Eta"),
     T_JetAKCHS_Tag_CombInclusiveSVtxV2(reader,"T_JetAKCHS_Tag_CombInclusiveSVtxV2"),
     T_JetAKCHS_CharHadEnergyFrac(reader,"T_JetAKCHS_CharHadEnergyFrac"),
     T_JetAKCHS_NeutHadEnergyFrac(reader,"T_JetAKCHS_NeutHadEnergyFrac"),
     T_JetAKCHS_CharEmEnergyFrac(reader,"T_JetAKCHS_CharEmEnergyFrac"),
     T_JetAKCHS_NeutEmEnergyFrac(reader,"T_JetAKCHS_NeutEmEnergyFrac"),
     T_JetAKCHS_ChargedMultiplicity(reader,"T_JetAKCHS_ChargedMultiplicity"),
     T_JetAKCHS_nDaughters(reader,"T_JetAKCHS_nDaughters")
{}


   // destructor
   virtual ~ConverterPHYS14()  { }

   // TSelector functions
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   bool isTightMuon(unsigned iMuon, float minPt = 20.);
   bool isVetoMuon(unsigned iMuon, float minPt = 20.);
   float muonIsolation(unsigned iElec);
   bool isTightElec(unsigned iElec, float minPt = 20.);
   bool isVetoElec(unsigned iElec, float minPt = 20.);
   float elecIsolation(unsigned iElec);
   bool passJetID(unsigned iJet);
   float getMT2( const mut::Lepton & lep0, const mut::Lepton & lep1, const mut::MET & met);


};

#endif

#ifdef ConverterPHYS14_cxx

// each new tree is opened
void ConverterPHYS14::Init(TTree *tree)
{

  reader.SetTree(tree);
  n_ev_tree_ = tree->GetEntries(); 
  n_ev_proc_ = 0;
  std::string tfile_name(tree->GetDirectory()->GetFile()->GetName());
  std::cout << "Processing TTree: " << tfile_name << std::endl;
  std::cout << "  total # events: " << n_ev_tree_ << std::endl;

}

// each new file is opened
Bool_t ConverterPHYS14::Notify()
{

   return kTRUE;
}

#endif // #ifdef ConverterPHYS14_cxx
