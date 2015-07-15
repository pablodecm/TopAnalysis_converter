
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

// mut_dataformats includes
#include "mut_framework/mut_dataformats/interface/EventInfo.h"
#include "mut_framework/mut_dataformats/interface/Lepton.h"
#include "mut_framework/mut_dataformats/interface/Jet.h"
#include "mut_framework/mut_dataformats/interface/MET.h"



class ConverterPHYS14 : public TSelector {
  public :

    // asociated with a TTree 
    TTreeReader reader;

    // TTreeReader  for event information
    TTreeReaderValue<int> event; 
    TTreeReaderValue<int> lumi;
    TTreeReaderValue<int> run; 
    TTreeReaderValue<int> nPV;
    
    // TreeReader pointer to MET variables
    TTreeReaderValue<float> pf_met_et;
    TTreeReaderValue<float> pf_met_phi;

    // output filename
    std::string o_filename;

    // mut objects to save to TTree
    mut::EventInfo * eventInfo = nullptr;
    std::vector<mut::Jet> * pfjets = nullptr;
    std::vector<mut::Lepton> * muons = nullptr;
    mut::MET * pfmet = nullptr;

    // TTree pointer
    TTree * ttree;

    // default constructor
    ConverterPHYS14(TTree * /*tree*/ =0) : 
     event(reader,"T_Event_EventNumber"), 
     lumi(reader, "T_Event_LuminosityBlock"),
     run(reader, "T_Event_RunNumber"),
     nPV(reader, "T_Event_nPU"),
     pf_met_et(reader,"T_MET_ET"),
     pf_met_phi(reader,"T_MET_Phi")
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


};

#endif

#ifdef ConverterPHYS14_cxx

// each new tree is opened
void ConverterPHYS14::Init(TTree *tree)
{
  reader.SetTree(tree);
}

// each new file is opened
Bool_t ConverterPHYS14::Notify()
{

   return kTRUE;
}

#endif // #ifdef ConverterPHYS14_cxx
