#define ConverterPHYS14_cxx
#include "../interface/ConverterPHYS14.h"

// start of query (executed on client)
void ConverterPHYS14::Begin(TTree * /*tree*/)
{

   std::string option = GetOption();

   std::size_t i_ofile = option.find("ofile="); 
   if (i_ofile != std::string::npos) {
     std::size_t length = (option.find(";", i_ofile) -  option.find("=", i_ofile) - 1);
     o_filename = option.substr(option.find("=", i_ofile)+1 , length );
   } else {
     o_filename = "output.root";
   }

   std::cout << "Output filename: " << o_filename << std::endl;
}

// right after begin (executed on slave)
void ConverterPHYS14::SlaveBegin(TTree * /*tree*/)
{

   TString option = GetOption();

  ttree = new TTree("tree","Physics Object based TTree");

  ttree->Branch("eventInfo","mut::EventInfo", &eventInfo, 64000,1);
  ttree->Branch("pfmet","mut::MET", &pfmet, 64000,1);
  ttree->Branch("pfjets","std::vector<mut::Jet>", &pfjets, 64000,1);
  ttree->Branch("muons","std::vector<mut::Lepton>", &muons, 64000,1);

  fOutput->Add(ttree);  

}

// for each entry of the TTree
Bool_t ConverterPHYS14::Process(Long64_t entry)
{

  // set TTreeReader entry
  reader.SetEntry(entry);

  // create and fill EventInfo
  eventInfo = new mut::EventInfo(*event, *lumi, *run);
  eventInfo->setNumPV(*nPV);
  std::vector<std::pair<std::string, bool>> filterPairs;
  eventInfo->setFilterPairs(filterPairs);
  std::vector<std::pair<std::string, float>> weightPairs;
  eventInfo->setWeightPairs(weightPairs);

  // create and fill mut MET
  pfmet = new mut::MET(*pf_met_et,0.0, *pf_met_phi, *pf_met_et);

  // create and fill Jet vector
  pfjets = new std::vector<mut::Jet>;
  /*for (int i=0;i < n_pfjet;i++) {*/
    //if (**pfjet_energy[i] > 0) {
      //pfjets->emplace_back(**pfjet_pt[i], **pfjet_eta[i], **pfjet_phi[i], **pfjet_energy[i]);
      //pfjets->back().setPartonFlavour(**jet_flavour[i]);
      //std::vector<std::pair<std::string, float>> disPairs;
      //for  (auto const &a_disc : pfjet_disc.at(i)) {
        //disPairs.emplace_back(a_disc.first,**a_disc.second);
      //}
      //pfjets->back().setDiscriminatorPairs(disPairs);
    //}
  //}

    // create and fill Muon vector
  muons = new std::vector<mut::Lepton>;
    //for (int i=0;i < n_muon;i++) {
    //if (**muon_energy[i] > 0) {
       //muons->emplace_back(**muon_pt[i], **muon_eta[i], **muon_phi[i], **muon_energy[i]);
    //} else {
      //return false;
    //}
    //// add isolation variables
    //std::vector<std::pair<std::string, float>> isoPairs;
    //isoPairs.emplace_back("relIso", **muon_relIso[i]);
    //muons->back().setLeptonIsoPairs(isoPairs);
  /*}*/

  ttree->Fill();

  delete eventInfo;
  delete pfjets;
  delete pfmet;
  delete muons;

  return true;
}

// all entries have been processed (executed in slave)
void ConverterPHYS14::SlaveTerminate()
{

}

// last function called (on client)
void ConverterPHYS14::Terminate()
{
  TTree* ttree = dynamic_cast<TTree *>(fOutput->FindObject(Form("tree"))); 

  TFile *o_file = new TFile(o_filename.c_str(), "RECREATE");
  if ( o_file->IsOpen() ) std::cout << "File is opened successfully" << std::endl;
  ttree->Write();
  fOutput->Clear(); 
}
