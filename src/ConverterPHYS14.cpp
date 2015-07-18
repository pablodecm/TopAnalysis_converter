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
  ttree->Branch("leptons","std::vector<mut::Lepton>", &leptons, 64000,1);
  ttree->Branch("ev_high","EventHighLevel", &ev_high, 64000,1);

  fOutput->Add(ttree);  

}

// for each entry of the TTree
Bool_t ConverterPHYS14::Process(Long64_t entry)
{

  // set TTreeReader entry
  reader.SetEntry(entry);
  n_ev_proc_++;

  if (n_ev_proc_%10000 == 0) 
   std::cout << " --- " << n_ev_proc_ << "/" << n_ev_tree_ << " --- " << std::endl;  

  // create and fill EventInfo
  eventInfo = new mut::EventInfo(*event, *lumi, *run);
  eventInfo->setNumPV(*nPV);
  std::vector<std::pair<std::string, bool>> filterPairs;
  eventInfo->setFilterPairs(filterPairs);
  std::vector<std::pair<std::string, float>> weightPairs;
  eventInfo->setWeightPairs(weightPairs);

  // create and fill mut MET
  pfmet = new mut::MET(*pf_met_et,0.0, *pf_met_phi, *pf_met_et);

  // create and fill Lepton vector
  leptons = new std::vector<mut::Lepton>;

    // Muon Selection
    int nMuon = T_Muon_Energy->size();
    for (int i = 0; i < nMuon; i++ ) {
        if(isTightMuon(i)) {
            leptons->emplace_back();
            leptons->back().SetPxPyPzE(T_Muon_Px->at(i), T_Muon_Py->at(i),
                                       T_Muon_Pz->at(i), T_Muon_Energy->at(i));
            leptons->back().setCharge(T_Muon_Charge->at(i));
            leptons->back().setPdgId(T_Muon_Charge->at(i)*13);
        } // end if tight muon
    } // end muon loop

    int nElec = T_Elec_Energy->size();
    for (int i=0; i < nElec; i++) {
        if(isTightElec(i)) {
            leptons->emplace_back();
            leptons->back().SetPxPyPzE(T_Elec_Px->at(i), T_Elec_Py->at(i),
                                       T_Elec_Pz->at(i), T_Elec_Energy->at(i));
            leptons->back().setCharge(T_Elec_Charge->at(i));
            leptons->back().setPdgId(T_Elec_Charge->at(i)*11);
        } // end if tight elec
    } // end electron loop
   
    // get number of loose leptons
    int nVetoLepton = 0;
    for (int i = 0; i < nMuon; i++ ) {
        if(isVetoMuon(i)) {
            nVetoLepton++;
        } // end if veto muon
    } // end muon loop
    for (int i = 0; i < nElec; i++ ) {
        if(isVetoElec(i)) {
            nVetoLepton++;
        } // end if veto elec
    } // end muon loop

    // create and fill Jet vector
    pfjets = new std::vector<mut::Jet>;
 
    // Jet Selection
    int nJet= T_JetAKCHS_Energy->size();
    for ( int i = 0 ; i < nJet; i++) {
        mut::Jet jet;
        jet.SetPxPyPzE( T_JetAKCHS_Px->at(i), T_JetAKCHS_Py->at(i),
                        T_JetAKCHS_Pz->at(i), T_JetAKCHS_Energy->at(i));
        // check if jets are cleaned of selected leptons
        bool isClean = true;
        for (auto lept : *leptons ) {
          if (VectorUtil::DeltaR(jet, lept) < 0.4) isClean = false;
        } // end cleaning check
        if (jet.Et() > 30. &&
            fabs(jet.eta()) < 2.4 &&
            isClean && passJetID(i) ) {
              pfjets->emplace_back(jet);
              std::vector<std::pair<std::string, float>> disPairs;
              disPairs.emplace_back("CombInclusiveSVtxV2", T_JetAKCHS_Tag_CombInclusiveSVtxV2->at(i));
              pfjets->back().setDiscriminatorPairs(disPairs);
     } // end good jet loop
  } // end jet loop

  // sort in pt order
  std::sort(leptons->begin(), leptons->end(), []( const mut::Lepton  & lhs,  const mut::Lepton & rhs) 
                                                  { return lhs.pt() < rhs.pt();  });

  ev_high = new EventHighLevel();

  // event selection 
  if (leptons->size() == 2) { // two tight leptons
    if (leptons->at(0).charge()*leptons->at(1).charge()  < 0) { // opposite sign
      double mll = VectorUtil::InvariantMass(leptons->at(0), leptons->at(1));
      if (mll > 20.0 && (nVetoLepton < 3) && (pfjets->size() > 1)) { // low mll and extralepton veto
        if ((std::abs(leptons->at(0).pdgId()) + std::abs(leptons->at(1).pdgId())) == 24 ) {
            
            ev_high->met = pfmet->Et();
            ev_high->dilept_inv_mass = mll;
            ev_high->dilepton_MT2 = getMT2( leptons->at(0), leptons->at(1), *pfmet);
            ev_high->d_phi_l_l = std::abs(VectorUtil::DeltaPhi(leptons->at(0), leptons->at(1)));
            ev_high->d_phi_min_met_j = 10;
            for ( const auto & jet : *pfjets) {
              float d_phi_min_met_j =  std::abs(VectorUtil::DeltaPhi(*pfmet, jet));
              if (d_phi_min_met_j < ev_high->d_phi_min_met_j) 
                ev_high->d_phi_min_met_j = d_phi_min_met_j;
            }
            ev_high->d_phi_met_l0 = std::abs(VectorUtil::DeltaPhi(*pfmet, leptons->at(0)));
            ev_high->d_phi_met_l1 = std::abs(VectorUtil::DeltaPhi(*pfmet, leptons->at(1)));
            ev_high->d_phi_met_metll = std::abs(VectorUtil::DeltaPhi(*pfmet, 
                                                    *pfmet+leptons->at(0)+leptons->at(1)));
            
            // only emu and mue channel
            ttree->Fill();
        }
      }
    }
  }

  delete eventInfo;
  delete pfjets;
  delete pfmet;
  delete leptons;
  delete ev_high;

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



bool ConverterPHYS14::isTightMuon( unsigned iMuon, float minPt) 
{
  TLorentzVector lep;
  lep.SetPxPyPzE(T_Muon_Px->at(iMuon), T_Muon_Py->at(iMuon),
                 T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  if (lep.Pt() < minPt)            return false;
  if (TMath::Abs(lep.Eta()) > 2.4) return false;
  
  // POG Tight Muons definition              
  if (!T_Muon_IsGlobalMuon->at(iMuon)                            ) return false;
  if (!T_Muon_IsPFMuon->at(iMuon)                                ) return false;
  if (T_Muon_NormChi2GTrk->at(iMuon)                       >= 10.) return false;
  if (T_Muon_NValidHitsGTrk->at(iMuon)                     <= 0  ) return false;
  if (T_Muon_NumOfMatchedStations->at(iMuon)               <= 1  ) return false;                     
  if (TMath::Abs(T_Muon_dxyInTrack->at(iMuon))             >= 0.2) return false; 
  if (TMath::Abs(T_Muon_dzInTrack ->at(iMuon))             >= 0.5) return false;
  if (T_Muon_NLayers->at(iMuon)                            <= 5  ) return false;
  if (T_Muon_NValidPixelHitsInTrk->at(iMuon)               <= 0  ) return false;

  float relIso = muonIsolation(iMuon);
  
  if (relIso > 0.12) return false;
  
  return true;
}

bool ConverterPHYS14::isVetoMuon( unsigned iMuon, float minPt) 
{
  TLorentzVector lep;
  lep.SetPxPyPzE(T_Muon_Px->at(iMuon), T_Muon_Py->at(iMuon),
                 T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  if (lep.Pt() < minPt)            return false;
  if (TMath::Abs(lep.Eta()) > 2.4) return false;
  float relIso = muonIsolation(iMuon);
  if (relIso > 0.20)               return false; 
  if (T_Muon_IsPFMuon->at(iMuon) == 0) return false;
  if (T_Muon_IsGlobalMuon->at(iMuon) == 0 && T_Muon_IsTrackerMuonArbitrated->at(iMuon) == 0) return false; 
  return true;
}

float ConverterPHYS14::muonIsolation(unsigned iMuon){

  TLorentzVector lep;
  lep.SetPxPyPzE(T_Muon_Px->at(iMuon), T_Muon_Py->at(iMuon),
                 T_Muon_Pz->at(iMuon), T_Muon_Energy->at(iMuon));
  
  return (T_Muon_chargedHadronIsoR04->at(iMuon) + std::max(0.0 , T_Muon_neutralHadronIsoR04->at(iMuon) + T_Muon_photonIsoR04->at(iMuon)- 0.5*T_Muon_sumPUPtR04->at(iMuon)))/lep.Pt();
}


bool ConverterPHYS14::isTightElec(unsigned int iElec, float minPt){
  
  TLorentzVector lep;
  lep.SetPxPyPzE( T_Elec_Px->at(iElec), T_Elec_Py->at(iElec),
                  T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));
  
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if (sceta > 1.4442 && sceta < 1.566) return false;
  if (lep.Pt() < minPt)                return false;
  if (TMath::Abs(lep.Eta()) > 2.5)     return false;

  float relIso =  elecIsolation(iElec);

  // Tight ID requirements
  bool passTightID = false;
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.010181 && //
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.006574 && //
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.022868 && //       
       T_Elec_HtoE->at(iElec)                               < 0.037553 && //
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) -
                   T_Elec_eSuperClusterOverP->at(iElec) /
                   T_Elec_ecalEnergy->at(iElec))            < 0.131191 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.009924 && //
       TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.015310 &&
       T_Elec_nLost->at(iElec)                              <= 1       && //
       T_Elec_passConversionVeto->at(iElec)                 > 0        &&
       relIso                                               < 0.074355)
      passTightID = true;
  }
  else {
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.028766 && //
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.005681 && //
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.032046 && //       
       T_Elec_HtoE->at(iElec)                               < 0.081902 && //
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) -
                   T_Elec_eSuperClusterOverP->at(iElec) /
                   T_Elec_ecalEnergy->at(iElec))            < 0.106055 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.027261 && //
       TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.147154 &&
       T_Elec_nLost->at(iElec)                              <= 1       && //
       T_Elec_passConversionVeto->at(iElec)                 > 0        &&
       relIso                                               < 0.090185)
       passTightID = true;
  }  
  
  // Medium ID requirements
  bool passMediumID = false;
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.010399 && 
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.007641 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.032643 &&    
       T_Elec_HtoE->at(iElec)                               < 0.060662 && 
       relIso                                               < 0.097213 &&
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) - 
                  T_Elec_eSuperClusterOverP->at(iElec) / 
                  T_Elec_ecalEnergy->at(iElec))             < 0.153897 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.011811 && 
       //TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.070775 &&
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.070775 &&
       T_Elec_nLost->at(iElec)                              <= 1       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
      passMediumID = true;
  }
  else {
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.029524 && 
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.009285 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.042447 &&    
       T_Elec_HtoE->at(iElec)                               < 0.104263 && 
       relIso                                               < 0.116708 &&
       TMath::Abs( 1.0/T_Elec_ecalEnergy->at(iElec) - 
       T_Elec_eSuperClusterOverP->at(iElec) /
       T_Elec_ecalEnergy->at(iElec))                        < 0.137468 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.051682 && 
       //TMath::Abs(T_Elec_vz->at(iElec) - T_Vertex_z->at(0)) < 0.180720 &&
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.180720 &&
       T_Elec_nLost->at(iElec)                              <= 1       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
       passMediumID = true;
  }

  if (!passTightID) {}
  if (!passMediumID) return false;
   
  return true;
}


bool ConverterPHYS14::isVetoElec(unsigned int iElec, float minPt){
  
  TLorentzVector lep;
  lep.SetPxPyPzE( T_Elec_Px->at(iElec), T_Elec_Py->at(iElec),
                  T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));
  
  float sceta = TMath::Abs(T_Elec_SC_Eta->at(iElec));
  if (sceta > 1.4442 && sceta < 1.566) return false;
  if (lep.Pt() < minPt)                return false;
  if (TMath::Abs(lep.Eta()) > 2.5)     return false;

  float relIso =  elecIsolation(iElec);

  bool passVetoID = false;
  // veto ID requirements
  if(TMath::Abs(T_Elec_SC_Eta->at(iElec)) < 1.479) {  //Barrel electron
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.011100 &&
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.016315 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.252044 &&       
       T_Elec_HtoE->at(iElec)                               < 0.345843 && 
       relIso                                               < 0.164369 &&
       ( 1.0/T_Elec_ecalEnergy->at(iElec) -
         T_Elec_eSuperClusterOverP->at(iElec) /
         T_Elec_ecalEnergy->at(iElec))                      < 0.248070 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.060279 && 
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.800538 &&
       T_Elec_nLost->at(iElec)                              <= 2       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
      passVetoID = true;
  }
  else {
    if(T_Elec_sigmaIetaIetaFull5by5->at(iElec)              < 0.033987 && 
       TMath::Abs(T_Elec_deltaEtaIn->at(iElec))             < 0.010671 && 
       TMath::Abs(T_Elec_deltaPhiIn->at(iElec))             < 0.245263 &&    
       T_Elec_HtoE->at(iElec)                               < 0.134691 && 
       relIso                                               < 0.212604 &&
       ( 1.0/T_Elec_ecalEnergy->at(iElec) -
         T_Elec_eSuperClusterOverP->at(iElec) /
         T_Elec_ecalEnergy->at(iElec))                      < 0.157160 &&
       TMath::Abs(T_Elec_IPwrtPV->at(iElec))                < 0.273097 && 
       TMath::Abs(T_Elec_dzwrtPV->at(iElec))                < 0.885860 &&
       T_Elec_nLost->at(iElec)                              <= 3       && 
       T_Elec_passConversionVeto->at(iElec)                 > 0 )
       passVetoID = true;
  }  
  
  if (!passVetoID) return false;
  
  return true;
}


float ConverterPHYS14::elecIsolation(unsigned iElec) {

  TLorentzVector lep( T_Elec_Px->at(iElec), T_Elec_Py->at(iElec),
                      T_Elec_Pz->at(iElec), T_Elec_Energy->at(iElec));

  float pt     = lep.Pt();
  float relIso = (T_Elec_sumChargedHadronPt->at(iElec) +
                  std::max(0.0 , T_Elec_sumNeutralHadronEt->at(iElec) +
                                 T_Elec_sumPhotonEt->at(iElec) -
                                 0.5*T_Elec_sumPUPt->at(iElec)))/pt;
  return relIso;

}

bool ConverterPHYS14::passJetID(unsigned iJet) {

  if ( !(T_JetAKCHS_nDaughters->at(iJet)        > 1   ) ) return false;
  if ( !(T_JetAKCHS_NeutHadEnergyFrac->at(iJet) < 0.99) ) return false;
  if ( !(T_JetAKCHS_NeutEmEnergyFrac ->at(iJet) < 0.99) ) return false;
  if (TMath::Abs(T_JetAKCHS_Eta->at(iJet)) < 2.5){
  if ( !(T_JetAKCHS_CharEmEnergyFrac->at(iJet)  < 0.99) ) return false;
  if ( !(T_JetAKCHS_CharHadEnergyFrac->at(iJet) > 0.00) ) return false;
  if ( !(T_JetAKCHS_ChargedMultiplicity->at(iJet) > 0 ) ) return false;
  return true;
  }
  return true;
}


float ConverterPHYS14::getMT2( const mut::Lepton & lep0, const mut::Lepton & lep1, const mut::MET & met) {
    double pa[3];
    double pb[3];
    double pmiss[3];

    pmiss[0] = 0.; // not required
    pmiss[1] = met.px();
    pmiss[2] = met.py();

    pa[0] = 0.;
    pa[1] = lep0.px();
    pa[2] = lep0.py();

    pb[0] = 0.;
    pb[1] = lep1.px();
    pb[2] = lep1.py();

    mt2bisect* MT2bisect = new mt2bisect();
    MT2bisect->set_verbose(0);
    MT2bisect->set_momenta(pa, pb, pmiss);
    MT2bisect->set_mn(0.); // test mass
    float mt2 = MT2bisect->get_mt2();
    delete MT2bisect;
    return mt2;
}
