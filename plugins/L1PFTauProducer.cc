// -*- C++ -*-
//
// Package:    L1PFTauProducer
// Class:      L1PFTauProducer
// 
/**\class L1PFTauProducer L1PFTauProducer.cc L1Trigger/phase2Demonstrator/plugins/L1PFTauProducer.cc

 Description: Level 1 L1PFTaus for the Demonstrator

 Implementation:
     [Notes on implementation]
*/


/*Implement Position Diff for Crystals
 *Position Diff for iEta from RECO
 *Position Diff for iPhi from RECO
 *Put new methods in header file
 *Make the couts for the number of taus per event
 *Make the couts for the content of the tau
 *Merge the electron Grid and the Photon Grid
 * --> Need a way to add electrons which are in the same grid position
 * IMPLEMENT THE PHOTON AND E/G FINDING IN CLUSTERS
 *
 * Finish the isolation variable calculation for type of object
 */

#include "L1Trigger/phase2Demonstrator/interface/L1PFTauProducer.hh"

L1PFTauProducer::L1PFTauProducer(const edm::ParameterSet& cfg) :
  debug(                cfg.getUntrackedParameter<bool>("debug", false)),
  input_EoH_cut_(       cfg.getUntrackedParameter<int>("EoH_cut", 2)), // LSB is 0.1 so 2 corresonds to 0.2
  input_HoE_cut_(       cfg.getUntrackedParameter<int>("HoE_cut", 8)), // LSB is 0.1 so 8 corresonds to 0.8
  input_min_n_stubs_(   cfg.getUntrackedParameter<int>("min_n_stubs", 4)), // LSB is 0.1 so 8 corresonds to 0.8
  input_max_chi2_(      cfg.getUntrackedParameter<int>("max_chi2", 15)), // LSB is 0.1 so 8 corresonds to 0.8
  //input_pt_seed_(       cfg.getUntrackedParameter<double>("pt_seed", 0.2)), // LSB is 0.1 so 8 corresonds to 0.8
  three_prong_delta_r_( cfg.getUntrackedParameter<double>("three_prong_dr", 0.2)), // LSB is 0.1 so 8 corresonds to 0.8
  three_prong_max_delta_Z_( cfg.getUntrackedParameter<double>("three_prong_max_dZ", 0.3)), // LSB is 0.1 so 8 corresonds to 0.8
  isolation_delta_r_(   cfg.getUntrackedParameter<double>("iso_dr", 0.5)), // LSB is 0.1 so 8 corresonds to 0.8
  L1TrackInputTag(      cfg.getParameter<edm::InputTag>("L1TrackInputTag")),
  L1ClustersToken_(     consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Clusters"))),
  ttTrackToken_(        consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag)),
  L1PFToken_(           consumes< L1PFObjectCollection >(cfg.getParameter<edm::InputTag>("L1PFObjects"))),
  L1NeutralToken_(      consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Neutrals")) )
{
  //produces three collections of taus, one that uses full FW, one that uses only L1 objects and one that uses only reco objects
  produces< L1PFTauCollection >( "L1PFTaus" ).setBranchAlias("L1PFTaus");

}

void L1PFTauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::unique_ptr<L1PFTauCollection> newL1PFTauCollection(new L1PFTauCollection);

  using namespace edm;
  edm::Handle< L1CaloClusterCollection > l1NeutralClusters;
  iEvent.getByToken( L1NeutralToken_, l1NeutralClusters);

  edm::Handle< L1CaloClusterCollection > l1CaloClusters;
  iEvent.getByToken( L1ClustersToken_, l1CaloClusters);

  edm::Handle< std::vector<L1PFObject> > l1PFChargedCandidates;
  iEvent.getByToken( L1PFToken_, l1PFChargedCandidates);

  L1PFTau tauCands[12]; //6 tau cands possible, each with as many as 3 signal candidates
  L1PFObject electronGrid[12][5][5];

  unsigned int nTauCands = 0;

  algo_config_t algo_config;
  
  algo_config.three_prong_seed = 7;	 
  algo_config.three_prong_delta_r = 0.1;	 
  algo_config.isolation_delta_r = 0.5; 
  algo_config.one_prong_seed = 15;
  algo_config.input_EoH_cut = 10;
  algo_config.max_neighbor_strip_dist = 0.1; 
  algo_config.min_strip = 2.5;
  algo_config.eg_strip_merge = 2.5;          

  //Find taus using charged hadrons
  for(unsigned int iCand = 0; iCand < l1PFChargedCandidates->size(); iCand++){

    //Find Tau Seed
    if(l1PFChargedCandidates->at(iCand).isChargedHadron() == true && l1PFChargedCandidates->at(iCand).p4().Pt() > 2){
      //std::cout<<"New Charged Hadron Candidate Search"<<std::endl;
      L1PFObject tau_cand[3];
      L1PFObject pfChargedHadron_Seed = l1PFChargedCandidates->at(iCand);
      if(debug)
	std::cout<<"PFChargedHadronSeed Tower Eta: "<<pfChargedHadron_Seed.towerEta()<<
	  " phi: "<< pfChargedHadron_Seed.towerPhi()<<" etaside: "<<pfChargedHadron_Seed.towerEtaSide() <<std::endl;
      tau_cand[0] = pfChargedHadron_Seed;
      int n_prongs_found = 1;
      float isolationSum = 0;
      L1PFObject electronGrid[5][5];
      float isolationSumChargedHadron = 0;
      if(debug)
	std::cout<<"Seed Cand Pt: "<<pfChargedHadron_Seed.p4().Pt() <<
	  "    Eta: "<<pfChargedHadron_Seed.p4().Eta() <<
	  "    Phi: "<<pfChargedHadron_Seed.p4().Phi()<< 
	  "     z0: "<< pfChargedHadron_Seed.z0() << 
	  " chi2_4: "<< pfChargedHadron_Seed.trkChi2_4()<< 
	  " chi2_5: "<<  pfChargedHadron_Seed.trkChi2_5() <<
	  " chi2_Red: "<<  pfChargedHadron_Seed.trkChi2_Red() <<
	  " nStubs: "<<  pfChargedHadron_Seed.nStubs() <<std::endl;

      if(pfChargedHadron_Seed.nStubs() < input_min_n_stubs_)
	continue;

      if(pfChargedHadron_Seed.trkChi2_Red() > input_max_chi2_)
	continue;

      for(unsigned int jCand = 0; jCand < l1PFChargedCandidates->size(); jCand++){

	if(jCand <= iCand)
	  continue;

	if(l1PFChargedCandidates->at(jCand).isChargedHadron()){

	  L1PFObject pfChargedHadron_signal_cand = l1PFChargedCandidates->at(jCand);

	  if(debug)
	    std::cout<<"   charged Hadron Found, pt: "<< pfChargedHadron_signal_cand.p4().Pt() <<
	      " eta: "<< pfChargedHadron_signal_cand.p4().Eta() <<
	      " phi: "<< pfChargedHadron_signal_cand.p4().Phi() <<
	      " Delta R: "<< (fabs(pfChargedHadron_signal_cand.p4().Eta()-pfChargedHadron_Seed.p4().Eta())
			      +fabs(pfChargedHadron_signal_cand.p4().Phi()-pfChargedHadron_Seed.p4().Phi()) ) << 
	      " Delta Z: "<< (pfChargedHadron_Seed.z0() - pfChargedHadron_signal_cand.z0()) <<
	      " chi2_4: "<< pfChargedHadron_signal_cand.trkChi2_4() << 
	      " chi2_5: "<< pfChargedHadron_signal_cand.trkChi2_5() << 
	      " chi2_Red: "<< pfChargedHadron_signal_cand.trkChi2_Red() <<
	      "   nStubs: "<< pfChargedHadron_signal_cand.nStubs() << std::endl;
	  

	  if(pfChargedHadron_signal_cand.nStubs() < input_min_n_stubs_)
	    continue;
	  
	  if(pfChargedHadron_signal_cand.trkChi2_Red() > input_max_chi2_)
	    continue;

	  float pfCharged_pt = pfChargedHadron_signal_cand.p4().pt();
	  float deltaR = 0.2;
	  if(pfCharged_pt>3){
	    deltaR = 0.18;
	  }
	  if(pfCharged_pt>5){
	    deltaR = 0.17;
	  }
	  if(pfCharged_pt>10){
	    deltaR = 0.1;
	  }
	  if(pfCharged_pt>18){
	    deltaR = 0.06;
	  }
	  if(pfCharged_pt>25){
	    deltaR = 0.03;
	  }
	  
	  
	  if(Delta_R(pfChargedHadron_signal_cand.p4().Eta(), pfChargedHadron_signal_cand.p4().Phi(), 
		     pfChargedHadron_Seed.p4().Eta(), pfChargedHadron_Seed.p4().Phi(), 
		     deltaR)
	     && fabs(pfChargedHadron_Seed.z0() - pfChargedHadron_signal_cand.z0()) < three_prong_max_delta_Z_){


	    n_prongs_found++;

	  if(n_prongs_found==2){
	    tau_cand[1] = pfChargedHadron_signal_cand;
	  }

	  if(n_prongs_found==3){
	    tau_cand[2] = pfChargedHadron_signal_cand;
	    //std::cout<<"3 charged hadrons found"<<std::endl;
	    if(debug)
	      std::cout<<"prong 1: "<<tau_cand[0].p4().Pt() <<
		" prong 2: " << tau_cand[1].p4().Pt() <<
		" prong 3: " << tau_cand[2].p4().Pt() << std::endl;
	  }
	  
	  }// Close If Greater than Three_prong_delta_r
	  else if(Delta_R(pfChargedHadron_signal_cand.p4().Eta(), pfChargedHadron_signal_cand.p4().Phi(), 
			  pfChargedHadron_Seed.p4().Eta(),        pfChargedHadron_Seed.p4().Phi(), 
			  isolation_delta_r_)){

	    // Sum the isolation here 
	    isolationSumChargedHadron += pfChargedHadron_signal_cand.p4().Pt();

	  }// Less than isolation_delta_r
	}// Is Charged Hadron

	/*
	else if(l1PFChargedCandidates->at(jCand).isElectron()){
	  
	  L1PFObject electronCand;
	  if(Delta_R(electronCand.p4().Eta(),         electronCand.p4().Phi(), 
		     pfChargedHadron_Seed.p4().Eta(), pfChargedHadron_Seed.p4().Phi(), 
		     isolation_delta_r_)){
	    //IMPLEMENT ME
	    unsigned int index_eta = ieta_diff(l1PFChargedCandidates->at(jCand),  pfChargedHadron_Seed);
	    unsigned int index_phi = iphi_diff(l1PFChargedCandidates->at(jCand),  pfChargedHadron_Seed);
	    electronGridTemp[index_eta][index_phi] = l1PFChargedCandidates->at(jCand);
	  }
	}
      }//Finished finding signal cands and creating Isolation and electron grid for strip finding
	*/
      }
      // Create the 1 Prong Taus
      if(n_prongs_found==1||n_prongs_found==2){
	if(pfChargedHadron_Seed.p4().Pt() > 2){
	  if(n_prongs_found==2){
	    //Add second prong to Isolation Cone
	    isolationSumChargedHadron += tau_cand[1].p4().Pt();	    
	  }
	  L1PFTau tempL1PFTau;

	  math::PtEtaPhiMLorentzVector tempP4(pfChargedHadron_Seed.p4().Pt(), 
					      pfChargedHadron_Seed.p4().Eta(), 
					      pfChargedHadron_Seed.p4().Phi(), 
					      pfChargedHadron_Seed.p4().Pt());
	  tempL1PFTau.setP4(tempP4);

	  tempL1PFTau.setTowerEtaSide(pfChargedHadron_Seed.towerEtaSide());
	  tempL1PFTau.setTowerEta(pfChargedHadron_Seed.towerEta());
	  tempL1PFTau.setTowerPhi(pfChargedHadron_Seed.towerPhi());
	  tempL1PFTau.setTauType(0);
	  tempL1PFTau.setRawIso(isolationSumChargedHadron);
	  tempL1PFTau.setRelIso(isolationSumChargedHadron/pfChargedHadron_Seed.p4().Pt());
	  tempL1PFTau.setChargedIso(isolationSumChargedHadron);
	  tempL1PFTau.setNeutralIso(0);

	  //electronGrid[nTauCands] = electronGridTemp;
	  tauCands[nTauCands] = tempL1PFTau;
	  nTauCands++;

	}
      }

      // Create the 3 Prong Taus
      if(n_prongs_found>2){
	L1PFTau tempL1PFTau;
	float totalPT = tau_cand[0].p4().Pt() + tau_cand[1].p4().Pt() + tau_cand[2].p4().Pt() ;
	float averageEta = (tau_cand[0].p4().Pt()*tau_cand[0].p4().Eta() + tau_cand[1].p4().Pt()*tau_cand[1].p4().Eta() + tau_cand[2].p4().Pt()*tau_cand[2].p4().Eta())/totalPT ;
	float averagePhi = (tau_cand[0].p4().Pt()*tau_cand[0].p4().Phi() + tau_cand[1].p4().Pt()*tau_cand[1].p4().Phi() + tau_cand[2].p4().Pt()*tau_cand[2].p4().Phi())/totalPT ;
	if(debug)
	  std::cout<<"NProngs Found > 2: total PT: "<<totalPT<<" averageEta: "<<averageEta<<" averagePhi: "<<averagePhi<< " tau number: "<< nTauCands<<std::endl<<std::endl;

	math::PtEtaPhiMLorentzVector tempP4(totalPT, 
					    averageEta, 
					    averagePhi, 
					    totalPT);
	tempL1PFTau.setP4(tempP4);
	tempL1PFTau.setTowerEtaSide(pfChargedHadron_Seed.towerEtaSide());
	tempL1PFTau.setTowerEta(pfChargedHadron_Seed.towerEta());
	tempL1PFTau.setTowerPhi(pfChargedHadron_Seed.towerPhi());
	tempL1PFTau.setTauType(10);
	tempL1PFTau.setRawIso(isolationSumChargedHadron);
	tempL1PFTau.setChargedIso(isolationSumChargedHadron);
	tempL1PFTau.setNeutralIso(0);
	tempL1PFTau.setRelIso(isolationSumChargedHadron/totalPT);

	//electronGrid[nTauCands] = electronGridTemp;	
	tauCands[nTauCands] = tempL1PFTau;
	nTauCands++;
      }
    }
    if(nTauCands > 11){
      	if(debug)
	  std::cout<<"nTauCands >12 breaking"<<std::endl;
      break;
    }
  }
  // For each 1prong Tau Candidate Create the Strip by grabbing the index of the nearby neutrals

  //find pi0's using electrons/gamma
  //stitch the strips by combining electron and gamma candidates
  for(unsigned int i = 0; i < nTauCands; i++){

    if(debug)
      std::cout<<"Performing strip algorithm tau cand: "<<i<<std::endl;

    if(tauCands[i].tauType()==10)
      continue;

     //needs tower eta and tower phi, as well as pt, eta, phi
     //strip_t tempStrip[5];    
     //strip_t finalStrip;
     //float crystal_distance_ = 2;
     pf_charged_t electron_grid[5][5];
     strip_alg(tauCands[i], electron_grid, l1CaloClusters, algo_config);

  }
  //now fill the new collection
   for(unsigned int i = 0; i < 12; i++){
     newL1PFTauCollection->push_back(tauCands[i]);
     if(debug){
       if(tauCands[i].p4().Et()>0){
	 std::cout<<"------ SUMMARY OF THE TAU CANDS ------"<<std::endl;
	 std::cout<<tauCands[i]<<std::endl;
	}
    }
  }


  iEvent.put( std::move(newL1PFTauCollection) , "L1PFTaus" );

}
/*
void Check_And_Merge_Strip(L1CaloCluster cluster1, L1CaloCluster cluster2, strip_t &tempStrip, int crystal_distance){
  if( cluster1.isPhoton() && cluster2.isPhoton() 
      && Position_Diff( cluster1, cluster2 ) < crystal_distance 
      && tempStrip.Pt() < (cluster1.p4().Pt() + cluster2.p4().Pt())){

    tempStrip.setEt(   cluster1.p4().Pt()                     + cluster2.p4().Pt()                                      );
    tempStrip.setEta( (cluster1.p4().Eta()*cluster1.p4().Pt() + cluster2.p4().Eta()*cluster2.p4().Pt()) / tempStrip.et());
    tempStrip.setPhi( (cluster1.p4().Phi()*cluster1.p4().Pt() + cluster2.p4().Phi()*cluster2.p4().Pt()) / tempStrip.et());

  }
}
*/

// Note ieta/iphi add is 0 by default
uint32_t L1PFTauProducer::findTheIndexFromReco( float eta, float phi, int iEta_add, int iPhi_add ){
  //uint32_t PFObjectProducer::findTheIndexFromReco(float trackET, int charge, float eta, float phi){
 
  uint32_t index;
  triggerGeometryTools tool;
  int iEta = tool.convertGenEta(eta) + iEta_add;
  int iPhi = tool.convertGenPhi(phi) + iPhi_add;

  if(eta < 0)
    index = (uint32_t) ((20 - iEta)*72 + iPhi - 2);
  else
    index = (uint32_t) (( iEta )*72 + iPhi) ;

  /*
  if(trackET<15 && charge == -1){
  if(eta < 0)
    index = (uint32_t) ((20 - iEta)*72 + iPhi - 2);
  else
    index = (uint32_t) (( iEta )*72 + iPhi) ;
  }
  */
  //std::cout<<"tower iEta "<<iEta <<" iPhi "<<iPhi<<" index: "<<index<<std::endl;
  //iEta 0 iPhi 57 index: 7257
  return index;
}


bool L1PFTauProducer::Delta_R( float eta1, float phi1, float eta2, float phi2, float max_dr)
{
  if(fabs(eta1-eta2)+fabs(phi1-phi2) < max_dr ){
    return true;
  }
  else{
    return false;
  }
}


// Check if each tau_phi slice is strip like
// This module looks at each strip starting from the top cluster and going to the bottom
//  -> Geometry outlined above for index
// grid for strip creation is below
// strips can be at most two towers in phi, potentially 2 towers in eta if crystals pass neighbors requirements.
//   -----------
//   |0|0|0|0|0|
//   |1|1|1|1|1|
//   |2|2|2|2|2|
//   |3|3|3|3|3|
//   |4|4|4|4|4|
//   -----------
void L1PFTauProducer::strip_alg(pftau_t &tau_cand, pf_charged_t electron_grid[5][5],  edm::Handle<std::vector<L1CaloCluster> >&neutral_clusters, algo_config_t algo_config){

  if(debug)
    std::cout<<"Tau Cand Tower Eta: "<<tau_cand.towerEta()<<" Phi: "<<tau_cand.towerEta()<<std::endl;

  //electron_grid[phi][eta] (equivalent to [y][x], typically loop over second index and then first)
  if(debug){
    std::cout<<"electron_grid "<<std::endl;
    for(int i = 0; i<5; i++){
      for(int j = 0; j<5; j++){
	std::cout<<electron_grid[i][j].p4().pt()<<" ";
      }
      std::cout<<std::endl;
    }
  }

  float tau_eta      = tau_cand.p4().Eta();
  //float tau_eta_side = tau_cand.eta_side;
  float tau_phi      = tau_cand.p4().Phi();
  unsigned tau_teta = tau_cand.towerEta();
  unsigned tau_tphi = tau_cand.towerPhi();
  int tau_etaside   = tau_cand.towerEtaSide();

  if(debug)
    std::cout<<"Tau Tower Eta: "<<tau_teta<<" phi: "<<tau_tphi<<" tau_etaside: "<<tau_etaside<<std::endl;
  //float index[5];
  cluster_t neutral_cluster_grid[5][5];

  //cluster_grid[phi][eta] (equivalent to [y][x], typically loop over second index and then first)
  /*
  for(unsigned int i = 0; i<5; i++){
    for(unsigned int j = 0; j<5; j++){
      float grid_eta = tau_eta + j*0.087 - 0.087*2; // subtract two clusters to center the grid
      float grid_phi = tau_phi + i*0.087 - 0.087*2; // subtract two clusters to center the grid
      //neutral_cluster_grid[i][j] = find_matched_cluster(neutral_clusters, grid_eta, grid_phi);
      neutral_cluster_grid[i][j] = find_matched_cluster_tower(neutral_clusters, grid_eta, grid_phi);

    }
  }
  */
  for( int i = -2; i<3; i++){
    for( int j = -2; j<3; j++){
      int grid_eta = tau_teta + i; // subtract two clusters to center the grid
      int grid_phi = tau_tphi + j; // subtract two clusters to center the grid
      int grid_etaside = tau_etaside;
      if(grid_eta < 0){
	grid_etaside = 0;
	grid_eta = abs(grid_eta);
      }

      if(grid_phi < 0){
	grid_phi = 72 + tau_tphi - i;
      }

      //neutral_cluster_grid[i][j] = find_matched_cluster(neutral_clusters, grid_eta, grid_phi);
      neutral_cluster_grid[i+2][j+2] = find_matched_cluster_tower(neutral_clusters, (unsigned)grid_eta, grid_phi, grid_etaside);

    }
  }

  if(debug){
    std::cout<<"neutral_grid "<<std::endl;
    for(int i = 0; i<5; i++){
      for(int j = 0; j<5; j++){
	std::cout<<neutral_cluster_grid[i][j].ecalEnergy()<<" ";
      }
      std::cout<<std::endl;
    }
  }


  strip_t temp_strip[5];    
  strip_t final_strip;
  cluster_t cluster;
  pf_charged_t charged;

  //float index1; float index2;

  for(unsigned int i = 0; i<4; i++){
    for(unsigned int j = 0; j<5; j++){
      unsigned int ip = i+1;
      ///////// merge_strip_algo //////////
      merge_strip_algo(neutral_cluster_grid[i][j], electron_grid[i][j], neutral_cluster_grid[ip][j], electron_grid[ip][j], temp_strip[j], algo_config);
      if(debug)
	std::cout<<"cluster i: "<< i <<" j: "<<j<<" et: "<<neutral_cluster_grid[i][j].ecalEnergy()<<" et2: "<< neutral_cluster_grid[ip][j].ecalEnergy()<<" tempStrip et: "<< temp_strip[i].et <<std::endl;
      //std::cout<<"temp_strip"
    }
  }

  final_strip = temp_strip[0];
  for(unsigned int j = 1; j < 5; j++){
    
    //first check if strip j is greater than final strip
    if(temp_strip[j].et > final_strip.et){
      final_strip = temp_strip[j];
    }

    if(debug)
      if(temp_strip[j].et>0 && temp_strip[j-1].et>0){
	std::cout<<"Check Strip To Merge: et1: "<<temp_strip[j].et<<" et2: "<<temp_strip[j-1].et<< " deltaR: " << delta_r_strip( temp_strip[j], temp_strip[j-1])<<std::endl;
      }

    //merge strips if two are close by
    if(delta_r_strip( temp_strip[j], temp_strip[j-1]) < algo_config.max_neighbor_strip_dist){
      if(debug)
	std::cout<<"Strips are close by et1: "<<temp_strip[j].et<<" et2: "<<temp_strip[j-1].et<<std::endl;
      if(temp_strip[j].et + temp_strip[j-1].et > final_strip.et){

	float et  = temp_strip[j].et + temp_strip[j-1].et;
	final_strip.et  = et ;
	final_strip.eta = weighted_avg_eta_s_s(temp_strip[j], temp_strip[j-1]);
	final_strip.phi = weighted_avg_phi_s_s(temp_strip[j], temp_strip[j-1]);

      }
    }
  }

  if(debug)
    std::cout<<"final strip et: "<<final_strip.et<<" eta: "<<final_strip.eta<<" phi: "<<final_strip.phi<<std::endl;
  //find 1 prong pi0's by combining the 1 prong and the pi0's
  if(tau_cand.tauType() == 0 && final_strip.et > algo_config.min_strip){
    //take care of eta_side
    float temp_tau_et;
    float temp_tau_eta;
    float temp_tau_phi;
    temp_tau_et  = tau_cand.p4().Pt() + final_strip.et;
    temp_tau_eta = weighted_avg_eta_t_s(tau_cand,final_strip);
    temp_tau_phi = weighted_avg_phi_t_s(tau_cand,final_strip);

    math::PtEtaPhiMLorentzVector tempP4(temp_tau_et,temp_tau_eta,temp_tau_phi,temp_tau_et);
    tau_cand.setP4(tempP4);

    //tau_cand.et       = temp_tau_et;
    //tau_cand.eta      = temp_tau_eta;
    //tau_cand.phi      = temp_tau_phi;
    tau_cand.setTauType(1);
  }

  if(tau_cand.tauType() == 10 && final_strip.et > algo_config.min_strip){
    //take care of eta_side
    float temp_tau_et;
    float temp_tau_eta;
    float temp_tau_phi;
    temp_tau_et  = tau_cand.p4().Pt() + final_strip.et;
    temp_tau_eta = weighted_avg_eta_t_s(tau_cand,final_strip);
    temp_tau_phi = weighted_avg_phi_t_s(tau_cand,final_strip);

    math::PtEtaPhiMLorentzVector tempP4(temp_tau_et,temp_tau_eta,temp_tau_phi,temp_tau_et);
    tau_cand.setP4(tempP4);

    //tau_cand.et       = temp_tau_et;
    //tau_cand.eta      = temp_tau_eta;
    //tau_cand.phi      = temp_tau_phi;
    tau_cand.setTauType(11);
  }
}


cluster_t L1PFTauProducer::find_matched_cluster(edm::Handle<std::vector<L1CaloCluster> >&neutral_clusters, float eta, float phi){

  cluster_t cluster_to_return;
  //cluster_to_return.setPtEtaPhiE(0,0,0,0);
  cluster_to_return.setEcalEnergy(0);
  cluster_to_return.setHcalEnergy(0);

  float minDistance = 10;
  float minPt = 0;
  for(unsigned int i = 0; i < neutral_clusters->size(); i++){
    float tempDistance = fabs(neutral_clusters->at(i).p4().Eta()-eta) + fabs(neutral_clusters->at(i).p4().Phi() - phi);
    float tempPt = neutral_clusters->at(i).ecalEnergy();
 
    if(fabs(neutral_clusters->at(i).p4().Eta() - eta)< 0.087/2 && 
       fabs(neutral_clusters->at(i).p4().Phi() - phi)< 0.087/2 && 
       tempDistance < minDistance && tempPt >= minPt ){
	 
	 minPt = tempPt;
	 minDistance = tempDistance;
	 cluster_to_return = neutral_clusters->at(i);
	 
	 if(debug)
	   if(neutral_clusters->at(i).ecalEnergy()>0){
	     std::cout<<"Matched ecal cluster et greater than 0: "<<neutral_clusters->at(i).ecalEnergy()<<std::endl;
	     //std::cout<<"   Matched Cluster Tower Eta: "<< neutral_clusters->at(i).towerEta()<<" Tower phi: "<<neutral_clusters->at(i).towerPhi()<<std::endl;

	     unsigned cluster_teta = neutral_clusters->at(i).towerEta();
	     unsigned cluster_tphi = neutral_clusters->at(i).towerPhi();
	     int cluster_etaside   = neutral_clusters->at(i).towerEtaSide();
	     
	     if(debug)
	       std::cout<<"Cluster Tower Eta: "<<cluster_teta<<" phi: "<<cluster_tphi<<" cluster_etaside: "<<cluster_etaside<<std::endl;
	   }
       }
    }
    return cluster_to_return;
    
}


cluster_t L1PFTauProducer::find_matched_cluster_tower(edm::Handle<std::vector<L1CaloCluster> >&neutral_clusters, unsigned eta, unsigned phi, int eta_side){

  cluster_t cluster_to_return;
  //cluster_to_return.setPtEtaPhiE(0,0,0,0);
  cluster_to_return.setEcalEnergy(0);
  cluster_to_return.setHcalEnergy(0);

  float minDistance = 10;
  float minPt = 0;
  //if(debug)
    //std::cout<<"Find Matched Cluster teta: "<<eta<<" tphi: "<<phi<<" teta_side: "<<eta_side<<std::endl;
 
  for(unsigned int i = 0; i < neutral_clusters->size(); i++){

    unsigned cluster_teta = neutral_clusters->at(i).towerEta();
    unsigned cluster_tphi = neutral_clusters->at(i).towerPhi();
    int cluster_etaside   = neutral_clusters->at(i).towerEtaSide();

    //if(debug && neutral_clusters->at(i).et()>0)
      //std::cout<<"              cluster eta: "<<eta<<"  phi: "<<phi<<"  eta_side: "<<eta_side<<std::endl;
   
    if(cluster_teta == eta && cluster_tphi == phi && cluster_etaside == eta_side){
      cluster_to_return = neutral_clusters->at(i);
      if(debug)
	std::cout<<"Matched Cluster Tower Eta: "<<cluster_teta<<" phi: "<<cluster_tphi<<" cluster_etaside: "<<cluster_etaside<<std::endl;
      break;
   }
  }
  return cluster_to_return;
    
}


/* Start by checking if the new potetntial strip is greater in pt than the old and if both clusters are pi0
 * Next
 *
 */

void L1PFTauProducer::merge_strip_algo(cluster_t cluster_1, pf_charged_t electron_1, cluster_t cluster_2, pf_charged_t electron_2, strip_t &strip, algo_config_t algo_config){

  // Requirement that it is a photon is checking if E/H is high for the cluster
  if(//cluster_1.isPi0() > 0 && cluster_2.isPi0() > 0 &&
     delta_r_cluster(cluster_1, cluster_2) < algo_config.eg_strip_merge
     && strip.et < (cluster_1.ecalEnergy() + cluster_2.ecalEnergy() )){

    strip.et  = cluster_1.ecalEnergy() + cluster_2.ecalEnergy();
    strip.eta = weighted_avg_eta_c_c(cluster_1, cluster_2);
    strip.phi = weighted_avg_phi_c_c(cluster_1, cluster_2);
  } //if clusters ddon't merge then take the higher pt cluster if it is greater than current strip pt and it is pi0 like
  else if(strip.et < cluster_1.ecalEnergy() && cluster_1.ecalEnergy() > 0 && cluster_1.ecalEnergy() > cluster_2.ecalEnergy()){
    strip.et  = cluster_1.ecalEnergy();
    strip.eta = cluster_1.p4().Eta();
    strip.phi = cluster_1.p4().Phi();
  }    
  else if(strip.et < cluster_2.ecalEnergy() && cluster_2.ecalEnergy()>0 && cluster_2.ecalEnergy() > cluster_1.ecalEnergy()){
    strip.et  = cluster_2.ecalEnergy();
    strip.eta = cluster_2.p4().Eta();
    strip.phi = cluster_2.p4().Phi();
  }    
}
/*
void L1PFTauProducer::merge_strip_algo(cluster_t cluster_1, pf_charged_t electron_1, cluster_t cluster_2, pf_charged_t electron_2, strip_t &strip, algo_config_t algo_config){
  cluster_t temp_cluster_1;
  float temp1_pt  = cluster_1.p4().Et() + electron_1.p4().Et();
  float temp1_phi = weighted_avg_phi_c_p(cluster_1, electron_1);
  float temp1_eta = weighted_avg_eta_c_p(cluster_1, electron_1);
  temp_cluster_1.setPtEtaPhiE(temp1_pt, temp1_eta, temp1_phi, temp1_pt);

  cluster_t temp_cluster_2;
  float temp_2_pt  = cluster_2.p4().Et() + electron_2.p4().Et();
  float temp_2_phi = weighted_avg_phi_c_p(cluster_2, electron_2);
  float temp_2_eta = weighted_avg_eta_c_p(cluster_2, electron_2);
  temp_cluster_2.setPtEtaPhiE(temp_2_pt, temp_2_eta, temp_2_phi, temp_2_pt);

  // Requirement that it is a photon is checking if E/H is high for the cluster
  if(temp_cluster_1.isPi0() > 0 && temp_cluster_2.isPi0() > 0
     && delta_r_cluster(temp_cluster_1, temp_cluster_2) < algo_config.eg_strip_merge
     && strip.et < (temp_cluster_1.p4().Et() + temp_cluster_2.p4().Et())){
    
    strip.et  = temp_cluster_1.p4().Et() + temp_cluster_2.p4().Et();
    strip.eta = weighted_avg_eta_c_c(temp_cluster_1, temp_cluster_2);
    strip.phi = weighted_avg_phi_c_c(temp_cluster_1, temp_cluster_2);
  }
}
*/
float L1PFTauProducer::weighted_avg_phi_c_p(cluster_t cluster_1, pf_charged_t pf_charged_1){
  float total_pt = (cluster_1.p4().Pt()+pf_charged_1.p4().Pt());
  float avg_phi = (cluster_1.p4().Phi()*cluster_1.p4().Pt() + pf_charged_1.p4().Phi()*pf_charged_1.p4().Pt())/total_pt;
  return avg_phi;
}

float L1PFTauProducer::weighted_avg_eta_c_p(cluster_t cluster_1, pf_charged_t pf_charged_1){
  float total_pt = (cluster_1.p4().Pt()+pf_charged_1.p4().Pt());
  float avg_eta = (cluster_1.p4().Eta()*cluster_1.p4().Pt() + pf_charged_1.p4().Eta()*pf_charged_1.p4().Pt())/total_pt;
  return avg_eta;
}

float L1PFTauProducer::weighted_avg_phi_c_c(cluster_t cluster_1, cluster_t cluster_2){
  float total_pt = (cluster_1.p4().Pt()+cluster_2.p4().Pt());
  float avg_phi = (cluster_1.p4().Phi()*cluster_1.p4().Pt() + cluster_2.p4().Phi()*cluster_2.p4().Pt())/total_pt;
  return avg_phi;
}

float L1PFTauProducer::weighted_avg_eta_c_c(cluster_t cluster_1, cluster_t cluster_2){
  float total_pt = (cluster_1.p4().Pt()+cluster_2.p4().Pt());
  float avg_eta = (cluster_1.p4().Eta()*cluster_1.p4().Pt() + cluster_2.p4().Eta()*cluster_2.p4().Pt())/total_pt;
  return avg_eta;
}

float L1PFTauProducer::weighted_avg_phi_s_s(strip_t strip_1, strip_t strip_2){
  float total_pt = (strip_1.et+strip_2.et);
  float avg_phi  = (strip_1.phi*strip_1.et + strip_2.phi*strip_2.et)/total_pt;
  return avg_phi;
}

float L1PFTauProducer::weighted_avg_eta_s_s(strip_t strip_1, strip_t strip_2){
  float total_pt = (strip_1.et+strip_2.et);
  float avg_eta = (strip_1.eta*strip_1.et + strip_2.eta*strip_2.et)/total_pt;
  return avg_eta;
}

float L1PFTauProducer::weighted_avg_phi_t_s(pftau_t tau_1, strip_t strip_1){
  float total_pt = (tau_1.p4().Pt()+strip_1.et);
  float t_phi = tau_1.p4().Phi();
  float s_phi = strip_1.phi;
  float pi = 3.14159;
  if(t_phi < 0)
    t_phi=(2*pi+t_phi);

  if(s_phi < 0)
    s_phi = (2*pi+s_phi);

  float avg_phi  = (t_phi*tau_1.p4().Pt() + s_phi*strip_1.et)/total_pt;
  if(avg_phi>pi)
    avg_phi = avg_phi - 2*pi;
  return avg_phi;
}

float L1PFTauProducer::weighted_avg_eta_t_s(pftau_t tau_1, strip_t strip_1){
  float total_pt = (tau_1.p4().Pt()+strip_1.et);
  float avg_eta  = (tau_1.p4().Eta()*tau_1.p4().Pt() + strip_1.eta*strip_1.et)/total_pt;
  return avg_eta;
}


float L1PFTauProducer::delta_r_cluster(cluster_t cluster_1, cluster_t cluster_2){
  float delta_r = 20;
  delta_r = fabs(cluster_1.p4().Eta() - cluster_2.p4().Eta()) + fabs(cluster_1.p4().Phi() - cluster_2.p4().Phi());

  return delta_r;
}

float L1PFTauProducer::delta_r_strip(strip_t strip_1, strip_t strip_2){
  float delta_r = 20;
  delta_r = fabs(strip_1.eta - strip_2.eta) + fabs(strip_1.phi - strip_2.phi);

  return delta_r;
}


/////////////
// DESTRUCTOR
L1PFTauProducer::~L1PFTauProducer()
{
}  


//////////
// END JOB
void L1PFTauProducer::endRun(const edm::Run& run, const edm::EventSetup& iSetup)
{
}

////////////
// BEGIN JOB
void L1PFTauProducer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup )
{
}


DEFINE_FWK_MODULE(L1PFTauProducer);

//  LocalWords:  PFChargedCandidates
