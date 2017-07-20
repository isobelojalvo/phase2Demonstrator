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
  three_prong_delta_r_( cfg.getUntrackedParameter<double>("three_prong_dr", 0.1)), // LSB is 0.1 so 8 corresonds to 0.8
  isolation_delta_r_(   cfg.getUntrackedParameter<double>("three_prong_dr", 0.5)), // LSB is 0.1 so 8 corresonds to 0.8
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

  //Find taus using charged hadrons
  for(unsigned int iCand = 0; iCand < l1PFChargedCandidates->size(); iCand++){

    //Find Tau Seed
    if(l1PFChargedCandidates->at(iCand).isChargedHadron() == true && l1PFChargedCandidates->at(iCand).p4().Pt() > 6){
      std::cout<<"New Charged Hadron Candidate Search"<<std::endl;
      L1PFObject tau_cand[3];
      L1PFObject pfChargedHadron_Seed = l1PFChargedCandidates->at(iCand);
      tau_cand[0] = pfChargedHadron_Seed;
      int n_prongs_found = 1;
      float isolationSum = 0;
      L1PFObject electronGrid[5][5];
      float isolationSumChargedHadron;
      std::cout<<"Seed Cand Pt: "<<pfChargedHadron_Seed.p4().Pt()<<" Eta: "<<pfChargedHadron_Seed.p4().Eta()<<" Phi: "<<pfChargedHadron_Seed.p4().Phi()<<std::endl;
      for(unsigned int jCand = 0; jCand < l1PFChargedCandidates->size(); jCand++){

	if(jCand <= iCand)
	  continue;

	if(l1PFChargedCandidates->at(jCand).isChargedHadron()){

	  L1PFObject pfChargedHadron_signal_cand = l1PFChargedCandidates->at(jCand);

	  std::cout<<"charged Hadron Found, Delta R: "<<(fabs(pfChargedHadron_signal_cand.p4().Eta()-pfChargedHadron_Seed.p4().Eta())+fabs(pfChargedHadron_signal_cand.p4().Phi()-pfChargedHadron_Seed.p4().Phi()))<<std::endl;	  

	  if(Delta_R(pfChargedHadron_signal_cand.p4().Eta(), pfChargedHadron_signal_cand.p4().Phi(), 
		     pfChargedHadron_Seed.p4().Eta(), pfChargedHadron_Seed.p4().Phi(), 
		     three_prong_delta_r_)){
	    n_prongs_found++;

	  if(n_prongs_found==2){
	    tau_cand[1] = pfChargedHadron_signal_cand;
	  }

	  if(n_prongs_found==3){
	    tau_cand[2] = pfChargedHadron_signal_cand;
	    std::cout<<"3 charged hadrons found"<<std::endl;
	    std::cout<<"prong 1: "<<tau_cand[0].p4().Pt()<<" prong 2: "<<tau_cand[1].p4().Pt()<<" prong 3: "<<tau_cand[2].p4().Pt()<<std::endl;
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
	if(pfChargedHadron_Seed.p4().Pt() > 15){
	  if(n_prongs_found==2){
	    //Add second prong to Isolation Cone
	    isolationSumChargedHadron += tau_cand[2].p4().Pt();	    
	  }
	  L1PFTau tempL1PFTau;
	  tempL1PFTau.setPtEtaPhiE(pfChargedHadron_Seed.p4().Pt(), 
				   pfChargedHadron_Seed.p4().Eta(), 
				   pfChargedHadron_Seed.p4().Phi(), 
				   pfChargedHadron_Seed.p4().Pt());

	  tempL1PFTau.setTauType(0);
	  tempL1PFTau.setRawIso(isolationSum);
	  tempL1PFTau.setRelIso(isolationSum/pfChargedHadron_Seed.p4().Pt());

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
	std::cout<<"total PT: "<<totalPT<<" averageEta: "<<averageEta<<" averagePhi: "<<averagePhi<<std::endl;
	tempL1PFTau.setPtEtaPhiE(totalPT, 
				 averageEta,
				 averagePhi, 
				 totalPT);
	
	tempL1PFTau.setTauType(10);
	tempL1PFTau.setRawIso(isolationSum);
	tempL1PFTau.setRelIso(isolationSum/totalPT);

	//electronGrid[nTauCands] = electronGridTemp;	
	tauCands[nTauCands] = tempL1PFTau;
	nTauCands++;
      }
    }
    if(nTauCands > 12){
      break;
    }
  }
  // For each Tau Candidate Create the Strip by grabbing the index of the nearby neutrals

  //find pi0's using electrons/gamma
  //stitch the strips by combining electron and gamma candidates
  for(unsigned int i = 0; i < nTauCands; i++){
     if(tauCands[i].tauType()==10)
       continue;

     //needs tower eta and tower phi, as well as pt, eta, phi
     //strip_t tempStrip[5];    
     //strip_t finalStrip;
     //float crystal_distance_ = 2;

     /*
     // Create Grid of 5x5 neutral Clusters
     uint32_t index_cen = findTheIndexFromReco(trackEta, trackPhi,  0, 0);
     uint32_t index_p   = findTheIndexFromReco(trackEta, trackPhi,  1, 0);
     uint32_t index_pp  = findTheIndexFromReco(trackEta, trackPhi,  2, 0);
     uint32_t index_m   = findTheIndexFromReco(trackEta, trackPhi, -1, 0);
     uint32_t index_mm  = findTheIndexFromReco(trackEta, trackPhi, -2, 0);

     // Check if each iPhi is strip like

     Check_And_Merge_Strip(neutralCluster[index_mm + iPhi - 2],  neutralCluster[index_mm  + iPhi - 1], tempStrip[0], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_mm + iPhi - 1],  neutralCluster[index_mm  + iPhi - 0], tempStrip[0], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_mm + iPhi + 0],  neutralCluster[index_mm  + iPhi + 1], tempStrip[0], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_mm + iPhi + 1],  neutralCluster[index_mm  + iPhi + 2], tempStrip[0], crystal_distance_);

     Check_And_Merge_Strip(neutralCluster[index_m   + iPhi - 2], neutralCluster[index_m   + iPhi - 1], tempStrip[1], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_m   + iPhi - 1], neutralCluster[index_m   + iPhi - 0], tempStrip[1], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_m   + iPhi + 0], neutralCluster[index_m   + iPhi + 1], tempStrip[1], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_m   + iPhi + 1], neutralCluster[index_m   + iPhi + 2], tempStrip[1], crystal_distance_);

     Check_And_Merge_Strip(neutralCluster[index_cen + iPhi - 2], neutralCluster[index_cen + iPhi - 1], tempStrip[2], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_cen + iPhi - 1], neutralCluster[index_cen + iPhi - 0], tempStrip[2], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_cen + iPhi + 0], neutralCluster[index_cen + iPhi + 1], tempStrip[2], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_cen + iPhi + 1], neutralCluster[index_cen + iPhi + 2], tempStrip[2], crystal_distance_);

     Check_And_Merge_Strip(neutralCluster[index_p   + iPhi - 2], neutralCluster[index_p   + iPhi - 1], tempStrip[3], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_p   + iPhi - 1], neutralCluster[index_p   + iPhi - 0], tempStrip[3], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_p   + iPhi + 0], neutralCluster[index_p   + iPhi + 1], tempStrip[3], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_p   + iPhi + 1], neutralCluster[index_p   + iPhi + 2], tempStrip[3], crystal_distance_);

     Check_And_Merge_Strip(neutralCluster[index_pp  + iPhi - 2], neutralCluster[index_pp  + iPhi - 1], tempStrip[4], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_pp  + iPhi - 1], neutralCluster[index_pp  + iPhi - 0], tempStrip[4], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_pp  + iPhi + 0], neutralCluster[index_pp  + iPhi + 1], tempStrip[4], crystal_distance_);
     Check_And_Merge_Strip(neutralCluster[index_pp  + iPhi + 1], neutralCluster[index_pp  + iPhi + 2], tempStrip[4], crystal_distance_);

     for(unsigned int j = 0; j < 4; j++){
       if(Position_Diff(tempStrip[j],tempStrip[j+1]) < crystal_distance){
	 if(tempStrip[j].et() + tempStrip[j+1] > finalStrip.et()){
	   float et =  tempStrip[j].et() + tempStrip[j+1].et();
	   finalStrip.setEt(  et );
	   finalStrip.setEta( eta );
	   finalStrip.setPhi( (tempStrip[j].et()*tempStrip[j].Phi() + tempStrip[j+1].et()*tempStrip[j+1].Phi())/ et );
	 }
       }
     }

   //find 1 prong pi0's by combining the 1 prong and the pi0's
     if(tauCands[i].tauType()==0 && finalStrip.et() > 0){
       float et = tauCands[i].p4().Pt() + finalStrip.et();
       tauCands[i].setPtEtaPhiE( et, 
				(tauCands[i].p4().Pt()*tauCands[i].p4().Eta() + finalStrip.et()*finalStrip.eta())/et, 
				(tauCands[i].p4().Pt()*tauCands[i].p4().Phi() + finalStrip.et()*finalStrip.phi())/et, 
				et);
       tempL1PFTau.setTauType(1);
     }
  */
  }
  //now fill the new collection
   for(unsigned int i = 0; i < 12; i++){
     newL1PFTauCollection->push_back(tauCands[i]);
     //if(debug){
       std::cout<<"------ SUMMARY OF THE TAU CANDS ------"<<std::endl;
       if(tauCands[i].p4().Et()>0){
	std::cout<<tauCands[i]<<std::endl;
	//}
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
    std::cout<<"DeltaR Returning true"<<std::endl;
    return true;
  }
  else{
    std::cout<<"DeltaR Returning false"<<std::endl;
    return false;
  }
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
