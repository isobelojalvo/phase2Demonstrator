// -*- C++ -*-
//
// Package:    PFObjectProducer
// Class:      PFObjectProducer
// 
/**\class PFObjectProducer PFObjectProducer.cc L1Trigger/phase2Demonstrator/plugins/PFObjectProducer.cc

 Description: Level 1 PFObjects for the Demonstrator

 Implementation:
     [Notes on implementation]
*/


#include "L1Trigger/phase2Demonstrator/interface/PFObjectProducer.hh"


PFObjectProducer::PFObjectProducer(const edm::ParameterSet& cfg) :
  debug(cfg.getUntrackedParameter<bool>("debug", false)),
  input_EoH_cut_(cfg.getUntrackedParameter<int>("EoH_cut", 2)), // LSB is 0.1 so 2 corresonds to 0.2
  input_HoE_cut_(cfg.getUntrackedParameter<int>("HoE_cut", 8)), // LSB is 0.1 so 8 corresonds to 0.8
  L1ClustersToken_(consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Clusters")))
{
  L1TrackInputTag = cfg.getParameter<edm::InputTag>("L1TrackInputTag");
  ttTrackToken_ = consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag);  

  produces< L1CaloClusterCollection >( "L1Phase2NeutralClusters" ).setBranchAlias("L1Phase2NeutralClusters");
  produces< L1PFObjectCollection >( "L1Phase2PFObjects" ).setBranchAlias("L1Phase2PFObjects");
}


void PFObjectProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::unique_ptr<L1CaloClusterCollection> newL1NeutralClusters( new L1CaloClusterCollection );
  std::unique_ptr<L1PFObjectCollection>    newL1PFObjects( new L1PFObjectCollection );

  using namespace edm;
  edm::Handle< std::vector<L1CaloCluster> > l1CaloClusters;
  iEvent.getByToken( L1ClustersToken_, l1CaloClusters);


  // L1 tracks  
  std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > l1Tracks;
  edm::Handle< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > l1trackHandle;
  iEvent.getByToken(ttTrackToken_, l1trackHandle);
  l1Tracks.clear();

  //Find and sort the tracks
  for(size_t track_index=0; track_index<l1trackHandle->size(); ++track_index)
    {

      edm::Ptr< TTTrack< Ref_Phase2TrackerDigi_ > > ptr(l1trackHandle, track_index);
      double pt  = ptr->getMomentum().perp();
      double eta = ptr->getMomentum().eta();

      //only using tracks with eta less than 1.5 and pt greater than 2.5 GeV
      if(abs(eta)<1.5 && pt > 2.5)
	l1Tracks.push_back(l1trackHandle->at(track_index));       

    }

  if(debug){
    std::cout<<"PFObject Producer "<<l1trackHandle->size()<<" tracks found"<<std::endl;}

  std::sort(l1Tracks.begin(), l1Tracks.end(), [](TTTrack< Ref_Phase2TrackerDigi_ > i,TTTrack< Ref_Phase2TrackerDigi_ > j){return(i.getMomentum().perp() > j.getMomentum().perp());});   

  if(debug){
    std::cout<<"PFObject Producer "<<l1CaloClusters->size()<<" caloclusters found"<<std::endl;}

  int nElectrons = 0;
  int nPhotons = 0;
  int nChargedHadrons = 0;
  int nNeutralHadrons = 0;
  
  for(unsigned int i = 0; i < l1CaloClusters->size(); i++){
    L1CaloCluster cluster = l1CaloClusters->at(i);
    L1CaloCluster newNeutralCluster;

    newNeutralCluster.setEt( cluster.et() );
    newNeutralCluster.setTowerEta( cluster.towerEta() );
    newNeutralCluster.setTowerPhi( cluster.towerPhi() );

    for(auto l1Track : l1Tracks){
      L1PFObject newL1PFObject;

      
      //change me to ieta iphi
      if((l1Track.getMomentum().eta() - cluster.p4().Eta() < 0.087) && (l1Track.getMomentum().phi() - cluster.p4().Phi() < 0.087)){

	//note: implement me
	newL1PFObject.setTrackRef(l1Track);

	// Take the cluster, use h/e to determine if hadron or electron/pi0 
	// note: implement me
	newL1PFObject.setHcalEnergy(cluster.hcalEnergy());
	newL1PFObject.setHcalEnergy(cluster.ecalEnergy());
	newL1PFObject.setCaloEnergy(cluster.caloEnergy());
	newL1PFObject.setHoE(cluster.HoE());
	newL1PFObject.setEoH(cluster.EoH());

	// Electron ID
	if(cluster.EoH() > input_EoH_cut_){
	  newL1PFObject.setIsElectron(true);
	  nElectrons++;
	}
	if(cluster.HoE() > input_HoE_cut_){
	  newL1PFObject.setIsChargedHadron(true);
	  nChargedHadrons++;
	}

	/// Subtract track pt to create Neutral Photons and Hadrons
	if(l1Track.getMomentum().perp() < cluster.et()){
	  newNeutralCluster.setEt(cluster.et() - l1Track.getMomentum().perp());

	  if(cluster.EoH() > input_EoH_cut_){
	    newNeutralCluster.setIsPhoton(true);
	    nPhotons++;
	  }
	  else{
	    newNeutralCluster.setIsNeutralHadron(true);
	    nNeutralHadrons++;
	  }
	}
	else if(l1Track.getMomentum().perp() > cluster.et()){
	  
	  if(debug){
	    //note implement the tower eta for l1Track and cluster p4 (if not already there)
	    std::cout<<"---------------------------------"<<std::endl;
	    std::cout<<"The Track Momentum is greater than the Cluster Momentum"<<std::endl;
	    std::cout<<"Cluster iEta: "<< cluster.towerEta() << " iPhi: "<< cluster.towerPhi() <<" Reco eta: "<< cluster.p4().Eta() <<" Reco Phi "<< cluster.p4().Phi()<<std::endl;
	    std::cout<<"Track Reco eta: "<< l1Track.getMomentum().eta() <<" Reco Phi "<< l1Track.getMomentum().phi()<<std::endl; //   iEta: "<< l1Track.towerEta()   << " iPhi: "<< l1Track.towerPhi()   
	      
	    std::cout<<"---------------------------------"<<std::endl;
	  }
	}
	newL1PFObjects->push_back(newL1PFObject);
      }
    }

    newL1NeutralClusters->push_back(newNeutralCluster);
  }
  std::cout<<"------------ Event Summary From PF Object Producer ------------"<<std::endl;
  std::cout<<"nElectrons:      "<<nElectrons<<std::endl;
  std::cout<<"nPhotons:        "<<nPhotons<<std::endl;
  std::cout<<"nChargedHadrons: "<<nChargedHadrons<<std::endl;
  std::cout<<"nNeutralHadrons: "<<nNeutralHadrons<<std::endl;



  iEvent.put( std::move(newL1PFObjects) ,       "L1PFObjects" );
  iEvent.put( std::move(newL1NeutralClusters) , "L1NeutralClusters" );

}



/////////////
// DESTRUCTOR
PFObjectProducer::~PFObjectProducer()
{
}  


//////////
// END JOB
void PFObjectProducer::endRun(const edm::Run& run, const edm::EventSetup& iSetup)
{
}

////////////
// BEGIN JOB
void PFObjectProducer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup )
{
}

DEFINE_FWK_MODULE(PFObjectProducer);
