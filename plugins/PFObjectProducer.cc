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

  std::sort(l1Tracks.begin(), l1Tracks.end(), [](TTTrack< Ref_Phase2TrackerDigi_ > i,TTTrack< Ref_Phase2TrackerDigi_ > j){return(i.getMomentum().perp() > j.getMomentum().perp());});   
  
  
  for(unsigned int i = 0; i < l1CaloClusters->size(); i++){
    L1CaloCluster cluster = l1CaloClusters->at(i);
    L1CaloCluster newNeutralCluster;
    L1PFObject newL1PFObject;
    newNeutralCluster.setEt(cluster.et());
    newNeutralCluster.setTowerEta(cluster.towerEta());
    newNeutralCluster.setTowerPhi(cluster.towerPhi());

    for(auto l1Track : l1Tracks){
      if((l1Track.getMomentum().eta() - cluster.p4().Eta() < 0.087) && (l1Track.getMomentum().phi() - cluster.p4().Phi() < 0.087)){

	/// Take the cluster, use h/e to determine if hadron or electron/pi0 

	/// Subtract track pt to create Neutrals
	newNeutralCluster.setEt(cluster.et()-l1Track.getMomentum().perp());
      }
    }

    newL1PFObjects->push_back(newL1PFObject);
    newL1NeutralClusters->push_back(newNeutralCluster);
  }

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
