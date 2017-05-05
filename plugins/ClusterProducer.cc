// -*- C++ -*-
//
// Package:    ClusterProducer
// Class:      ClusterProducer
// 
/**\class ClusterProducer ClusterProducer.cc L1Trigger/ClusterProducer/plugins/ClusterProducer.cc

 Description: Level 1 Clusters for the Demonstrator

 Implementation:
     [Notes on implementation]
*/


#include "L1Trigger/phase2Demonstrator/interface/ClusterProducer.hh"


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//


ClusterProducer::ClusterProducer(const edm::ParameterSet& cfg) :
  //ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  debug(cfg.getUntrackedParameter<bool>("debug", false)),
  useECalEndcap(cfg.getParameter<bool>("useECalEndcap")),
  ecalRecHitEBToken_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ecalRecHitEB"))),
  ecalRecHitEEToken_(consumes<EcalRecHitCollection>(cfg.getParameter<edm::InputTag>("ecalRecHitEE"))),
  hcalSrc_(consumes<HcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("hcalDigis")))
{
  produces< L1CaloClusterCollection >( "L1Phase2CaloClusters" ).setBranchAlias("L1Phase2CaloClusters");

}

//////////
// PRODUCE

void ClusterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup){
  std::unique_ptr<L1CaloClusterCollection> newClusters( new L1CaloClusterCollection );


  edm::Handle<EcalRecHitCollection> pcalohits;
  iEvent.getByToken(ecalRecHitEBToken_,pcalohits);

  edm::Handle<EcalRecHitCollection> pcalohitsEndcap;
  if(useECalEndcap)
    iEvent.getByToken(ecalRecHitEEToken_,pcalohitsEndcap);

  //take input ecal crystals and sort by eta/phi
  // ---> Patterns produed should use sorted ecal crystals
  L1CaloCluster cluster;
  newClusters->push_back(cluster);

  //For each ieta/iphi produce a grid of 5x5 (make configurable to 3x5)
  //Dataformat calocluster should have a member which is the internal crystals in an array

  //for each region create the ID bits and central max (and what else?)

  //push out the collection
  iEvent.put( std::move(newClusters) , "L1Phase2CaloClusters" );
}


/////////////
// DESTRUCTOR
ClusterProducer::~ClusterProducer()
{
}  


//////////
// END JOB
void ClusterProducer::endRun(const edm::Run& run, const edm::EventSetup& iSetup)
{
}

////////////
// BEGIN JOB
void ClusterProducer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup )
{
}

DEFINE_FWK_MODULE(ClusterProducer);
