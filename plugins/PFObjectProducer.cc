// -*- C++ -*-
//
// Package:    PFObjectProducer
// Class:      PFObjectProducer
// 
/**\class PFObjectProducer PFObjectProducer.cc L1Trigger/phase2Demonstrator/plugins/PFObjectProducer.cc

 Description: Level 1 PFObjects for the Demonstrator

 Implementation:

   Caloclusters are filled using following format
   *                           2880 TOTAL CLUSTERS 
   *       iEta
   *       |-20                    -1 || 1                        20|
   *       |-------------------------------------------------------- 
   * iPhi 1|   0   72  144            ||1440                    2808
   *      2|   1   73  145            ||1441
   *      .|   .    .    .            || 
   *      .|   .    .    .            ||
   *      .|   .    .    .            ||                        2878
   *     72|  71  143  215            ||                        2879
   
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

  produces< L1CaloClusterCollection >( "L1NeutralClusters" ).setBranchAlias("L1NeutralClusters");
  produces< L1PFObjectCollection >( "L1PFObjects" ).setBranchAlias("L1PFObjects");
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
  
  // Counters for Debugging the Event
  int nElectrons = 0;
  int nPhotons = 0;
  int nChargedHadrons = 0;
  int nNeutralHadrons = 0;
 
  for(unsigned int i = 0; i < l1CaloClusters->size(); i++){
    L1CaloCluster cluster = l1CaloClusters->at(i);
    L1CaloCluster newNeutralCluster = cluster;
    newL1NeutralClusters->push_back(newNeutralCluster);
  }

  for(auto l1Track : l1Tracks){
    L1PFObject newL1PFObject;
    float trackEta = l1Track.getMomentum().eta();
    float trackPhi = l1Track.getMomentum().phi();

    if(fabs(trackEta)>1.74)continue;

    uint32_t index = findTheIndexFromReco(trackEta, trackPhi);
    if(index > newL1NeutralClusters->size()){
      std::cout<<"-----------------------------------FATAL ERROR PFOBjectProducer-----------------------------------"<<std::endl;
      std::cout<<"Error the calculated index is greater than the size of the clusters. I am skipping this track"<<std::endl;
      std::cout<<" track Eta = "<<trackEta<<" Phi = "<<trackPhi<<std::endl;
      continue;
    }

    //find the cluster that matches
    L1CaloCluster cluster =  newL1NeutralClusters->at(index);
    cluster.setEt(newL1NeutralClusters->at(index).et());

    //note: implement me
    newL1PFObject.setTrackRef(l1Track);
    newL1PFObject.setPtEtaPhiE(l1Track.getMomentum().perp(),trackEta, trackPhi,l1Track.getMomentum().perp());

    // Take the cluster, use h/e to determine if hadron or electron/pi0 
    newL1PFObject.setHcalEnergy(cluster.hcalEnergy());
    newL1PFObject.setEcalEnergy(cluster.ecalEnergy());
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
      newL1NeutralClusters->at(index).setEt(cluster.et() - l1Track.getMomentum().perp());
      
      //Fix me
      if(cluster.EoH() > input_EoH_cut_){
	newL1NeutralClusters->at(index).setIsPhoton(true);
	if(debug)
	  std::cout<<"Found Photon Deposit"<<std::endl;

	nPhotons++;
      }
      else{
	if(debug)
	  std::cout<<"Found Hadron Deposit"<<std::endl;
	newL1NeutralClusters->at(index).setIsNeutralHadron(true);
	nNeutralHadrons++;
      }
    }
    
    //else if(l1Track.getMomentum().perp() > cluster.et()){
      
      if(debug){
	//note implement the tower eta for l1Track and cluster p4 (if not already there)
	std::cout<<"---------------------------------"<<std::endl;
	//std::cout<<"The Track Momentum is greater than the Cluster Momentum"<<std::endl;
	std::cout<<"Cluster "<< l1CaloClusters->at(index) <<std::endl;
	std::cout<<"PFObject "<< newL1PFObject<<std::endl;
	std::cout<<"---------------------------------"<<std::endl;
	//}
	}
    
    newL1PFObjects->push_back(newL1PFObject);
  }


  std::cout<<"------------ Event Summary From PF Object Producer ------------"<<std::endl;
  std::cout<<"nElectrons:      "<<nElectrons<<std::endl;
  std::cout<<"nPhotons:        "<<nPhotons<<std::endl;
  std::cout<<"nChargedHadrons: "<<nChargedHadrons<<std::endl;
  std::cout<<"nNeutralHadrons: "<<nNeutralHadrons<<std::endl;



  iEvent.put( std::move(newL1PFObjects) ,       "L1PFObjects" );
  iEvent.put( std::move(newL1NeutralClusters) , "L1NeutralClusters" );

}

uint32_t PFObjectProducer::findTheIndexFromReco( float eta, float phi){
  //uint32_t PFObjectProducer::findTheIndexFromReco(float trackET, int charge, float eta, float phi){
 
  uint32_t index;
  triggerGeometryTools tool;
  uint32_t cEta = (tool.getCrystalIEta(eta));
  uint32_t cPhi = (tool.getCrystalIPhi(phi));
  std::cout<<"crystal iEta "<<cEta <<" iPhi "<<cPhi<<" index: "<<index<<std::endl;
  int iEta = tool.convertGenEta(eta);
  int iPhi = tool.convertGenPhi(phi);

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
