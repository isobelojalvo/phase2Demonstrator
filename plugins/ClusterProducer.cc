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


ClusterProducer::ClusterProducer(const edm::ParameterSet& cfg) :
  //ecalSrc_(consumes<EcalTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalDigis"))),
  debug(cfg.getUntrackedParameter<bool>("debug", false)),
  ecalTPGBToken_(consumes<EcalEBTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalTPGsBarrel"))),
  hcalSrc_(consumes<edm::SortedCollection<HcalTriggerPrimitiveDigi> >(cfg.getParameter<edm::InputTag>("hcalTPGs")))
{
  produces< L1CaloClusterCollection >( "L1Phase2CaloClusters" ).setBranchAlias("L1Phase2CaloClusters");

}

//////////
// PRODUCE

void ClusterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::unique_ptr<L1CaloClusterCollection> newClusters( new L1CaloClusterCollection );

  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

  edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgCollection;
  iEvent.getByToken(ecalTPGBToken_,ecaltpgCollection);
  vector<ecalCrystal> ecalCrystals;
  std::cout<<"getting ecal crystals"<<std::endl;
  getEcalCrystals(ecaltpgCollection, ecalCrystals);
  std::cout<<"got ecal crystals"<<std::endl;
  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcalTPGCollection;
  vector<TLorentzVector> *hcalTPGs = new std::vector<TLorentzVector>;
  std::cout<<" handle getting tpgs"<<std::endl;
  if(!iEvent.getByToken(hcalSrc_, hcalTPGCollection))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;

  edm::ESHandle<L1CaloHcalScale> hcalScale;
  iSetup.get<L1CaloHcalScaleRcd>().get(hcalScale);
  std::cout<<"getting tpgs"<<std::endl;
  getHcalTPGs(hcalTPGCollection, hcalScale, hcalTPGs);
  //take input ecal crystals and sort by eta/phi
  // ---> Patterns produed should use sorted ecal crystals

  L1CaloCluster cluster;
  newClusters->push_back(cluster);

  
  //For each ieta/iphi produce a grid of 5x5 (make configurable to 3x5)
  //Dataformat calocluster should have a member which is the internal crystals in an array
  // loop over all towers in eta
  for(int tEta = -20; tEta < 20; tEta++)
    //loop over all towers in phi
    for(int tPhi = 0; tPhi < 36; tPhi++){
      //build cluster object
      
      float crystals[5][5] = {{0}};
      // build ecal crystals
      for(int cEta = 0; cEta < 5; cEta++){
	int crystalEtaIndex = tEta * 4 + cEta;
	//find crystal Phi
	for(int cPhi = 0; cPhi < 5 ; cPhi++ ){
	  int crystalPhiIndex = tPhi * 4 + cPhi;
	}
      }
    }

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

void ClusterProducer::getEcalCrystals(edm::Handle<EcalEBTrigPrimDigiCollection> ecalTPGs, vector<ecalCrystal> ecalCrystals){
  
  for(auto& tpg : *ecalTPGs.product())
    {
      if(tpg.encodedEt() > 0) 
	{

	  GlobalVector position;
	  auto cell = ebGeometry->getGeometry(tpg.id());

	  //position = GlobalVector(cell->getPosition().x(), cell->getPosition().y(), cell->getPosition().z());

	  float et = tpg.encodedEt()/8.;

	  if(et<0.5) continue;
	  float energy = et / sin(position.theta());
	  float eta = cell->getPosition().eta();
	  float phi = cell->getPosition().phi();
	  EBDetId detID = tpg.id();
	  float iEta = detID.ieta();
	  float iPhi = detID.iphi();

	  if(et>2){
	    std::cout<<"ET "<<et<<std::endl;
	    std::cout<<"eta  "<< eta<< " phi  "<< phi<<std::endl;
	    std::cout<<"iEta"<<iEta<<" iphi "<<iPhi<<std::endl;
	  }
	  ecalCrystal tempCrystal;
	  tempCrystal.p4.SetPtEtaPhiE(et, eta, phi,et);
	  tempCrystal.iEta = iEta;
	  tempCrystal.iPhi = iPhi;
	  tempCrystal.id = detID;

	  ecalCrystals.push_back(tempCrystal);
	}
    }


}

void ClusterProducer::getHcalTPGs( edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcaltpgCollection, edm::ESHandle<L1CaloHcalScale> &hcalScale, vector<TLorentzVector> *allHcalTPGs){
  for (size_t i = 0; i < hcaltpgCollection->size(); ++i) {

    HcalTriggerPrimitiveDigi tpg = (*hcaltpgCollection)[i];
    int cal_ieta = tpg.id().ieta();
    int cal_iphi = tpg.id().iphi();
    if(cal_ieta>28)continue; 
    if(cal_ieta<-28)continue; 
    int ieta = TPGEtaRange(cal_ieta);
    short absieta = std::abs(tpg.id().ieta());
    short zside = tpg.id().zside();
    double et = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
    //if(et>0)
    //std::cout<<"HCAL ET "<<et<<std::endl;
    if(ieta<0){
      std::cout<<"sorry, ieta less than 1 :("<<std::endl;
      std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
  }
    float eta = getRecoEta(ieta, zside);
    float phi = getRecoPhi(cal_iphi);    
    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);
    allHcalTPGs->push_back(temp);
  }
  
}


DEFINE_FWK_MODULE(ClusterProducer);
