// -*- C++ -*-
//
// Package:    ClusterProducer
// Class:      ClusterProducer
// 
/**\class ClusterProducer ClusterProducer.cc L1Trigger/ClusterProducer/plugins/ClusterProducer.cc

 Description: Level 1 Clusters for the Demonstrator

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


#include "L1Trigger/phase2Demonstrator/interface/ClusterProducer.hh"
#define LSB 10

ClusterProducer::ClusterProducer(const edm::ParameterSet& cfg) :
  debug(cfg.getUntrackedParameter<bool>("debug", false)),
  ecalTPGBToken_(consumes<EcalEBTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalTPGsBarrel"))),
  hcalSrc_(consumes<edm::SortedCollection<HcalTriggerPrimitiveDigi> >(cfg.getParameter<edm::InputTag>("hcalTPGs")))
{
  produces< L1CaloClusterCollection >( "L1Phase2CaloClusters" ).setBranchAlias("L1Phase2CaloClusters");

}

void ClusterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Output Collecction
  std::unique_ptr<L1CaloClusterCollection> newClusters( new L1CaloClusterCollection );

  // Set up the ECAL Crystals
  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

  edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgCollection;
  iEvent.getByToken( ecalTPGBToken_, ecaltpgCollection);
  vector<ecalCrystal_t> ecalCrystals;


  // Get the ECAL Crystals as ecalCrystal_t
  getEcalCrystals(ecaltpgCollection, ecalCrystals);
  

  edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcalTPGCollection;
  vector<TLorentzVector> hcalTPGs;
  
  if(!iEvent.getByToken(hcalSrc_, hcalTPGCollection))
    std::cout<<"ERROR GETTING THE HCAL TPGS"<<std::endl;

  edm::ESHandle<L1CaloHcalScale> hcalScale;
  iSetup.get<L1CaloHcalScaleRcd>().get(hcalScale);


  // Get the HCAL TPGs as TLorentzVector 
  getHcalTPGs(hcalTPGCollection, hcalScale, hcalTPGs);

  int nPrints = 0;
  // For each ieta/iphi produce a grid of 5x5 (make configurable to 3x5?)
  // The collection is produced in the ordered format as seen above. 
  // Loop over all towers in eta
  for(int tEta = -20; tEta < 21; tEta++){

    //skip tower eta == 0
    if(tEta == 0)
      continue;

    //loop over all towers in phi
    for(int tPhi = 0; tPhi < 72; tPhi++){

      //build cluster object
      //tRecoEta is from -1.74 to 1.74 add 0.0435 to get the center of the bin
      float tRecoEta = tEta * 0.087 + 0.087/2;      
      //tRecoEta is from 0 to 2pi (bin is 2pi/72) add 0.0436 to get the center of the bin
      float tRecoPhi = tPhi * 0.0871380 + 0.0871380/2;
      
      L1CaloCluster tempCluster;
      
      // Intialize member data
      float HCALEt = 0;
      float pt     = 0;
      float central_eta  = tRecoEta; //central cluster Reco Eta
      float central_phi  = tRecoPhi; //central cluster Reco Phi
      float sumCrystals  = 0; 
      float maxCrystalEt = 0;
      float crystals[5][5] = {{0}};

      unsigned int maxCrystalEta = 2; 
      unsigned int maxCrystalPhi = 2;
      
      //Find matching HcalTPG
      for(auto hcalTPG : hcalTPGs){
	//fixme
	if( abs(hcalTPG.Eta() - tRecoEta)+abs(hcalTPG.Phi() - tRecoPhi)< 0.08727 ){
	  HCALEt = hcalTPG.Pt(); 
	  break;
	}
      }

      // find matching ecal crystals
      for(unsigned int cEta = 0; cEta < 5; cEta++){
	//find crystal Phi
	for(unsigned int cPhi = 0; cPhi < 5 ; cPhi++ ){
	  ecalCrystal_t foundCrystal;

	  crystals[cEta][cPhi] =  findEcalCrystal(tEta, tPhi, cEta, cPhi, ecalCrystals, foundCrystal);

	  sumCrystals += crystals[cEta][cPhi];

	  if(crystals[cEta][cPhi] > maxCrystalEt){
	    maxCrystalEt = crystals[cEta][cPhi];
	    maxCrystalEta = cEta;
	    maxCrystalPhi = cPhi;
	    central_eta = foundCrystal.p4.Eta();
	    central_phi = foundCrystal.p4.Phi();
	  }
	  //DEBUG statement
	  //if(crystals[cEta][cPhi]>0.5){
	    //float outputPhi =  (tPhi*5 + cPhi)*0.0174276;
	    //if(outputPhi > 3.14159) outputPhi = outputPhi - 3.14159*2;
	    //std::cout<<"recoEta "<<(tEta*5 + cEta)*0.0174<<" recoPhi "<< outputPhi<<std::endl;
	    //std::cout<<"tEta "<<tEta<<" tPhi "<<tPhi<<std::endl;
	    //std::cout<<"crystal et "<< crystals[cEta][cPhi]<< " cEta, cPhi "<< cEta<< ", " << cPhi<<std::endl;
	    //}
	}
      } // End ecal crystal matching

      pt = sumCrystals + HCALEt;
      clusterAlgoMax(crystals, maxCrystalEta, maxCrystalPhi);
      //clusterAlgoMaxStrips(crystals, maxCrystalEta, maxCrystalPhi);
      //clusterAlgoMaxInternet(crystals, maxCrystalEta, maxCrystalPhi);

      unsigned HoE = 0;
      unsigned EoH = 0;
      if(sumCrystals > 0 && HCALEt == 0){
	EoH = (unsigned) (10*1); // NOTE WE MULTIPLY BY EoH LSB HERE
      }
      else{
	EoH = (unsigned) (10*sumCrystals/HCALEt); // NOTE WE MULTIPLY BY EoH LSB HERE
      }

      if(sumCrystals == 0 && HCALEt > 0){
	HoE = (unsigned) (10*1);// NOTE WE MULTIPLY BY HoE LSB HERE
      }
      else{
	HoE = (unsigned) (10*HCALEt/sumCrystals); // NOTE WE MULTIPLY BY HoE LSB HERE
      }

      unsigned etaSide = 0;
      if(tEta>0)
	etaSide = 1;

      // RECO p4 only for debugging
      tempCluster.setPtEtaPhiE(pt, central_eta, central_phi, pt);
      
      tempCluster.setEt(            (unsigned) (LSB * pt)      ); // NOTE WE MULTIPLY BY LSB HERE
      tempCluster.setEcalEnergy(    sumCrystals                );
      tempCluster.setHcalEnergy(    HCALEt                     );

      tempCluster.setTowerPhi(      (unsigned) (tPhi*5+maxCrystalPhi)     );
      tempCluster.setTowerEta(      (unsigned) abs(tEta*5+maxCrystalEta)  );
      tempCluster.setTowerEtaSide(  (unsigned) etaSide                    );
      tempCluster.setEoH( EoH );
      tempCluster.setHoE( HoE );

      if(debug)
	if(pt>2){
	  std::cout<<tempCluster<<std::endl;
	  for(unsigned int cPhi = 0; cPhi < 5 ; cPhi++ ){
	    for(unsigned int cEta = 0; cEta < 5; cEta++){
	      std::cout<<" "<<std::setw(8)<< std::setprecision(5)<<crystals[cEta][cPhi];
	    }
	    std::cout<<std::endl;
	  }
	}
      if(debug)
	if(nPrints < 3 && tempCluster.et() == 0){
	  std::cout<<"-------------- printing a 0 cluster --------------"<<std::endl;
	  std::cout<<tempCluster<<std::endl;
	  nPrints++;
	}
      //DEBUG statement
      //std::cout<<"eta "<< central_eta <<" tempCluster.p4().Eta() "<<tempCluster.p4().Eta()<<std::endl;

      newClusters->push_back(tempCluster);
    }
  }

  //push out the collection
  iEvent.put( std::move(newClusters) , "L1Phase2CaloClusters" );
}


ClusterProducer::~ClusterProducer(){}  

void ClusterProducer::endRun(const edm::Run& run, const edm::EventSetup& iSetup){}

void ClusterProducer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup ){}

void ClusterProducer::getEcalCrystals(edm::Handle<EcalEBTrigPrimDigiCollection> ecalTPGs, vector<ecalCrystal_t> &ecalCrystals)
{
  if(debug)
    std::cout<<"-------------- ECAL Crystals --------------"<<std::endl;
    
  for(auto& tpg : *ecalTPGs.product())
    {

      if(tpg.encodedEt() > 0) 
	{

	  GlobalVector position;
	  auto cell = ebGeometry->getGeometry(tpg.id());

	  float et = tpg.encodedEt()/8.;

	  if(et<0.1) continue;// LSB is 0.1
	  //float energy = et / sin(position.theta());
	  float eta = cell->getPosition().eta();
	  float phi = cell->getPosition().phi();
	  EBDetId detID = tpg.id();
	  float iEta = detID.ieta();
	  float iPhi = detID.iphi();

	  //DEBUG STATMENT
	  if(debug)
	    if(et>1){
	      std::cout<<"ET "<<et<<std::endl;
	      std::cout<<" eta  "<< eta<< " phi  "<< phi<<std::endl;
	      std::cout<<" iEta"<<iEta<<" iphi "<<iPhi<<std::endl;
	    }
	  ecalCrystal_t tempCrystal;
	  tempCrystal.p4.SetPtEtaPhiE(et, eta, phi,et);
	  tempCrystal.iEta = iEta;
	  tempCrystal.iPhi = iPhi;
	  tempCrystal.id = detID;

	  ecalCrystals.push_back(tempCrystal);
	}
    }
}

/*
  eta  0.871948 phi  0.672644
  iEta50 iphi 49
 */
float ClusterProducer::findEcalCrystal(int tEta, int tPhi, int cEta, int cPhi, 
				       vector<ecalCrystal_t> ecalCrystals, 
				       ecalCrystal_t &foundCrystal){
  float crystalET = 0;
  for(auto ecalCrystal : ecalCrystals){

    int findEcalIEta = tEta * 5 + cEta;
    int findEcalIPhi = tPhi * 5 + cPhi + 10; // 10 crystal offset from what is provided to match to hcal

    if(ecalCrystal.iEta == findEcalIEta && 
       ecalCrystal.iPhi == findEcalIPhi){
      crystalET    = ecalCrystal.p4.Pt();
      foundCrystal = ecalCrystal;
    }

  }

  return crystalET;
}

void ClusterProducer::getHcalTPGs( edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcaltpgCollection, 
				   edm::ESHandle<L1CaloHcalScale> &hcalScale, 
				   vector<TLorentzVector> &allHcalTPGs){

  if(debug)
    std::cout<<"-------------- HCAL TPGs --------------"<<std::endl;

  for (size_t i = 0; i < hcaltpgCollection->size(); ++i) {

    HcalTriggerPrimitiveDigi tpg = (*hcaltpgCollection)[i];
    int cal_ieta = tpg.id().ieta();
    int cal_iphi = tpg.id().iphi();
    if(cal_ieta> 28)continue; 
    if(cal_ieta<-28)continue; 

    int ieta      = TPGEtaRange(cal_ieta);
    short absieta = std::abs(tpg.id().ieta());
    short zside   = tpg.id().zside();
    double et     = hcalScale->et(tpg.SOI_compressedEt(), absieta, zside); 
    //DEBUG STATEMENT
    if(debug)
      if(et>0)
	std::cout<<"HCAL ET "<<et<<std::endl;

    if(ieta<0){
      std::cout<<"sorry, ieta less than 1 :("<<std::endl;
      std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
    }

    triggerGeometryTools tool;
    float eta = tool.getRecoEta(ieta, zside);
    float phi = tool.getRecoPhi(cal_iphi);    

    TLorentzVector temp ;
    temp.SetPtEtaPhiE(et,eta,phi,et);

    allHcalTPGs.push_back(temp);
  }
  
}

void ClusterProducer::clusterAlgoMax(float crystals[5][5], unsigned int &maxCrystalEta, unsigned int &maxCrystalPhi){
  unsigned int maxCrystalEt = 0;
  maxCrystalEta = 2;
  maxCrystalPhi = 2;
  for(unsigned int cEta = 0; cEta<5; cEta++){
    for(unsigned int cPhi = 0; cPhi<5; cPhi++){
      if(crystals[cEta][cPhi]>maxCrystalEt){
	maxCrystalEt = crystals[cEta][cPhi];
	maxCrystalEta = cEta;
	maxCrystalPhi = cPhi;
      }
      else
	continue;
    }
  }
}

void ClusterProducer::clusterAlgoMaxInternet(float crystals[25], unsigned int &maxCrystalEta, unsigned int &maxCrystalPhi){

  unsigned int maxCrystalEt = 0;
  unsigned int maxCrystalIndex = 0;  // 2*5 + 2 

  unsigned int big = 0;
  unsigned int bigIndex = 0;

  maxCrystalEta = 2;
  maxCrystalPhi = 2;

  unsigned int start = 0;
  unsigned int end = 25;

  unsigned int index = start + 1;
  unsigned int n = end - start + 1;//n: the number of elements to be sorted, assuming n>0

  maxCrystalEt = crystals[start];

  for (unsigned int i = index; i < n-1; i = i+2 ){
    if ( crystals[i] < crystals[i+1] ){ //one comparison
      big = crystals[i+1];
      bigIndex = i+1;
    }
    else{
      big = crystals[i];
      bigIndex = i;
    }
    if ( maxCrystalEt < big ){ //one comparison
      maxCrystalEt = big;
      maxCrystalIndex = bigIndex;
    }
  }

  maxCrystalEta = maxCrystalIndex/5;
  maxCrystalPhi = maxCrystalIndex%5;
}


DEFINE_FWK_MODULE(ClusterProducer);
