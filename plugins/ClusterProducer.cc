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
#define activityFractionPi0 0.035;

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
  vector<hcalTPG_t> hcalTPGs;
  
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
  for(int tEta = -20; tEta < 20; tEta++){

    //skip tower eta == 0
    //if(tEta == 0)
      //continue;

    //loop over all towers in phi
    //check me, was 73, putting 72
    for(int tPhi = 0; tPhi < 72; tPhi++){

      //build cluster object
      //tRecoEta is from -1.74 to 1.74 add 0.0435 to get the center of the bin
      float tRecoEta = tEta * 0.087;
      //if(tEta<0)
      //tRecoEta += 0.087/2;      
      //else
      //tRecoEta -= 0.087/2;      

      //tRecoEta is from 0 to 2pi (bin is 2pi/72) add 0.0436 to get the center of the bin
      float tRecoPhi = tPhi * 0.0871380;// + 0.0871380/2;
      if(tRecoPhi > 3.14159)
	tRecoPhi = tRecoPhi - 2*3.14159;

      L1CaloCluster tempCluster;
      
      // Intialize member data
      float HCALEt = 0;
      float pt     = 0;
      float central_eta  = tRecoEta + 0.087/2; //central cluster Reco Eta
      float central_phi  = tRecoPhi + 0.087/2; //central cluster Reco Phi
      float ecal_eta = central_eta;
      float ecal_phi = central_phi;
      float sumCrystals  = 0; 
      float crystals[5][5] = {{0}};
      float foundCrystalsiEta[5][5] = {{0}};//remove me
      float foundCrystalsiPhi[5][5] = {{0}};//remove me
      float foundCrystalsEta[5][5] = {{0}};//remove me
      float foundCrystalsPhi[5][5] = {{0}};//remove me

      unsigned int maxCrystalEta = 2; 
      unsigned int maxCrystalPhi = 2;
      
      //Find matching HcalTPG
      for(auto hcalTPG : hcalTPGs){
	//fixme change to unigned
	//Figure out where to put this tower eta offset
	if((((tEta == hcalTPG.iEta - 1)&&(tEta>0))||((tEta == hcalTPG.iEta)&&(tEta<0)))&& tPhi == hcalTPG.iPhi){
	  //DEBUG fix me
	  //if(tRecoPhi > -2 && tRecoPhi < -1.7 && tRecoEta > 0.42 && tRecoEta < 0.65)
	  //if(tEta + 1 == 6 && tPhi == 51)
	  //HCALEt = 10;

	  HCALEt = hcalTPG.p4.Pt(); 
	  if(debug && hcalTPG.p4.Pt()>5){
	    std::cout<<"Reco Eta: "<<tRecoEta<<" Phi: "<<tRecoPhi<<std::endl;
	    std::cout<<"Matched HCAL Eta: "<<hcalTPG.p4.Eta()<<" Phi: "<<hcalTPG.p4.Phi()<<std::endl;
	    std::cout<<"finding tower eta: "<<tEta<<" tower phi: "<<tPhi<<std::endl;
	  }
	  
	  break;
	}
      }

      float maxCrystalEt = 0;
      // find matching ecal crystals
      for(int cEta = 0; cEta < 5; cEta++){
	//find crystal Phi
	for(int cPhi = 0; cPhi < 5 ; cPhi++ ){
	  ecalCrystal_t foundCrystal;

	  //if(tRecoPhi > -2 && tRecoPhi < -1.7 && tRecoEta > 0.42 && tRecoEta < 0.65)
	  crystals[cEta][cPhi] =  findEcalCrystal(tEta, tPhi, cEta, cPhi, ecalCrystals, foundCrystal);
	  //remove me
	  float tCrysRecoEta = tRecoEta + cEta*0.0174 + 0.0174/2;
	  float tCrysRecoPhi = tRecoPhi + cPhi*0.01745 + 0.01745;
	  foundCrystalsiEta[cEta][cPhi] = foundCrystal.id.ieta();
	  foundCrystalsiPhi[cEta][cPhi] = foundCrystal.id.iphi();
	  foundCrystalsEta[cEta][cPhi] = foundCrystal.p4.Eta();
	  foundCrystalsPhi[cEta][cPhi] = foundCrystal.p4.Phi();

	  //if(jEta > 0.4 && jEta < 0.44 && jPhi > -1.74 && jPhi <  -1.7 )
	  /*
	  if(tCrysRecoEta > 0.4 && tCrysRecoEta < 0.52)
	    if(tCrysRecoPhi > -1.9 && tCrysRecoPhi <  -1.7 ){
	      crystals[cEta][cPhi] = 3;
	      std::cout<<"setting crystal to 3"<<std::endl;
	    }	  
	  */
	  sumCrystals += crystals[cEta][cPhi];

	  if(crystals[cEta][cPhi] > maxCrystalEt){
	    maxCrystalEt = crystals[cEta][cPhi];
	    maxCrystalEta = cEta;
	    maxCrystalPhi = cPhi;
	    //central_eta = tCrysRecoEta;
	    //central_phi = tCrysRecoPhi;
	    //debug put me back
	    ecal_eta = foundCrystal.p4.Eta();
	    ecal_phi = foundCrystal.p4.Phi();
	  }
	}
      } // End ecal crystal matching

      if(HCALEt < sumCrystals){
	central_eta = ecal_eta;
	central_phi = ecal_phi;
	
      }

      bitset<5> activeTowerEtaPattern = 0;
      bitset<5> activeTowerPhiPattern = 0;
      float activeTowerLevel = sumCrystals * activityFractionPi0;

      // find matching ecal crystals
      for(int cEta = 0; cEta < 5; cEta++){
	//find crystal Phi
	bool activeStrip = false;
	for(int cPhi = 0; cPhi < 5 ; cPhi++ ){
	  if(crystals[cEta][cPhi] > activeTowerLevel) activeStrip = true;
	}
	if(activeStrip == true) activeTowerEtaPattern |= (0x1 << cEta);
      }

      // find matching ecal crystals
      for(int cPhi = 0; cPhi < 5 ; cPhi++ ){
	//find crystal Phi
	bool activeStrip = false;
	for(int cEta = 0; cEta < 5; cEta++){
	  if(crystals[cEta][cPhi] > activeTowerLevel) activeStrip = true;
	}
	if(activeStrip == true) activeTowerPhiPattern |= (0x1 << cPhi);
      }

      bool pi0Like = false;
      if(sumCrystals>2)
	pi0Like = pi0BitSet(activeTowerEtaPattern, activeTowerPhiPattern);
      if(pi0Like && debug){
	std::cout<<"Pi0 like cluster!"<<std::endl;
      }
      
      pt = sumCrystals + HCALEt;
      clusterAlgoMax(crystals, maxCrystalEta, maxCrystalPhi);

      // If a large majority of the energy is in HCAL then leave the ieta/iphi at the cenral value
      // FIX ME PERHAPS change this to getHitTowerLocation using averages
      //if(HCALEt > sumCrystals*3){
      //maxCrystalEta = 0;
      //maxCrystalPhi = 0;
      //central_eta = tRecoEta;
      //central_phi = tRecoPhi;
      //}
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

      // to be tuned
      //if(pi0Like == true && HoE > 3)
      //pi0Like = false;

      unsigned etaSide = 0;
      if(tEta>0)
	etaSide = 1;

      // RECO p4 only for debugging
      tempCluster.setPtEtaPhiE(pt, central_eta, central_phi, pt);
      
      tempCluster.setEt(            (unsigned) (LSB * pt)      ); // NOTE WE MULTIPLY BY LSB HERE
      tempCluster.setEcalEnergy(    sumCrystals                );
      tempCluster.setHcalEnergy(    HCALEt                     );

      tempCluster.setTowerPhi(      (unsigned) tPhi       );
      tempCluster.setTowerEta(      (unsigned) abs(tEta)  );
      tempCluster.setTowerEtaSide(  (unsigned) etaSide    );

      tempCluster.setCrystalPhi(      (unsigned) (tPhi*5+maxCrystalPhi)     );
      tempCluster.setCrystalEta(      (unsigned) abs(tEta*5+maxCrystalEta)  );

      tempCluster.setEoH( EoH );
      tempCluster.setHoE( HoE );
      tempCluster.setIsPi0( pi0Like );


      if(debug){
	if(pt>0){
	  std::cout<<tempCluster<<std::endl;
	  for(int cPhi = 0; cPhi < 5 ; cPhi++ ){
	    for(int cEta = 0; cEta < 5; cEta++){
	      std::cout<<" "<<std::setw(8)<< std::setprecision(5)<<crystals[cEta][cPhi];
	    }
	    std::cout<<std::endl;
	  }
	  std::cout<<"tEta "<<tEta<<" tPhi "<<tPhi<<std::endl;
	}
      }

      ///delete these
      /*
      if((tEta == 4 && tPhi == 51) || (tEta == 5 && tPhi == 51) || (tEta == 6 && tPhi == 51)){
	std::cout<<"print ieta value"<<std::endl;
	  for(int cPhi = 0; cPhi < 5 ; cPhi++ ){
	    for(int cEta = 0; cEta < 5; cEta++){
	      std::cout<<" "<<std::setw(10)<<foundCrystalsiEta[cEta][cPhi];
	    }
	    std::cout<<std::endl;
	  }
      }

      if((tEta == 4 && tPhi == 51) || (tEta == 5 && tPhi == 51) || (tEta == 6 && tPhi == 51)){
	std::cout<<"print reco eta value"<<std::endl;
	  for(int cPhi = 0; cPhi < 5 ; cPhi++ ){
	    for(int cEta = 0; cEta < 5; cEta++){
	      std::cout<<" "<<std::setw(10)<<foundCrystalsEta[cEta][cPhi];
	    }
	    std::cout<<std::endl;
	  }
      }

      if((tEta == 4 && tPhi == 51) || (tEta == 5 && tPhi == 51) || (tEta == 6 && tPhi == 51)){
	std::cout<<"print iphi value"<<std::endl;
	  for(int cPhi = 0; cPhi < 5 ; cPhi++ ){
	    for(int cEta = 0; cEta < 5; cEta++){
	      std::cout<<" "<<std::setw(10)<<foundCrystalsiPhi[cEta][cPhi];
	    }
	    std::cout<<std::endl;
	  }
      }

      if((tEta == 4 && tPhi == 51) || (tEta == 5 && tPhi == 51) || (tEta == 6 && tPhi == 51)){
	std::cout<<"print iphi value"<<std::endl;
	  for(int cPhi = 0; cPhi < 5 ; cPhi++ ){
	    for(int cEta = 0; cEta < 5; cEta++){
	      std::cout<<" "<<std::setw(10)<<foundCrystalsPhi[cEta][cPhi];
	    }
	    std::cout<<std::endl;
	  }
      }
      */
	  /// delete these

      newClusters->push_back(tempCluster);
    }
  }
  //now prune the new clusters

  //loop through all clusters
  //for(auto centralCluster:newClusters){
  for(int i = 0; i < newClusters->size(); i++ ){
  //get upper left neighbor
    int indexNW = getNeighbor(newClusters->at(i), -1, 1);
    int indexW  = getNeighbor(newClusters->at(i), -1, 0);
    int indexSW = getNeighbor(newClusters->at(i), -1,-1);
    int indexS  = getNeighbor(newClusters->at(i),  0,-1);
    int indexSE = getNeighbor(newClusters->at(i),  1,-1);
    int indexE  = getNeighbor(newClusters->at(i),  1, 0);
    int indexNE = getNeighbor(newClusters->at(i),  1, 1);
    int indexN  = getNeighbor(newClusters->at(i),  0, 1);

    if(indexNW >= -20 && indexNW < 20)
      checkAndMergeCluster(newClusters->at(i), newClusters->at(indexNW));

    if(indexW >= -20  && indexW  < 20)
      checkAndMergeCluster(newClusters->at(i), newClusters->at(indexW) );

    if(indexSW >= -20 && indexSW < 20)
      checkAndMergeCluster(newClusters->at(i), newClusters->at(indexSW));

    if(indexS >= -20  && indexS < 20)
      checkAndMergeCluster(newClusters->at(i), newClusters->at(indexS) );

    if(indexSE >= -20 && indexSE < 20)
      checkAndMergeCluster(newClusters->at(i), newClusters->at(indexSE));

    if(indexE >= -20  && indexE  < 20)
      checkAndMergeCluster(newClusters->at(i), newClusters->at(indexE) );

    if(indexNE >= -20 && indexNE < 20)
      checkAndMergeCluster(newClusters->at(i), newClusters->at(indexNE));

    if(indexN >= -20  && indexN  < 20)    
      checkAndMergeCluster(newClusters->at(i), newClusters->at(indexN) );

  }

  if(debug)
    std::cout<<"NClusters = "<<newClusters->size()<<std::endl;

  //push out the collection
  iEvent.put( std::move(newClusters) , "L1Phase2CaloClusters" );
}


ClusterProducer::~ClusterProducer(){}  

void ClusterProducer::endRun(const edm::Run& run, const edm::EventSetup& iSetup){}

void ClusterProducer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup ){}

///fix me to make generic
int ClusterProducer::getNeighbor(L1CaloCluster centralCluster, int dEta, int dPhi){
  int maxEta = 20;
  int minEta = -20;

  int towerEta = centralCluster.towerEta();
  int towerPhi = centralCluster.towerPhi();
  int index = -99;

  if( (towerEta + dEta) > minEta && (towerEta + maxEta) < maxEta ){
    towerEta = towerEta + dEta;
  }
  else 
    return -99;

  if(towerPhi > -1 && towerPhi < 72){
    towerPhi = towerPhi + dPhi;
  }

  if(towerPhi == -1){
    towerPhi = 71;
  }

  if(towerPhi == 72){
    towerPhi = 0;
  }

  triggerGeometryTools trigTools;

  return ( trigTools.getIndex(towerEta,towerPhi));
}

bool ClusterProducer::checkAndMergeCluster(L1CaloCluster &centralCluster, L1CaloCluster &neigborCluster){
  if(abs(centralCluster.crystalEta() - neigborCluster.crystalEta()) < 2 && 
     abs(centralCluster.crystalPhi() - neigborCluster.crystalPhi()) < 2 && 
       centralCluster.p4().Pt() > neigborCluster.p4().Pt() ){
      //if passes the criteria then merge the cluster
      int cluster1Et = centralCluster.p4().Pt();
      int cluster2Et = neigborCluster.p4().Pt();

      centralCluster.setEt((float)(cluster1Et+cluster1Et));
      neigborCluster.setEt((float)0);


      //std::cout<<"Merging Cluster, old Cluster ET = "<<cluster1Et<<" new Cluster ET = "<<centralCluster.p4().Pt()<<std::endl;
      return true;
    }
    return false;

}


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

	  if(et<0.001) continue;//
	  //float energy = et / sin(position.theta());
	  float eta = cell->getPosition().eta();
	  float phi = cell->getPosition().phi();
	  EBDetId detID = tpg.id();
	  //float iEta = detID.ieta();
	  float iEta = detID.ieta() - 1;
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
  int findEcalIEta = tEta * 5 + cEta;
  int findEcalIPhi = tPhi * 5 + cPhi + 11; // 10 crystal offset from what is provided to match to hcal

  for(auto ecalCrystal : ecalCrystals){

    if(ecalCrystal.iEta == findEcalIEta && 
       ecalCrystal.iPhi == findEcalIPhi){
      crystalET    = ecalCrystal.p4.Pt();
      foundCrystal = ecalCrystal;
      
      if(debug&&crystalET>5){
	std::cout<<"ecalCrystal pt: "<<ecalCrystal.p4.Pt()<<" eta: "<<ecalCrystal.p4.Eta()<<" phi: "<<ecalCrystal.p4.Phi()<< " ieta "<<ecalCrystal.iEta<< " iphi "<<findEcalIPhi <<std::endl;
      }
      break;
    }

  }

  return crystalET;
}

void ClusterProducer::getHcalTPGs( edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcaltpgCollection, 
				   edm::ESHandle<L1CaloHcalScale> &hcalScale, 
				   vector<hcalTPG_t> &allHcalTPGs){

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
    if(ieta<0){
      std::cout<<"sorry, ieta less than 1 :("<<std::endl;
      std::cout<<"cal_ieta "<<cal_ieta<<" ieta "<<ieta<<std::endl;
    }

    triggerGeometryTools tool;
    float eta = tool.getRecoEta(ieta, zside);
    float phi = tool.getRecoPhi(cal_iphi);    

    //DEBUG STATEMENT
    if(debug)
      if(et>0){
	std::cout<<"HCAL ET "<<et<< " eta: "<<eta <<" phi: "<<phi <<std::endl;
	std::cout<<"ieta: "<<cal_ieta<<" iphi: "<<cal_iphi<<std::endl;
      }
    
    hcalTPG_t temp;
    temp.p4.SetPtEtaPhiE(et,eta,phi,et);
    temp.iEta = cal_ieta;
    temp.iPhi = cal_iphi - 1;
    temp.id = tpg.id();

    /*
    if(cal_iphi>12)
      temp.iPhi = cal_iphi-12;
    else
      temp.iPhi = 72-cal_iphi;
    */
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

bool ClusterProducer::pi0BitSet(bitset<5> etaPattern, bitset<5> phiPattern){
  // 1 + 1's patterns
  bitset<5> badPattern5( std::string("00101"));
  bitset<5> badPattern10(std::string("01010"));
  bitset<5> badPattern20(std::string("10100"));
  bitset<5> badPattern9( std::string("01001"));
  bitset<5> badPattern18(std::string("10010"));
  bitset<5> badPattern17(std::string("10001"));

  // 2 + 1's patterns
  bitset<5> badPattern19(std::string("10011"));
  bitset<5> badPattern22(std::string("10110"));
  bitset<5> badPattern25(std::string("11001"));
  bitset<5> badPattern26(std::string("11010"));
  // 3 + 1's patterns
  bitset<5> badPattern29(std::string("11101"));
  bitset<5> badPattern23(std::string("10111"));
  // 3's patterns
  bitset<5> badPattern7( std::string("00111"));
  bitset<5> badPattern14(std::string("01110"));
  bitset<5> badPattern28(std::string("11100"));
  // 4's patterns
  bitset<5> badPattern15(std::string("01111"));
  bitset<5> badPattern30(std::string("11110"));
  bitset<5> badPattern31(std::string("11111"));

  bool answer = false;

  if(
     //eta patterns are not allowed to be disjoint or spread
     etaPattern != badPattern5  &&
     etaPattern != badPattern10 &&
     etaPattern != badPattern20 &&
     etaPattern != badPattern9  &&
     etaPattern != badPattern18 &&
     etaPattern != badPattern17 &&
     etaPattern != badPattern19 &&
     etaPattern != badPattern22 &&
     etaPattern != badPattern25 &&
     etaPattern != badPattern26 &&
     etaPattern != badPattern29 &&
     etaPattern != badPattern23 &&
     etaPattern != badPattern7  &&
     etaPattern != badPattern14 &&
     etaPattern != badPattern28 &&
     etaPattern != badPattern15 &&
     etaPattern != badPattern30 &&
     etaPattern != badPattern31 &&
     //phi patterns are not allowed to be disjoint
     phiPattern != badPattern5  &&
     phiPattern != badPattern10 &&
     phiPattern != badPattern20 &&
     phiPattern != badPattern9  &&
     phiPattern != badPattern18 &&
     phiPattern != badPattern17 &&
     phiPattern != badPattern19 &&
     phiPattern != badPattern22 &&
     phiPattern != badPattern25 &&
     phiPattern != badPattern26 &&
     phiPattern != badPattern29 &&
     phiPattern != badPattern23 
     ){
    answer = true;
  }
  return answer;
}

/*
bool vetoBit(bitset<4> etaPattern, bitset<4> phiPattern) {

  bitset<4> badPattern5(string("0101"));
  bitset<4> badPattern7(string("0111"));
  bitset<4> badPattern9(string("1001"));
  bitset<4> badPattern10(string("1010"));
  bitset<4> badPattern11(string("1011"));
  bitset<4> badPattern13(string("1101"));
  bitset<4> badPattern14(string("1110"));
  bitset<4> badPattern15(string("1111"));

  bool answer = true;

  if(etaPattern != badPattern5 && etaPattern != badPattern7 &&
     etaPattern != badPattern10 && etaPattern != badPattern11 &&
     etaPattern != badPattern13 && etaPattern != badPattern14 &&
     etaPattern != badPattern15 && phiPattern != badPattern5 &&
     phiPattern != badPattern7 && phiPattern != badPattern10 &&
     phiPattern != badPattern11 && phiPattern != badPattern13 &&
     phiPattern != badPattern14 && phiPattern != badPattern15 &&
     etaPattern != badPattern9 && phiPattern != badPattern9){
    answer = false;
  }
  return answer;

}

uint32_t getHitTowerLocation(uint32_t *et) {
  uint32_t etSum = et[0] + et[1] + et[2] + et[3];
  uint32_t iEtSum =
    (et[0] >> 1)                +  // 0.5xet[0]
    (et[1] >> 1) + et[1]        +  // 1.5xet[1]
    (et[2] >> 1) + (et[2] << 1) +  // 2.5xet[2]
    (et[3] << 2) - (et[3] >> 1) ;  // 3.5xet[3]
  uint32_t iAve = 0xDEADBEEF;
  if(     iEtSum <= etSum) iAve = 0;
  else if(iEtSum <= (etSum << 1)) iAve = 1;
  else if(iEtSum <= (etSum + (etSum << 1))) iAve = 2;
  else iAve = 3;
  return iAve;
}

bitset<4> activeTowerEtaPattern = 0;
for(uint32_t iEta = 0; iEta < nEta; iEta++) {
  bool activeStrip = false;
  for(uint32_t iPhi = 0; iPhi < nPhi; iPhi++) {
    if(activeTower[iEta][iPhi]) activeStrip = true;
  }
  if(activeStrip) activeTowerEtaPattern |= (0x1 << iEta);
 }
bitset<4> activeTowerPhiPattern = 0;
for(uint32_t iPhi = 0; iPhi < nPhi; iPhi++) {
  bool activeStrip = false;
  for(uint32_t iEta = 0; iEta < nEta; iEta++) {
    if(activeTower[iEta][iPhi]) activeStrip = true;
  }
  if(activeStrip) activeTowerPhiPattern |= (0x1 << iPhi);
 }

*/


DEFINE_FWK_MODULE(ClusterProducer);

