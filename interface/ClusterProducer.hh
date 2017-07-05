#include <vector>
#include <iostream>

#ifndef L1CLUSTER_PRDC_H
#define L1CLUSTER_PRDC_H

#include <iostream>
#include <math.h>
#include <vector>
#include <list>
#include <TLorentzVector.h>

////////////////////
// FRAMEWORK HEADERS
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"

//L1 TPG Legacy
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CondFormats/L1TObjects/interface/L1CaloHcalScale.h"
#include "CondFormats/DataRecord/interface/L1CaloHcalScaleRcd.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"

//Vertex and gen particle
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"

//Geometry
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HcalTowerAlgo/interface/HcalTrigTowerGeometry.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

#include <bitset>
using std::bitset;

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "L1Trigger/phase2Demonstrator/interface/triggerGeometryTools.hh"

using std::vector;

class ClusterProducer : public edm::EDProducer {

public:
  explicit ClusterProducer(const edm::ParameterSet&);
  ~ClusterProducer();

private:



  /// ///////////////// ///
  /// MANDATORY METHODS ///
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );
  
  bool debug;
  bool useECalEndcap;

  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPGBToken_;
  edm::EDGetTokenT< edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcalSrc_;

  edm::ESHandle<CaloGeometry> caloGeometry_;
  const CaloSubdetectorGeometry * ebGeometry;
  const CaloSubdetectorGeometry * hbGeometry;
  edm::ESHandle<HcalTopology> hbTopology;
  const HcalTopology * hcTopology_;

  struct ecalCrystal_t{
    TLorentzVector p4;
    float iEta;
    float iPhi;
    EBDetId id;
  } ;

  struct hcalTPG_t{
    TLorentzVector p4;
    int iEta;
    int iPhi;
    HcalTrigTowerDetId id;
  } ;
  
  void getEcalCrystals(edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgs, 
		       vector<ecalCrystal_t> &ecalCrystals);

  float findEcalCrystal(int tEta, int tPhi, int cEta, int cPhi, 
			vector<ecalCrystal_t> ecalCrystals, 
			ecalCrystal_t &foundCrystal);

  void getHcalTPGs( edm::Handle<edm::SortedCollection<HcalTriggerPrimitiveDigi> > hcaltpgCollection, 
		    edm::ESHandle<L1CaloHcalScale> &hcalScale, 
		    vector<hcalTPG_t> &allHcalTPGs);

  void clusterAlgoMaxInternet(float crystals[25], 
			      unsigned int &maxCrystalEta, 
			      unsigned int &maxCrystalPhi);

  void clusterAlgoMax(float crystals[5][5], 
		      unsigned int &maxCrystalEta, 
		      unsigned int &maxCrystalPhi);

  bool pi0BitSet(bitset<5> etaPattern, bitset<5> phiPattern);

  int getNeighbor(L1CaloCluster centralCluster, int dEta, int dPhi);
  bool checkAndMergeCluster(L1CaloCluster &centralCluster, L1CaloCluster &newClusters);

  int TPGEtaRange(int ieta){
    int iEta = 0;
    // So here, -28 becomes 0.  -1 be comes 27.  +1 becomes 28. +28 becomes 55.
    // And we have mapped [-28, -1], [1, 28] onto [0, 55]   
    if(ieta < 0)
      iEta = ieta + 28;
    else if(ieta > 0)
      iEta = ieta + 27;
    if(ieta==0){
      std::cout<<"Error! ieta is 0, ieta: "<<ieta<<" iEta "<<iEta<<std::endl;
      exit(0);
    }
    return iEta;
  }


};


#endif
