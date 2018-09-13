// -*- C++ -*-
//
// Package:    L1Trigger/phase2L1Clusters
// Class:      phase2L1Clusters
// 
/**\class phase2L1Clusters phase2L1Clusters.cc L1Trigger/phase2L1Clusters/plugins/phase2L1Clusters.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Isobel Ojalvo
//         Created:  Fri, 26 May 2017 15:44:30 GMT
//
//


// system include files
#include <memory>

// Math Include
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <TLorentzVector.h>
#include <memory>
#include <math.h>
#include <vector>
#include <list>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "L1Trigger/phase2Demonstrator/interface/triggerGeometryTools.hh"
#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"

//Vertex and gen particle
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"


//#include "L1Trigger/phase2L1Clusters/plugins/helpers.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Phase2L1CaloTrig/interface/L1CaloCluster.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.


class phase2L1Clusters : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit phase2L1Clusters(const edm::ParameterSet&);
      ~phase2L1Clusters();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
  typedef std::vector<reco::GenParticle> GenParticleCollectionType;

  const CaloSubdetectorGeometry * ebGeometry;
  edm::ESHandle<CaloGeometry> caloGeometry_;

  edm::EDGetTokenT< L1CaloClusterCollection > L1ClustersToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_;
  edm::EDGetTokenT<std::vector<pat::PackedCandidate> > PackedCands_;
  edm::EDGetTokenT<EcalEBTrigPrimDigiCollection> ecalTPGBToken_;
  edm::InputTag genSrc_;
  edm::InputTag L1TrackInputTag;

  TH1F* nEvents;  
  TTree* efficiencyTree;

  double genPt, genEta, genPhi;
  double recoPt, recoEta, recoPhi; 
  int run, lumi, event;

  TH1F* track_pt;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
phase2L1Clusters::phase2L1Clusters(const edm::ParameterSet& cfg):
  L1ClustersToken_( consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Clusters"))),
  PackedCands_(     consumes< std::vector<pat::PackedCandidate> >(cfg.getParameter<edm::InputTag>("packedCandidates"))),
  ecalTPGBToken_(   consumes<EcalEBTrigPrimDigiCollection>(cfg.getParameter<edm::InputTag>("ecalTPGsBarrel"))),
  genSrc_ ((        cfg.getParameter<edm::InputTag>( "genParticles")))
{
   //now do what ever initialization is needed
  usesResource("TFileService");
  genToken_ =     consumes<std::vector<reco::GenParticle> >(genSrc_);
  
  edm::Service<TFileService> fs;
  efficiencyTree = fs->make<TTree>("efficiencyTree", "Crystal cluster individual crystal pt values");
  efficiencyTree->Branch("run",    &run,     "run/I");
  efficiencyTree->Branch("lumi",   &lumi,    "lumi/I");
  efficiencyTree->Branch("event",  &event,   "event/I");
  
  efficiencyTree->Branch("genPt",  &genPt,   "genPt/D");
  efficiencyTree->Branch("genEta", &genEta,   "genEta/D");
  efficiencyTree->Branch("genPhi", &genPhi,   "genPhi/D");

  efficiencyTree->Branch("recoPt",  &recoPt,   "recoPt/D");
  efficiencyTree->Branch("recoEta", &recoEta,   "recoEta/D");
  efficiencyTree->Branch("recoPhi", &recoPhi,   "recoPhi/D");

  nEvents              = fs->make<TH1F>( "nEvents"  , "nEvents", 2,  0., 1. );
  track_pt             = fs->make<TH1F>( "track_pt"  , "p_{t}", 200,  0., 200. );

}


phase2L1Clusters::~phase2L1Clusters()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)


}


//
// member functions
//

// ------------ method called for each event  ------------
void
phase2L1Clusters::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  edm::Handle< std::vector<L1CaloCluster> > l1CaloClusters;
  iEvent.getByToken( L1ClustersToken_, l1CaloClusters);
  
  run   = iEvent.id().run();
  lumi  = iEvent.id().luminosityBlock();
  event = iEvent.id().event();
  
  edm::Handle< std::vector<L1CaloCluster> > l1Clusters;
  iEvent.getByToken( L1ClustersToken_, l1Clusters);
  
  Handle<std::vector<pat::PackedCandidate> > packedcands;
  iEvent.getByToken(PackedCands_, packedcands);
  
  //control plot for ecal crystals
  edm::Handle<EcalEBTrigPrimDigiCollection> ecaltpgCollection;
  iEvent.getByToken( ecalTPGBToken_, ecaltpgCollection);

  edm::ESHandle<CaloGeometry> caloGeometryHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeometryHandle);
  const CaloGeometry* caloGeometry_ = caloGeometryHandle.product();  
  ebGeometry = caloGeometry_->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

  for(unsigned int i = 0; i < l1CaloClusters->size(); i++){
    L1CaloCluster cluster = l1CaloClusters->at(i);
  }

  for(auto& tpg : *ecaltpgCollection.product())
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
	}
    }

}


// ------------ method called once each job just before starting event loop  ------------
void 
phase2L1Clusters::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
phase2L1Clusters::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
phase2L1Clusters::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(phase2L1Clusters);
