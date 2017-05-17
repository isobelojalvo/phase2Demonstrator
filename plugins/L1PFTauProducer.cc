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


#include "L1Trigger/phase2Demonstrator/interface/L1PFTauProducer.hh"

L1PFTauProducer::L1PFTauProducer(const edm::ParameterSet& cfg) :
  debug(            cfg.getUntrackedParameter<bool>("debug", false)),
  input_EoH_cut_(   cfg.getUntrackedParameter<int>("EoH_cut", 2)), // LSB is 0.1 so 2 corresonds to 0.2
  input_HoE_cut_(   cfg.getUntrackedParameter<int>("HoE_cut", 8)), // LSB is 0.1 so 8 corresonds to 0.8
  L1TrackInputTag(  cfg.getParameter<edm::InputTag>("L1TrackInputTag")),
  L1ClustersToken_( consumes< L1CaloClusterCollection >(cfg.getParameter<edm::InputTag>("L1Clusters"))),
  ttTrackToken_(    consumes< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > >(L1TrackInputTag)),
  L1PFToken_(       consumes< L1PFObjectCollection >(cfg.getParameter<edm::InputTag>("l1PFObjects"))),
  L1NeutralToken_(  consumes< std::vector< L1CaloClusterCollection> >(cfg.getParameter<edm::InputTag>("l1Neutrals")) )
{
  //produces three collections of taus, one that uses full FW, one that uses only L1 objects and one that uses only reco objects
  produces< L1PFTauCollection >( "L1Phase2PFTaus" ).setBranchAlias("L1Phase2PFTaus");

}

void L1PFTauProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  std::unique_ptr<L1PFTauCollection> newL1PFTauCollection(new L1PFTauCollection);

  using namespace edm;
  edm::Handle< std::vector<L1CaloCluster> > l1NeutralClusters;
  iEvent.getByToken( L1NeutralToken_, l1NeutralClusters);


  edm::Handle< std::vector<L1CaloCluster> > l1PFChargedCandidates;
  iEvent.getByToken( L1PFToken_, l1PFChargedCandidates);

  L1PFTau tauCands[6]; //6 tau cands possible, each with as many as 3 signal candidates

  int nTauCands = 0;

  //Find taus using charged hadrons
  for(int iCand = 0; iCand < l1PFChargedCandidates->size(); iCand++){

    //Find 1 Prong Candidates
    if(l1PFChargedCandidates->at(iCand).isChargedHadron() == true && l1PFChargedCandidates->at(iCand).p4().Pt() > 6){
      L1PFObject tau_cand[3];

      int n_prongs_found = 0;
      float isolationSum = 0;

      for(int jCand = 0; jCand < l1PFChargedCandidates->size(); jCand++){

	if(jCand <= iCand)
	  continue;

	if(Delta_R(l1PFChargedCandidates->at(jCand).p4().Eta(), l1PFChargedCandidates->at(jCand).p4().Phi(), 
		   l1PFChargedCandidates->at(iCand).p4().Eta(), l1PFChargedCandidates->at(iCand).p4().Phi(), 
		   three_prong_delta_r_)){

	  if(n_prongs_found<3){
	    tau_cand[2] = l1PFChargedCandidates->at(jCand);
	  }

	  if(n_prongs_found<2){
	    tau_cand[1] = l1PFChargedCandidates->at(jCand);
	  }

	}// greater than three_prong_delta_r
	else if(Delta_R(l1PFChargedCandidates->at(jCand).p4().Eta(), l1PFChargedCandidates->at(jCand).p4().Phi(), 
			l1PFChargedCandidates->at(iCand).p4().Eta(), l1PFChargedCandidates->at(iCand).p4().Phi(), 
			isolation_delta_r_)){
	  // Sum the isolation here
	  isolationSum += l1PFChargedCandidates->at(iCand).p4().Pt();
	}// Less than isolation_delta_r
      }

      // Create the 1 Prong Objects
      if(n_prongs_found<3){
	if(l1PFChargedCandidates->at(iCand).p4().Pt() > 15){
	  if(n_prongs_found>1){
	    //Add second prong to Isolation Cone
	    isolationSum += tau_cand[2].P4().Pt();	    
	  }
	  L1PFTau tempL1PFTau;
	  tempL1PFTau.setPtEtaPhiE(l1PFChargedCandidates->at(iCand).p4().Pt(), 
				   l1PFChargedCandidates->at(iCand).p4().Eta(), 
				   l1PFChargedCandidates->at(iCand).p4().Phi(), 
				   l1PFChargedCandidates->at(iCand).p4().Pt());

	  tempL1PFTau.setTauType(0);
	  tempL1PFTau.setRawIso(isolationSum);
	  tempL1PFTau.setRelIso(isolationSum/pt);

	  tauCands[nTauCands] = tempL1PFTau;
	}

      }

    }
  }

  // For each Tau Candidate Create the Strip by grabbing the index of the nearby neutrals



  //find pi0's using electrons/gamma
  //stitch the strips by combining electron and gamma candidates

  //find 1 prong pi0's by combining the 1 prong and the pi0's
  
  
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
