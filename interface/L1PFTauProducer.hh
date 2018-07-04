#ifndef L1PFTAU_PRDC_H
#define L1PFTAU_PRDC_H

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

#include "DataFormats/L1Trigger/interface/L1PFTau.h"
#include "DataFormats/Phase2L1CaloTrig/interface/L1CaloCluster.h"
#include "DataFormats/L1Trigger/interface/L1PFObject.h"
#include "L1Trigger/phase2Demonstrator/interface/triggerGeometryTools.hh"

//#include "L1Trigger/phase2Demonstrator/plugins/strip_alg.hh"
using namespace l1t;

typedef l1t::L1PFTau pftau_t;
typedef l1t::L1PFObject pf_charged_t;
typedef L1CaloCluster cluster_t;

typedef struct{
  float three_prong_seed;
  float three_prong_delta_r;
  float three_prong_max_delta_Z;
  float isolation_delta_r;
  float one_prong_seed;
  float dummy;
  float input_EoH_cut;
  float max_neighbor_strip_dist;
  float min_strip;
  float eg_strip_merge;
} algo_config_t;


typedef struct
{
  float et = 0;
  float eta = 0;
  float phi = 0;
} strip_t;



typedef L1PFTau pftau_t;
typedef L1PFObject pf_charged_t;
typedef L1CaloCluster cluster_t;

using std::vector;
using namespace l1t;

class L1PFTauProducer : public edm::EDProducer {

public:
  explicit L1PFTauProducer(const edm::ParameterSet&);
  ~L1PFTauProducer();

private:

  bool Delta_R( float eta1, float eta2, float phi1, float phi2, float max_dr);
  uint32_t findTheIndexFromReco( float eta, float phi, int iEta_add = 0, int iPhi_add = 0);
  void strip_alg(pftau_t &tau_cand, pf_charged_t electron_grid[5][5],  edm::Handle<std::vector<L1CaloCluster> >&neutral_clusters, algo_config_t algo_config);
  cluster_t find_matched_cluster(edm::Handle<std::vector<L1CaloCluster> >&neutral_clusters, float eta_1, float phi_1);
  cluster_t find_matched_cluster_tower(edm::Handle<std::vector<L1CaloCluster> >&neutral_clusters, unsigned eta, unsigned phi, int eta_side);
  void merge_strip_algo(cluster_t cluster_1, pf_charged_t electron_1, cluster_t cluster_2, pf_charged_t electron_2, strip_t &strip, algo_config_t algo_config);
  float weighted_avg_phi_c_p(cluster_t cluster_1, pf_charged_t pf_charged_1);
  float weighted_avg_eta_c_p(cluster_t cluster_1, pf_charged_t pf_charged_1);
  float weighted_avg_phi_c_c(cluster_t cluster_1, cluster_t cluster_2);
  float weighted_avg_eta_c_c(cluster_t cluster_1, cluster_t cluster_2);
  float weighted_avg_phi_s_s(strip_t strip_1, strip_t strip_2);
  float weighted_avg_eta_s_s(strip_t strip_1, strip_t strip_2);
  float weighted_avg_phi_t_s(pftau_t tau_1, strip_t strip_1);
  float weighted_avg_eta_t_s(pftau_t tau_1, strip_t strip_1);
  float delta_r_cluster(cluster_t cluster_1, cluster_t cluster_2);
  float delta_r_strip(strip_t strip_1, strip_t strip_2);
  /// ///////////////// ///
  /// MANDATORY METHODS ///
  virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
  virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );
  /// ///////////////// ///

  bool debug;
  int input_EoH_cut_;
  int input_HoE_cut_;
  int input_min_n_stubs_;
  int input_max_chi2_; 
  float three_prong_delta_r_;
  float three_prong_max_delta_Z_;
  float isolation_delta_r_;
  edm::InputTag L1TrackInputTag;
  edm::EDGetTokenT< L1CaloClusterCollection > L1ClustersToken_;
  edm::EDGetTokenT< std::vector< TTTrack< Ref_Phase2TrackerDigi_ > > > ttTrackToken_;
  edm::EDGetTokenT< L1PFObjectCollection > L1PFToken_;
  edm::EDGetTokenT< L1CaloClusterCollection > L1NeutralToken_;

};


#endif
