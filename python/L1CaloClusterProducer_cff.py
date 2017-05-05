import FWCore.ParameterSet.Config as cms

L1CaloClusterProducer = cms.EDProducer("ClusterProducer",
                                       debug         = cms.untracked.bool(False),
                                       useECalEndcap = cms.bool(False),
                                       ecalRecHitEB  = cms.InputTag("ecalRecHit","EcalRecHitsEB","RECO"),
                                       ecalRecHitEE  = cms.InputTag("ecalRecHit","EcalRecHitsEE","RECO"),
                                       hcalDigis     = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT"),
                                       )
