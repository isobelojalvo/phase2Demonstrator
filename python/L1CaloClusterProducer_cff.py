import FWCore.ParameterSet.Config as cms

L1CaloClusterProducer = cms.EDProducer("ClusterProducer",
                                       hcalDigis = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT"),
                                       )
