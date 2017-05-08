import FWCore.ParameterSet.Config as cms

L1CaloClusterProducer = cms.EDProducer("ClusterProducer",
                                       debug          = cms.untracked.bool(False),
                                       ecalTPGsBarrel = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
                                       hcalTPGs       = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT"),
                                       )
