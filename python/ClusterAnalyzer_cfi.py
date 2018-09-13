import FWCore.ParameterSet.Config as cms

L1ClusterAnalyzer = cms.EDAnalyzer('phase2L1Clusters',
                               L1Clusters       = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                               genParticles     = cms.InputTag("genParticles", "", "HLT"),
                               packedCandidates = cms.InputTag("packedPFCandidates","","RECO"),
                               ecalTPGsBarrel = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT")
)
