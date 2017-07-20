import FWCore.ParameterSet.Config as cms

L1PFTauProducer = cms.EDProducer("L1PFTauProducer",
                                 debug           = cms.untracked.bool(False),
                                 EoH_cut         = cms.untracked.int32(2),
                                 HoE_cut         = cms.untracked.int32(8),
                                 three_prong_dr  = cms.untracked.double(0.15),
                                 L1Clusters      = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                                 L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),            
                                 L1PFObjects     = cms.InputTag("L1PFProducer","L1PFObjects"),
                                 L1Neutrals      = cms.InputTag("L1PFProducer", "L1PFObjects")
                                 )
