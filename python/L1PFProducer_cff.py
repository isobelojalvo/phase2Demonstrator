import FWCore.ParameterSet.Config as cms

L1PFProducer = cms.EDProducer("PFObjectProducer",
                              debug           = cms.untracked.bool(False),
                              EoH_cut         = cms.untracked.int32(2),
                              HoE_cut         = cms.untracked.int32(8),
                              L1Clusters      = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                              L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),            

                                       )
