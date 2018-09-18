# define basic process
import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("L1CaloClusters")
 

# import standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_v1', '')

# input
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(30))

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/MINIAODSIM/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/C033607A-8E2C-E811-B4EF-0CC47A78A478.root'),
    secondaryFileNames = cms.untracked.vstring(
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/02306B8E-6D2C-E811-9625-0025905B858C.root',
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/163A5C8E-6D2C-E811-B990-0025905B85B2.root',
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/1A7E6D8F-6D2C-E811-8DFB-0025905B85BE.root',
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/5864268A-6D2C-E811-90F5-0025905A60EE.root',
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/627E087F-6D2C-E811-8F81-0CC47A4C8E34.root',
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/64AE2B7F-6D2C-E811-800E-0CC47A4D76A0.root',
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/820C0278-6D2C-E811-B95B-0CC47A4C8F0C.root',
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/BE6A3C80-6D2C-E811-9143-0CC47A4D7668.root',
        'root://cmsxrootd.fnal.gov///store/relval/CMSSW_9_3_7/RelValZTT_14TeV/GEN-SIM-DIGI-RAW/93X_upgrade2023_realistic_v5_2023D17noPU-v2/10000/DA7C408A-6D2C-E811-9BFF-0025905B858A.root'
        ),
    inputCommands = cms.untracked.vstring("keep *", 
        "drop l1tHGCalTowerMapBXVector_hgcalTriggerPrimitiveDigiProducer_towerMap_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT",
        "drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT",
        "drop l1tEMTFHit2016s_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT",
        "drop l1tEMTFTrack2016s_simEmtfDigis__HLT")
     #skipEvents = cms.untracked.uint32(80)
)

process.options = cms.untracked.PSet(

)

# L1 Cluster Producer
process.load("L1Trigger.phase2Demonstrator.L1CaloClusterProducer_cff")
process.L1CaloClusterProducer.debug = cms.untracked.bool(False)
process.L1Clusters = cms.Path(process.L1CaloClusterProducer)

# output module
process.out = cms.OutputModule( "PoolOutputModule",
                                fileName = cms.untracked.string("CaloClusters.root"),
                                fastCloning = cms.untracked.bool( False ),
                                outputCommands = cms.untracked.vstring(#'drop *',
                                                                       'keep *_*_L1Phase2CaloClusters_*', 
                                                                       )
)
process.FEVToutput_step = cms.EndPath(process.out)

# L1 Tau Analyzer
process.load("L1Trigger.phase2Demonstrator.ClusterAnalyzer_cfi")
process.analyzer = cms.Path(process.L1ClusterAnalyzer)


process.TFileService = cms.Service("TFileService", 
   fileName = cms.string("analyzer.root"), 
   closeFileFast = cms.untracked.bool(True)
)


process.schedule = cms.Schedule(process.L1Clusters,process.analyzer,process.FEVToutput_step)


