import FWCore.ParameterSet.Config as cms

process = cms.Process("GenLevelInfoAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
            input = cms.untracked.int32(-1)
            )

process.source = cms.Source("PoolSource",
            fileNames = cms.untracked.vstring(
                #'file:test.root'
                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1_ext1-v2/120000/00EDC2A7-74B1-324A-A65F-97072DCE3283.root'
                #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv9/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v16_L1v1-v1/30000/224CDDD8-1627-6343-B3AA-D9DFE9CF4FAF.root',
                #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16NanoAODv9/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_mcRun2_asymptotic_v17-v1/80000/AA7BEF20-415D-E24B-A9B3-CD1C03A7D8D6.root'
#                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/2810000/41C0B183-5307-6D4D-9240-623F38F64370.root',
#                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/2810000/C7D2492E-1209-7449-A305-DA5581ADEDB5.root',
#                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/2810000/69B88828-345C-F74D-ABB9-50C51779CD4A.root',
#                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/2810000/5BDB111A-DBE6-3449-8748-4DCE1AD27CF2.root',
#                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/2810000/41DB8960-2567-264A-9282-B9F5FAC8D8A3.root',
#                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/2810000/516831CB-DF99-CC40-9307-47FDF9C03F2E.root',
#                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/2810000/C5D4B698-004D-7445-B1A9-73947ABD98D7.root',
#                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/80000/B42385A2-8CE9-9943-B605-BC4D554775EC.root',
#                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/HPlusCharm_4FS_MuRFScaleDynX0p50_HToWWTo2L2Nu_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v2/80000/6566A299-6A2C-124D-8CCA-CA7C15258692.root',
                            )
            )

process.GenLevelInfoAnalyzer = cms.EDAnalyzer('GenLevelInfoAnalyzer',
            genParticles = cms.InputTag('prunedGenParticles'),
)

process.p = cms.Path(process.GenLevelInfoAnalyzer)
