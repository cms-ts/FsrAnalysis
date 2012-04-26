import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("JetValidation")

###################
##### Loading what we need!
###################

from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger(process,sequence='patDefaultSequence',hltProcess = '*')
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *
from RecoJets.JetProducers.FastjetParameters_cfi import *
from RecoJets.JetProducers.ak5TrackJets_cfi import *
from RecoJets.JetProducers.GenJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

process.load("CommonTools.ParticleFlow.pfElectrons_cff")
process.load("CommonTools.ParticleFlow.pfMuons_cff")
process.load("CommonTools.ParticleFlow.ParticleSelectors.pfSortByType_cff")
process.load("CommonTools.ParticleFlow.pfNoPileUp_cff")
process.load("CommonTools.ParticleFlow.TopProjectors.pfNoElectron_cfi")
process.load("CommonTools.ParticleFlow.TopProjectors.pfNoMuon_cfi")
process.load("CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi")
##-------------------- Import the JEC services -----------------------
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")
process.load("JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff")
##process.load("JetMETCorrections.Configuration.JetCorrectionProducers_cff")

##-------------------- Import the Jet RECO modules -----------------------
process.load('RecoJets.Configuration.RecoPFJets_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.load("MagneticField.Engine.uniformMagneticField_cfi") 
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("RecoJets.Configuration.GenJetParticles_cff")
process.load("RecoJets.Configuration.RecoGenJets_cff")

##-------------------- Turn-on the FastJet density calculation -----------------------
process.kt6PFJets.doRhoFastjet = True

##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
process.kt6PFJets.doAreaFastjet = True
process.ak5PFJets.doAreaFastjet = True

# to compute FastJet rho to correct isolation (note: EtaMax restricted to 2.5)
process.kt6PFJetsForIsolation = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True)
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)


process.options   = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True) )

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(True),
                                     makeTriggerResults=cms.untracked.bool(True),
                                     )

process.GlobalTag.globaltag = 'MC_44_V5D::All'

####################
#### Files
###################

# readFiles = cms.untracked.vstring()
# readFiles.extend([
# #"file:/gpfs/cms/data/2011/r9test/pythiaZ2tunesroot/FEF7EE7B-8780-E011-837F-E41F131816A8.root",
# "file:/gpfs/cms/data/2011/SynchTest/DYJetsToLL_TuneZ2_PU_S6_START44.root",
#     ])

#No FSR
# readFiles = cms.untracked.vstring()
# readFiles.extend([    
# "file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_nolepfsr_cff_py_RAW2DIGI_L1Reco_RECO_5_1_gjk.root",
# "file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_nolepfsr_cff_py_RAW2DIGI_L1Reco_RECO_2_1_zEi.root",
# "file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_nolepfsr_cff_py_RAW2DIGI_L1Reco_RECO_3_1_fL9.root",
# "file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_nolepfsr_cff_py_RAW2DIGI_L1Reco_RECO_1_1_BaG.root",
# "file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_nolepfsr_cff_py_RAW2DIGI_L1Reco_RECO_4_1_Rhn.root",
# "file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_nolepfsr_cff_py_RAW2DIGI_L1Reco_RECO_6_1_5wf.root",
# "file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_nolepfsr_cff_py_RAW2DIGI_L1Reco_RECO_7_1_jrE.root",
# "file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_nolepfsr_cff_py_RAW2DIGI_L1Reco_RECO_9_1_AeS.root",
# "file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_nolepfsr_cff_py_RAW2DIGI_L1Reco_RECO_8_1_tEX.root",
#     ])

# #Tauola
readFiles = cms.untracked.vstring()
readFiles.extend([
    #"file:/gpfs/grid/srm/cms/store/data/Fall11/DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola/AODSIM/PU_S6-START44_V5-v1/0001/B48B1A68-460A-E111-88DF-485B39800BAB.root"
"file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff_py_RAW2DIGI_L1Reco_RECO_2_1_vpH.root",
#"file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff_py_RAW2DIGI_L1Reco_RECO_1_1_5GF.root",
#"file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff_py_RAW2DIGI_L1Reco_RECO_6_1_JKT.root",
#"file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff_py_RAW2DIGI_L1Reco_RECO_7_1_BcD.root",
#"file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff_py_RAW2DIGI_L1Reco_RECO_4_1_LZx.root",
#"file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff_py_RAW2DIGI_L1Reco_RECO_3_1_s7f.root",
#"file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff_py_RAW2DIGI_L1Reco_RECO_8_1_kAs.root",
#"file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff_py_RAW2DIGI_L1Reco_RECO_5_1_lkN.root",
#"file:/gpfs/cms/data/2012/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph/20120410-RAW2DIGI_L1Reco_RECO-CMSSW_4_4_2_patch8/Hadronizer_MgmMatchTuneZ2_7TeV_madgraph_tauola_cff_py_RAW2DIGI_L1Reco_RECO_9_1_0zt.root",
    ])
    
process.MessageLogger.cerr.FwkReport  = cms.untracked.PSet(
     #reportEvery = cms.untracked.int32(10),
     reportEvery = cms.untracked.int32(500),
 )

process.maxEvents = cms.untracked.PSet(  input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = readFiles,
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )


####################
#### Trigger...
###################

trigger2011v3 = cms.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9", "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10")

trigger2011RunB= cms.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8", "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9", "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10")

trigger2011v1  = cms.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3","HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3","HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3","HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3","HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2","HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v3","HLT_Ele45_CaloIdVT_TrkIdT_v3","HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v4","HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v4","HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v4")

trigger2010    = cms.vstring("HLT_Ele17_CaloIdl_Ele8_CaloIsoIdL_CaloIsoVL_v3","HLT_Ele15_SW_L1R","HLT_Ele15_SW_CaloEleId_L1R","HLT_Ele17_SW_CaloEleId_L1R","HLT_Ele17_SW_TightEleId_L1R","HLT_Ele17_SW_TightEleId_L1R_v2","HLT_Ele17_SW_TightEleId_L1R_v3","HLT_Photon10_L1R","HLT_Photon15_L1R","HTL_Photon15_Cleaned_L1R")

alltriggers    = cms.vstring() # In this way, the HLT string is empty and it will trigger every event

trigger2011v2 = cms.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8")

# from dav HLT analysis
triggersMay10Jul05 = cms.vstring("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5","HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6")

triggersAug05 = cms.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7","HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v6")

triggersOct03 = cms.vstring("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7","HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v6","HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8","HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v7","HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_v1")


####################
#### Lepton Selection
###################


process.Selection = cms.EDFilter('ZpatFilterPf',
                                 electronCollection = cms.InputTag("patElectronsWithTrigger"),
                                 triggerCollectionTag = cms.InputTag("TriggerResults","","HLT"),
                                 UseCombinedPrescales = cms.bool(False),
                                 doTheHLTAnalysis = cms.bool(False),
                                 removePU=  cms.bool(False),
                                 #TriggerNames = triggersMay10Jul05+triggersAug05+triggersOct03+trigger2011v2+trigger2010,
                                 TriggerNames = trigger2011v3,
                                 secondEleEnThrhold   = cms.double(20.0),
                                 firstEleEnThrhold    = cms.double(20.0),
                                 lowZmassLimit        = cms.double(71.0),
                                 highZmassLimit       = cms.double(111.0),
                                 maxEtaForElectron    = cms.double(2.4),
                                 )

process.goodEPair = cms.EDProducer('pfAnalyzer',
                                   electronCollection = cms.InputTag("patElectronsWithTrigger"),
                                   pflowEleCollection = cms.untracked.InputTag("pfIsolatedElectrons"),
                                   removePU=  cms.bool(False),
                                   secondEleEnThrhold   = cms.double(20.0),
                                   firstEleEnThrhold    = cms.double(20.0),
                                   lowZmassLimit        = cms.double(71.0),
                                   highZmassLimit       = cms.double(111.0),
                                   maxEtaForElectron    = cms.double(2.4),
                                   )

process.fsrValidation = cms.EDAnalyzer('fsrValidation',
                                       goodEPair = cms.InputTag("goodEPair"),
                                       usingMC = cms.untracked.bool(True),
                                       usingPF = cms.untracked.bool(True),
                                       #useSecondChoice = cms.untracked.bool(True),
                                       deltaRConeGen         = cms.double(0.1),
                                       deltaRConeGenloose    = cms.double(0.2),                             
                                       )
                                   
####################
#### HLT Analysis, MC reweight, and other stuff
###################

process.demo = cms.EDProducer('HistoProducer',
                              electronCollection = cms.InputTag('patElectronsWithTrigger'),
                              triggerCollection = cms.InputTag("TriggerResults","","HLT"),
                              UseCombinedPrescales = cms.bool(False),
                              #TriggerNames = triggersMay10Jul05+triggersAug05+triggersOct03,
                              TriggerNames = trigger2011v3,
                              removePU=  cms.bool(True),
                              usingMC=  cms.bool(True),
                              doTheHLTAnalysis = cms.bool(False),
                              VertexCollectionTag = cms.InputTag('offlinePrimaryVertices'),              
                              TotalNEventTag = cms.vstring('TotalEventCounter'),
                              WhichRun = cms.string("Run2011AB"), ##Select which datasets you wonna use to reweight
                              RootuplaName = cms.string("treeVJ_"),
                              eventWeightsCollection= cms.string("EventWeight"),
                              giveEventWeightEqualToOne= cms.bool(True)
)



######################
#                    #
#     pfNoPileUP     #
#                    #
######################

process.pfPileUp.Vertices = 'goodOfflinePrimaryVertices'  # recipe 15th March JEC
process.pfPileUp.checkClosestZVertex = cms.bool(False) # recipe 15th March JEC
process.pfPileUp.PFCandidates = cms.InputTag("particleFlow")
process.pfNoPileUp.bottomCollection = cms.InputTag("particleFlow")

######################
#                    #
#    pfElectrons     #
#                    #
######################

process.patElectrons.useParticleFlow=True
#process.pfAllElectrons.src = "particleFlow"
process.pfAllElectrons.src = "pfNoPileUp"
process.isoValElectronWithNeutral.deposits[0].deltaR = 0.4
process.isoValElectronWithCharged.deposits[0].deltaR = 0.4
process.isoValElectronWithPhotons.deposits[0].deltaR = 0.4
process.pfIsolatedElectrons.isolationCut = 0.2 


######################
#                    #
#  TRG MATCHING -ON- #
#                    #
######################

### Standard trigger matching:
process.eleTriggerMatchHLT = cms.EDProducer( "PATTriggerMatcherDRLessByR",
                                             src     = cms.InputTag( "patElectrons" ),
                                             matched = cms.InputTag( "patTrigger"),
                                             matchedCuts = cms.string('(path("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*",0,0) && filter("hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter")) || (path("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*",0,0) && filter("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter"))'),
                                             maxDPtRel = cms.double( 5 ),
                                             maxDeltaR = cms.double( 0.3 ),
                                             resolveAmbiguities    = cms.bool( True ),
                                             resolveByMatchQuality = cms.bool( True )
                                             )

process.patElectronsWithTrigger = cms.EDProducer("PATTriggerMatchElectronEmbedder",
                                                    src     = cms.InputTag("patElectrons"),
                                                    matches = cms.VInputTag(cms.InputTag('eleTriggerMatchHLT'))
                                                 )

switchOnTriggerMatching( process, ['eleTriggerMatchHLT' ],sequence ='patDefaultSequence', hltProcess = '*' )


######################
#                    #
#  PROCESS OUT       #
#                    #
######################

process.out.fileName = cms.untracked.string('test-filtering.root')
process.out.outputCommands =  cms.untracked.vstring(
    'drop *',
    )
process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('JetValidation'))

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('fsrValidation.root')
                                   )


#####################
#                   #
#    Counting       #
#                   #
#####################


process.TotalEventCounter = cms.EDProducer("EventCountProducer")

######################
#                    #
#  SEQUENCE          #
#                    #
######################

process.ToolInizialization = cms.Path(
    process.kt6PFJetsForIsolation*
    process.kt6PFJets*
    process.ak5PFJets*
    process.ak5PFJetsL1FastL2L3*
    process.goodOfflinePrimaryVertices*
    process.pfNoPileUpSequence*
    process.pfAllNeutralHadrons*
    process.pfAllChargedHadrons*
    process.pfAllPhotons*
    process.pfElectronSequence*
    process.patTrigger*
    process.patDefaultSequence*
    process.eleTriggerMatchHLT*
    process.patElectronsWithTrigger
    #process.goodElec
    #process.goodEPair*
   #  process.pfNoElectron*
#     process.ak5PFJetsRC*
#     process.ak5PFchsJetsRCL1FastL2L3*
#     process.ak5PFJetsOLD*
#     process.ak5PFJetsOLDL1FastL2L3*
#     process.ak5PFJetsPU*
#     process.ak5PFchsJetsPUL1FastL2L3
    )


# process.TAPAnalysisWP80 = cms.Path(
#     process.goodOfflinePrimaryVertices*
#     process.trgmatchPatElectronsEle8*
#     process.TAPwp80
#     )

# process.TAPAnalysisWP80newHE = cms.Path(
#     process.goodOfflinePrimaryVertices*
#     process.trgmatchPatElectronsEle8*
#     process.TAPwp80newHE
#     )

# process.TAPAnalysisHLTele8NOTele17 = cms.Path(
#     process.goodOfflinePrimaryVertices*
#     process.trgmatchPatElectronsNOTEle17*
#     process.trgmatchPatElectronsEle8NOTEle17*
#     process.TAPhltele8NOTele17
#     )

# process.TAPAnalysisHLTele17 = cms.Path(
#     process.goodOfflinePrimaryVertices*
#     process.trgmatchPatElectronsEle17*
#     process.TAPhltele17
#     )

# process.TAPAnalysisRECO = cms.Path(
#     process.goodOfflinePrimaryVertices*
#     process.trgmatchPatElectronsEle8*
#     process.TAPreco
#     )

process.JetValidation = cms.Path(
    process.TotalEventCounter* 
    process.eleTriggerMatchHLT*
    process.patElectronsWithTrigger*  
    process.goodOfflinePrimaryVertices*   
    process.Selection*
    process.demo*
    process.goodEPair*
    process.fsrValidation
    )

#####################
#                   #
#    Outpath        #
#                   #
#####################

process.outpath = cms.EndPath(
        #process.out
        )
