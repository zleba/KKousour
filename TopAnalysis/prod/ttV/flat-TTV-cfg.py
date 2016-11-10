import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "/store/mc/RunIISpring16MiniAODv2/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/60000/14E4CB3E-D61B-E611-BEE0-141877412793.root"
    #"/store/mc/RunIISpring16MiniAODv2/TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/50000/5AF9DA75-2E1D-E611-AD17-D4AE527EE013.root"
    #"/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext3-v1/00000/000B9244-4B27-E611-91D2-7845C4FC3C6B.root"
    )
)
#############   Format MessageLogger #################
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#############   JEC #################
#process.load("CondCore.DBCommon.CondDBCommon_cfi")
#from CondCore.DBCommon.CondDBSetup_cfi import *
#process.jec = cms.ESSource("PoolDBESSource",
#      DBParameters = cms.PSet(
#        messageLevel = cms.untracked.int32(0)
#        ),
#      timetype = cms.string('runnumber'),
#      toGet = cms.VPSet(
#        cms.PSet(
#            record = cms.string('JetCorrectionsRecord'),
#            tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK4PFchs'),
#            label  = cms.untracked.string('AK4PFchs')
#        ),  
#      ),
#      connect = cms.string('sqlite:Fall15_25nsV2_MC.db')
#)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

#--- first re-apply JEC from the GT -------------------------
process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.patJetCorrFactorsReapplyJEC = process.updatedPatJetCorrFactors.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet','L2Relative','L3Absolute'],
  payload = 'AK4PFchs' 
) 

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.patJetsReapplyJEC = process.updatedPatJets.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)

#--- clean jets from leptons -----------------------------------
process. cleanedJets = cms.EDProducer('JetCleanedProducer',
  jets        = cms.InputTag('patJetsReapplyJEC'),
  rho         = cms.InputTag('fixedGridRhoFastjetAll'),
  muons       = cms.InputTag('slimmedMuons'),
  electrons   = cms.InputTag('slimmedElectrons'),
  vertices    = cms.InputTag('offlineSlimmedPrimaryVertices'),
  minMuPt     = cms.double(20),
  minElPt     = cms.double(20)
)

#--- finally define the good jets -------------------------------
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
process.goodJets = selectedPatJets.clone(src='cleanedJets',cut='pt>30 & abs(eta)<2.4')

##-------------------- User analyzers  --------------------------------
process.load('TopQuarkAnalysis.TopKinFitter.TtFullHadKinFitProducer_cfi')

process.kinFitTtFullHadEvent.jets                = 'goodJets'
process.kinFitTtFullHadEvent.bTagAlgo            = 'pfCombinedInclusiveSecondaryVertexV2BJetTags'
process.kinFitTtFullHadEvent.minBTagValueBJet    = 0.8
process.kinFitTtFullHadEvent.maxBTagValueNonBJet = 0.8
process.kinFitTtFullHadEvent.bTags               = 1
process.kinFitTtFullHadEvent.maxNJets            = 8

process.resolved = cms.EDAnalyzer('TTVFlatTreeProducer',
  jets             = cms.InputTag('goodJets'),
  muons            = cms.InputTag('slimmedMuons'),
  electrons        = cms.InputTag('slimmedElectrons'),
  met              = cms.InputTag('slimmedMETs'),
  vertices         = cms.InputTag('offlineSlimmedPrimaryVertices'),
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  ptMin            = cms.double(30),
  probMin          = cms.double(-10),
  etaMax           = cms.double(2.4),
  minMuPt          = cms.double(20),
  minElPt          = cms.double(20),
  kinfit           = cms.string('kinFitTtFullHadEvent'),
  btagMin          = cms.double(0.8),
  btagger          = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
  xmlFileTTW       = cms.string("ttW_mva_BDT.weights.xml"),
  xmlFileTTZ       = cms.string("ttZ_mva_BDT.weights.xml"),
  isMC             = cms.untracked.bool(True),
  saveWeights      = cms.untracked.bool(False),
  triggerNames     = cms.vstring(
    'HLT_IsoMu20_v'
  ),
  triggerResults   = cms.InputTag('TriggerResults','','HLT'),
  triggerPrescales = cms.InputTag('patTrigger')
)

process.eventCounter = cms.EDAnalyzer("EventCounter")

process.reapplyjec = cms.Sequence(
   process.patJetCorrFactorsReapplyJEC +
   process.patJetsReapplyJEC
)

process.selectjets = cms.Sequence(
   process.cleanedJets +
   process.goodJets
)

process.kinfit = cms.Sequence(
   process.kinFitTtFullHadEvent
)

process.resolvedanalyzer = cms.Sequence(
   process.resolved 
)

process.p = cms.Path(
   process.eventCounter *
   process.reapplyjec *     
   process.selectjets *
   process.kinfit *
   process.resolvedanalyzer
)








