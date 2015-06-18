import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'PHYS14_25_V1'
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/08B36E8F-5E7F-E411-9D5A-002590200AE4.root",
    "/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/00000/FC4E6E16-5C7F-E411-8843-002590200AE4.root",
    "/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/10000/629CE4D5-687F-E411-BF71-001E673969FA.root",
    "/store/mc/Phys14DR/TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola/MINIAODSIM/PU20bx25_tsg_PHYS14_25_V1-v2/20000/14587980-CB7E-E411-A0F4-001E67397701.root")
)
#############   Format MessageLogger #################
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
process.goodJets = selectedPatJets.clone(src='slimmedJets',cut='pt>30 & abs(eta)<2.4')

process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag('slimmedJets')
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

##-------------------- User analyzer  --------------------------------
process.hadtop = cms.EDAnalyzer('TTHFlatTreeProducer',
  jets             = cms.InputTag('slimmedJets'),
  muons            = cms.InputTag('slimmedMuons'),
  electrons        = cms.InputTag('slimmedElectrons'),
  met              = cms.InputTag('slimmedMETs'),
  vertices         = cms.InputTag('offlineSlimmedPrimaryVertices'),
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  nJetsMin         = cms.int32(6),
  nBJetsMin        = cms.int32(2),
  ptMin            = cms.double(30),
  htMin            = cms.double(400),
  etaMax           = cms.double(2.4),
  kinfit           = cms.string('kinFitTtFullHadEvent'),
  btagMinThreshold = cms.double(0.814),#0.423/0.814/0.941/ vs 0.244/0.679/0.898
  btagMaxThreshold = cms.double(1.1),
  btagger          = cms.string('combinedInclusiveSecondaryVertexV2BJetTags'),
  qgtagger         = cms.InputTag('QGTagger','qgLikelihood'),
  pu               = cms.untracked.string("addPileupInfo"),
  genparticles     = cms.untracked.InputTag('prunedGenParticles'),
  
  triggerNames     = cms.vstring('HLT_PFJet260_v','HLT_PFHT900_v'),
  triggerResults   = cms.InputTag('TriggerResults','','HLT')
)

process.load('TopQuarkAnalysis.TopKinFitter.TtFullHadKinFitProducer_cfi')

process.kinFitTtFullHadEvent.jets                = 'goodJets'
process.kinFitTtFullHadEvent.bTagAlgo            = 'combinedInclusiveSecondaryVertexV2BJetTags'
process.kinFitTtFullHadEvent.minBTagValueBJet    = 0.814
process.kinFitTtFullHadEvent.maxBTagValueNonBJet = 0.814
process.kinFitTtFullHadEvent.bTags               = 2
process.kinFitTtFullHadEvent.maxNJets            = 8

process.ttHFilter = cms.EDFilter('AllHadronicPartonFilter',
  genparticles    = cms.InputTag('prunedGenParticles'),
  forceTopDecay   = cms.bool(True),
  forceHiggsDecay = cms.bool(True)
)

process.p = cms.Path(
   process.ttHFilter * 
   process.goodJets *
   process.QGTagger * 
   process.kinFitTtFullHadEvent * 
   process.hadtop
)








