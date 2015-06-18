import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'PHYS14_25_V1'
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root",
     "/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00D3EAF1-3174-E411-A5B2-0025904B144E.root",
     "/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02EF3EFC-0475-E411-A9DB-002590DB9166.root",
     "/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/0434E222-1C75-E411-B4D4-0025907FD34C.root")
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
  btagMinThreshold = cms.double(0.814),
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

process.ttFilter = cms.EDFilter('AllHadronicPartonFilter',
  genparticles    = cms.InputTag('prunedGenParticles'),
  forceTopDecay   = cms.bool(True),
  forceHiggsDecay = cms.bool(False)
)

process.kinFitTtFullHadEventNoW = process.kinFitTtFullHadEvent.clone(constraints = cms.vuint32(5))
process.hadtopNoW = process.hadtop.clone(kinfit = 'kinFitTtFullHadEventNoW')

process.p = cms.Path(
   process.ttFilter * 
   process.goodJets * 
   process.QGTagger * 
   #process.kinFitTtFullHadEventNoW *
   #process.hadtopNoW *
   process.kinFitTtFullHadEvent * 
   process.hadtop
   
)








