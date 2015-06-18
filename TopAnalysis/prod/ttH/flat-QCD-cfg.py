import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'PHYS14_25_V1'
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "/store/mc/Phys14DR/QCD_HT-500To1000_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1_ext1-v1/00000/12FCC45C-5D78-E411-BEB2-00261894388B.root",
    "/store/mc/Phys14DR/QCD_HT-500To1000_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1_ext1-v1/00000/1603BD77-5D78-E411-8344-002618943926.root",
    "/store/mc/Phys14DR/QCD_HT-500To1000_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1_ext1-v1/00000/186824AE-5D78-E411-9926-00259059391E.root",
    "/store/mc/Phys14DR/QCD_HT-500To1000_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1_ext1-v1/00000/1A086774-5D78-E411-9210-0026189438F4.root")
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

process.kinFitTtFullHadEventNoBtag = process.kinFitTtFullHadEvent.clone(bTags = 0)
process.hadtopNoBtag = process.hadtop.clone(btagMinThreshold = 0, btagMaxThreshold = 0.814, kinfit = 'kinFitTtFullHadEventNoBtag')
process.hadtopL = process.hadtop.clone(btagMinThreshold = 0.423, btagMaxThreshold = 0.814, kinfit = '')  

process.p = cms.Path(
   process.goodJets * 
   process.QGTagger * 
   #process.kinFitTtFullHadEventNoBtag *
   #process.hadtopNoBtag *
   process.kinFitTtFullHadEvent * 
   process.hadtop *
   process.hadtopL
)








