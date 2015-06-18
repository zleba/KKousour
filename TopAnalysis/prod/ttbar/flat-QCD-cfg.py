import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "/store/mc/RunIISpring15DR74/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/10000/1C6830A2-410D-E511-AD22-848F69FD4EFB.root",
    "/store/mc/RunIISpring15DR74/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/10000/48B3F74E-4C0D-E511-A1A4-3417EBE6475F.root",
    "/store/mc/RunIISpring15DR74/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/10000/4A2005EF-640F-E511-B197-0025905A6056.root",
    "/store/mc/RunIISpring15DR74/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v2/10000/58D9DE5E-3D0D-E511-937B-00266CFAE268.root"
  )
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
process.hadtop = cms.EDAnalyzer('TTbarFlatTreeProducer',
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
  btagger          = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
  qgtagger         = cms.InputTag('QGTagger','qgLikelihood'),
  pu               = cms.untracked.string("addPileupInfo"),
  genparticles     = cms.untracked.InputTag('prunedGenParticles'),
  
  triggerNames     = cms.vstring('HLT_PFJet260_v','HLT_PFHT900_v'),
  triggerResults   = cms.InputTag('TriggerResults','','HLT')
)

process.load('TopQuarkAnalysis.TopKinFitter.TtFullHadKinFitProducer_cfi')

process.kinFitTtFullHadEvent.jets                = 'goodJets'
process.kinFitTtFullHadEvent.bTagAlgo            = 'pfCombinedInclusiveSecondaryVertexV2BJetTags'
process.kinFitTtFullHadEvent.minBTagValueBJet    = 0.814
process.kinFitTtFullHadEvent.maxBTagValueNonBJet = 0.814
process.kinFitTtFullHadEvent.bTags               = 2
process.kinFitTtFullHadEvent.maxNJets            = 8

process.kinFitTtFullHadEventNoBtag = process.kinFitTtFullHadEvent.clone(bTags = 0)
process.hadtopNoBtag = process.hadtop.clone(kinfit = "kinFitTtFullHadEventNoBtag",nBJetsMin=0)

process.p = cms.Path(
   process.goodJets * 
   process.QGTagger * 
   process.kinFitTtFullHadEvent * 
   process.kinFitTtFullHadEventNoBtag *
   process.hadtop *
   process.hadtopNoBtag
)








