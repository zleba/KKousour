import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_miniAODv2'
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "/store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0_ext3-v1/00000/000B9244-4B27-E611-91D2-7845C4FC3C6B.root"
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
#        cms.PSet(
#            record = cms.string('JetCorrectionsRecord'),
#            tag    = cms.string('JetCorrectorParametersCollection_Fall15_25nsV2_MC_AK8PFchs'),
#            label  = cms.untracked.string('AK8PFchs')
#        ) 
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

process.patJetCorrFactorsReapplyJECAK8 = process.updatedPatJetCorrFactors.clone(
  src = cms.InputTag("slimmedJetsAK8"),
  levels = ['L1FastJet','L2Relative','L3Absolute'],
  payload = 'AK8PFchs' 
) 

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.patJetsReapplyJECAK8 = process.updatedPatJets.clone(
  jetSource = cms.InputTag("slimmedJetsAK8"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK8"))
)

#--- then smear MC jets to match the JER in data --------------
process.smearedJets = cms.EDProducer('JetShiftProducer',
  jets        = cms.InputTag('patJetsReapplyJEC'),
  rho         = cms.InputTag('fixedGridRhoFastjetAll'),
  payload     = cms.untracked.string('AK4PFchs'),
  resSFFile   = cms.untracked.string('Summer15_25nsV6_MC_SF_AK4PFchs.txt'),
  shiftJES    = cms.untracked.double(0.0),
  shiftJER    = cms.untracked.double(0.0),
  doSmear     = cms.untracked.bool(True),
  doShift     = cms.untracked.bool(False)
)

process.smearedJetsUp   = process.smearedJets.clone(shiftJER = 1.0)
process.smearedJetsDown = process.smearedJets.clone(shiftJER = -1.0)

process.smearedJetsAK8 = cms.EDProducer('JetShiftProducer',
  jets        = cms.InputTag('patJetsReapplyJECAK8'),
  rho         = cms.InputTag('fixedGridRhoFastjetAll'),
  payload     = cms.untracked.string('AK8PFchs'),
  resSFFile   = cms.untracked.string('Summer15_25nsV6_MC_SF_AK4PFchs.txt'),
  shiftJES    = cms.untracked.double(0.0),
  shiftJER    = cms.untracked.double(0.0),
  doSmear     = cms.untracked.bool(True),
  doShift     = cms.untracked.bool(False)
)

process.smearedJetsAK8Up   = process.smearedJetsAK8.clone(shiftJER = 1.0)
process.smearedJetsAK8Down = process.smearedJetsAK8.clone(shiftJER = -1.0)

#--- JES variations -------------------------------------
process.shiftedJetsUp = cms.EDProducer('JetShiftProducer',
  jets        = cms.InputTag('patJetsReapplyJEC'),
  rho         = cms.InputTag('fixedGridRhoFastjetAll'),
  payload     = cms.untracked.string('AK4PFchs'),
  shiftJES    = cms.untracked.double(1.0),
  doShift     = cms.untracked.bool(True)
)
process.shiftedJetsDown = process.shiftedJetsUp.clone(shiftJES = -1.0)

process.shiftedJetsAK8Up = cms.EDProducer('JetShiftProducer',
  jets        = cms.InputTag('patJetsReapplyJECAK8'),
  rho         = cms.InputTag('fixedGridRhoFastjetAll'),
  payload     = cms.untracked.string('AK8PFchs'),
  shiftJES    = cms.untracked.double(1.0),
  doShift     = cms.untracked.bool(True)
)

process.shiftedJetsAK8Down = process.shiftedJetsAK8Up.clone(shiftJES = -1.0)

#--- finally define the good jets -------------------------------
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
process.goodJets            = selectedPatJets.clone(src='patJetsReapplyJEC',cut='pt>30 & abs(eta)<2.4')
process.goodSmearedJets     = selectedPatJets.clone(src='smearedJets',cut='pt>30 & abs(eta)<2.4')
process.goodSmearedJetsUp   = selectedPatJets.clone(src='smearedJetsUp',cut='pt>30 & abs(eta)<2.4')
process.goodSmearedJetsDown = selectedPatJets.clone(src='smearedJetsDown',cut='pt>30 & abs(eta)<2.4')
process.goodShiftedJets     = selectedPatJets.clone(src='shiftedJets',cut='pt>30 & abs(eta)<2.4')
process.goodShiftedJetsUp   = selectedPatJets.clone(src='shiftedJetsUp',cut='pt>30 & abs(eta)<2.4')
process.goodShiftedJetsDown = selectedPatJets.clone(src='shiftedJetsDown',cut='pt>30 & abs(eta)<2.4')

process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag('goodJets')
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

process.QGTaggerSmeared     = process.QGTagger.clone(srcJets = 'goodSmearedJets')
process.QGTaggerSmearedUp   = process.QGTagger.clone(srcJets = 'goodSmearedJetsUp')
process.QGTaggerSmearedDown = process.QGTagger.clone(srcJets = 'goodSmearedJetsDown')
process.QGTaggerShiftedUp   = process.QGTagger.clone(srcJets = 'goodShiftedJetsUp')
process.QGTaggerShiftedDown = process.QGTagger.clone(srcJets = 'goodShiftedJetsDown')

##-------------------- User analyzers  --------------------------------
process.boosted = cms.EDAnalyzer('BoostedTTbarFlatTreeProducer',
  jets             = cms.InputTag('patJetsReapplyJECAK8'),
  muons            = cms.InputTag('slimmedMuons'),
  electrons        = cms.InputTag('slimmedElectrons'),
  met              = cms.InputTag('slimmedMETs'),
  candidates       = cms.InputTag('packedPFCandidates'),
  vertices         = cms.InputTag('offlineSlimmedPrimaryVertices'),
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  massMin          = cms.double(50),
  ptMin            = cms.double(200),
  ptMinLeading     = cms.double(200),
  etaMax           = cms.double(2.4),
  minMuPt          = cms.double(20),
  minElPt          = cms.double(20),
  btagMin          = cms.double(0.8),
  btagger          = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
  xmlFile          = cms.string('boosted_mva_Fisher_new.weights.xml'),
  isMC             = cms.untracked.bool(True),
  saveWeights      = cms.untracked.bool(False),
  triggerNames     = cms.vstring(
    'HLT_AK8PFJet360_TrimMass30_v',
    'HLT_PFJet200_v'
  ),
  triggerResults   = cms.InputTag('TriggerResults','','HLT'),
  triggerPrescales = cms.InputTag('patTrigger')
)

process.boostedSmeared     = process.boosted.clone(jets = 'smearedJetsAK8', saveWeights = False)
process.boostedSmearedUp   = process.boosted.clone(jets = 'smearedJetsAK8Up', saveWeights = False)
process.boostedSmearedDown = process.boosted.clone(jets = 'smearedJetsAK8Down', saveWeights = False)
process.boostedShiftedUp   = process.boosted.clone(jets = 'shiftedJetsAK8Up', saveWeights = False)
process.boostedShiftedDown = process.boosted.clone(jets = 'shiftedJetsAK8Down', saveWeights = False)

process.load('TopQuarkAnalysis.TopKinFitter.TtFullHadKinFitProducer_cfi')

process.kinFitTtFullHadEvent.jets                = 'goodJets'
process.kinFitTtFullHadEvent.bTagAlgo            = 'pfCombinedInclusiveSecondaryVertexV2BJetTags'
process.kinFitTtFullHadEvent.minBTagValueBJet    = 0.8
process.kinFitTtFullHadEvent.maxBTagValueNonBJet = 0.8
process.kinFitTtFullHadEvent.bTags               = 2
process.kinFitTtFullHadEvent.maxNJets            = 8

process.kinFitTtFullHadEventSmeared     = process.kinFitTtFullHadEvent.clone(jets = 'goodSmearedJets')
process.kinFitTtFullHadEventSmearedUp   = process.kinFitTtFullHadEvent.clone(jets = 'goodSmearedJetsUp')
process.kinFitTtFullHadEventSmearedDown = process.kinFitTtFullHadEvent.clone(jets = 'goodSmearedJetsDown')
process.kinFitTtFullHadEventShiftedUp   = process.kinFitTtFullHadEvent.clone(jets = 'goodShiftedJetsUp')
process.kinFitTtFullHadEventShiftedDown = process.kinFitTtFullHadEvent.clone(jets = 'goodShiftedJetsDown')

process.resolved = cms.EDAnalyzer('TTbarFlatTreeProducer',
  jets             = cms.InputTag('goodJets'),
  muons            = cms.InputTag('slimmedMuons'),
  electrons        = cms.InputTag('slimmedElectrons'),
  met              = cms.InputTag('slimmedMETs'),
  vertices         = cms.InputTag('offlineSlimmedPrimaryVertices'),
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  nJetsMin         = cms.int32(6),
  nBJetsMin        = cms.int32(2),
  ptMin            = cms.double(30),
  htMin            = cms.double(400),
  probMin          = cms.double(0.0),
  etaMax           = cms.double(2.4),
  kinfit           = cms.string('kinFitTtFullHadEvent'),
  btagMin          = cms.double(0.8),
  btagger          = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
  qgtagger         = cms.InputTag('QGTagger','qgLikelihood'),
  isMC             = cms.untracked.bool(True),
  saveWeights      = cms.untracked.bool(True),
  triggerNames     = cms.vstring(
    'HLT_PFHT450_SixJet40_PFBTagCSV0p72_v',
    'HLT_PFHT450_SixJet40_v',
    'HLT_PFHT350_v'
  ),
  triggerResults   = cms.InputTag('TriggerResults','','HLT'),
  triggerPrescales = cms.InputTag('patTrigger')
)

process.resolvedSmeared     = process.resolved.clone(jets = 'goodSmearedJets',qgtagger = cms.InputTag('QGTaggerSmeared','qgLikelihood'), kinfit = "kinFitTtFullHadEventSmeared", saveWeights = False)
process.resolvedSmearedUp   = process.resolved.clone(jets = 'goodSmearedJetsUp',qgtagger = cms.InputTag('QGTaggerSmearedUp','qgLikelihood'), kinfit = "kinFitTtFullHadEventSmearedUp", saveWeights = False)
process.resolvedSmearedDown = process.resolved.clone(jets = 'goodSmearedJetsDown',qgtagger = cms.InputTag('QGTaggerSmearedDown','qgLikelihood'), kinfit = "kinFitTtFullHadEventSmearedDown", saveWeights = False)
process.resolvedShiftedUp   = process.resolved.clone(jets = 'goodShiftedJetsUp',qgtagger = cms.InputTag('QGTaggerShiftedUp','qgLikelihood'), kinfit = "kinFitTtFullHadEventShiftedUp", saveWeights = False)
process.resolvedShiftedDown = process.resolved.clone(jets = 'goodShiftedJetsDown',qgtagger = cms.InputTag('QGTaggerShiftedDown','qgLikelihood'), kinfit = "kinFitTtFullHadEventShiftedDown", saveWeights = False)

process.eventCounter = cms.EDAnalyzer("EventCounter")

process.reapplyjec = cms.Sequence(
   #process.patJetCorrFactorsReapplyJEC +
   process.patJetCorrFactorsReapplyJECAK8 + 
   process.patJetsReapplyJECAK8
   #process.patJetsReapplyJEC
)

process.smearjets = cms.Sequence(
   process.smearedJetsAK8 +
   process.smearedJetsAK8Up +
   process.smearedJetsAK8Down +
   process.smearedJets +
   process.smearedJetsUp +
   process.smearedJetsDown
)

process.shiftjets = cms.Sequence(
   process.shiftedJetsAK8Up +
   process.shiftedJetsAK8Down +
   process.shiftedJetsUp +
   process.shiftedJetsDown
)

process.selectjets = cms.Sequence(
   process.goodJets +
   process.goodSmearedJets + 
   process.goodSmearedJetsUp + 
   process.goodSmearedJetsDown + 
   process.goodShiftedJetsUp + 
   process.goodShiftedJetsDown
)

process.qgtagging = cms.Sequence(
   process.QGTagger +
   process.QGTaggerSmeared +
   process.QGTaggerSmearedUp +
   process.QGTaggerSmearedDown +
   process.QGTaggerShiftedUp +
   process.QGTaggerShiftedDown
)

process.kinfit = cms.Sequence(
   process.kinFitTtFullHadEvent + 
   process.kinFitTtFullHadEventSmeared +
   process.kinFitTtFullHadEventSmearedUp +
   process.kinFitTtFullHadEventSmearedDown +
   process.kinFitTtFullHadEventShiftedUp +
   process.kinFitTtFullHadEventShiftedDown
)

process.boostedanalyzer = cms.Sequence(
   process.boosted 
   #process.boostedSmeared +
   #process.boostedSmearedUp +
   #process.boostedSmearedDown +
   #process.boostedShiftedUp +
   #process.boostedShiftedDown
)

process.resolvedanalyzer = cms.Sequence(
   process.resolved +
   process.resolvedSmeared +
   process.resolvedSmearedUp +
   process.resolvedSmearedDown +
   process.resolvedShiftedUp +
   process.resolvedShiftedDown
)

process.p = cms.Path(
   process.eventCounter *
   process.reapplyjec *    
   #process.smearjets *
   #process.shiftjets * 
   #process.selectjets *
   #process.qgtagging *
   process.boostedanalyzer
   #process.kinfit *
   #process.resolvedanalyzer
)








