import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '80X_dataRun2_Prompt_ICHEP16JEC_v0'
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
    "/store/data/Run2016C/JetHT/MINIAOD/PromptReco-v2/000/275/782/00000/0AA1D765-D43C-E611-8937-02163E01463F.root"
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
#            tag    = cms.string('JetCorrectorParametersCollection_Spring16_25nsV3_DATA_AK4PFchs'),
#            label  = cms.untracked.string('AK4PFchs')
#        ),
#        cms.PSet(
#            record = cms.string('JetCorrectionsRecord'),
#            tag    = cms.string('JetCorrectorParametersCollection_Spring16_25nsV3_DATA_AK8PFchs'),
#            label  = cms.untracked.string('AK8PFchs')
#        )  
#      ),
#      connect = cms.string('sqlite:Spring16_25nsV3_DATA.db')
#)
## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
#process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.patJetCorrFactorsReapplyJEC = process.updatedPatJetCorrFactors.clone(
  src = cms.InputTag("slimmedJets"),
  levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual'],
  payload = 'AK4PFchs' 
) 

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.patJetsReapplyJEC = process.updatedPatJets.clone(
  jetSource = cms.InputTag("slimmedJets"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
)

process.patJetCorrFactorsReapplyJECAK8 = process.updatedPatJetCorrFactors.clone(
  src = cms.InputTag("slimmedJetsAK8"),
  levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual'],
  payload = 'AK8PFchs' 
) 

process.load("PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff")
process.patJetsReapplyJECAK8 = process.updatedPatJets.clone(
  jetSource = cms.InputTag("slimmedJetsAK8"),
  jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK8"))
)

from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
process.goodJets = selectedPatJets.clone(src='patJetsReapplyJEC',cut='pt>30 & abs(eta)<2.4')

############# QGL #################
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag('goodJets')
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

##-------------------- User analyzer  --------------------------------
process.boosted = cms.EDAnalyzer('BoostedTTbarFlatTreeProducer',
  jets             = cms.InputTag('patJetsReapplyJECAK8'),
  muons            = cms.InputTag('slimmedMuons'),
  electrons        = cms.InputTag('slimmedElectrons'),
  met              = cms.InputTag('slimmedMETs'),
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
  triggerNames     = cms.vstring(
    'HLT_AK8PFJet360_TrimMass30_v',
    'HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20_v',
    'HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v',
    'HLT_AK8DiPFJet250_200_TrimMass30_v',
    'HLT_AK8DiPFJet280_200_TrimMass30_v',
    'HLT_AK8PFJet140_v',
    'HLT_PFJet140_v',
    'HLT_PFJet200_v'
  ),
  triggerResults   = cms.InputTag('TriggerResults','','HLT'),
  triggerPrescales = cms.InputTag('patTrigger')
)

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
  
  triggerNames     = cms.vstring(
   'HLT_PFHT450_SixJet40_PFBTagCSV0p72_v',
   'HLT_PFHT450_SixJet40_v',
   'HLT_PFHT350_v'
  ),
  triggerResults   = cms.InputTag('TriggerResults','','HLT'),
  triggerPrescales = cms.InputTag('patTrigger')
)

process.resolvedNoBtag = process.resolved.clone(nBJetsMin = 0, kinfit = "kinFitTtFullHadEventNoBtag")

process.load('TopQuarkAnalysis.TopKinFitter.TtFullHadKinFitProducer_cfi')

process.kinFitTtFullHadEvent.jets                = 'goodJets'
process.kinFitTtFullHadEvent.jetCorrectionLevel  = 'L2L3Residual'
process.kinFitTtFullHadEvent.bTagAlgo            = 'pfCombinedInclusiveSecondaryVertexV2BJetTags'
process.kinFitTtFullHadEvent.minBTagValueBJet    = 0.8
process.kinFitTtFullHadEvent.maxBTagValueNonBJet = 0.8
process.kinFitTtFullHadEvent.bTags               = 2
process.kinFitTtFullHadEvent.maxNJets            = 8
process.kinFitTtFullHadEvent.jetEnergyResolutionScaleFactors = cms.vdouble(1.061,1.088,1.106,1.126,1.343)
process.kinFitTtFullHadEvent.jetEnergyResolutionEtaBinning = cms.vdouble(0.0,0.8,1.3,1.9,2.5,-1)

process.kinFitTtFullHadEventNoBtag = process.kinFitTtFullHadEvent.clone(bTags = 0)

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.filter = hlt.hltHighLevel.clone(
   HLTPaths = ['HLT_PFHT450_SixJet40_v*'],
   throw = False
)

process.p = cms.Path(
   #process.patJetCorrFactorsReapplyJEC +
   process.patJetCorrFactorsReapplyJECAK8 +
   #process.patJetsReapplyJEC +
   process.patJetsReapplyJECAK8 +
   #process.goodJets + 
   #process.QGTagger +
   process.boosted
   #process.kinFitTtFullHadEvent + 
   #process.resolved +
   #process.filter +
   #process.kinFitTtFullHadEventNoBtag +
   #process.resolvedNoBtag
)








