import FWCore.ParameterSet.Config as cms 
process = cms.Process('myprocess')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.TFileService=cms.Service("TFileService",fileName=cms.string('flatTree.root'))
from CMGTools.Production.datasetToSource import *
##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.source = datasetToSource('cmgtools','/VBF_HToBB_M-120_8TeV-powheg-pythia6_ext/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0','cmgTuple_.*.root')
#############   Format MessageLogger #################
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
##-------------------- User analyzer  --------------------------------
process.Hbb = cms.EDAnalyzer('VbfHbbFlatTreeProducer',
  jets             = cms.InputTag('cmgPFJetSel'),
  softJets         = cms.InputTag('ak5SoftPFJetsForVbfHbb'),
  met              = cms.InputTag('cmgPFMETRaw'),
  rho              = cms.InputTag('kt6PFJets','rho'),
  shiftJES         = cms.double(0.0),
  putag            = cms.string('full53x'),
  ptMin            = cms.double(20),
  dEtaMin          = cms.double(2),
  btagThresholds   = cms.vdouble(0.244,0.679,0.898),
  btagger          = cms.string('combinedSecondaryVertexBJetTags'),
  saveSoftJets     = cms.untracked.bool(False),
  ## MC ########################################
  genjets          = cms.untracked.InputTag('genJetSel'),
  genparticles     = cms.untracked.InputTag('genParticlesPruned'),
  pu               = cms.untracked.string('addPileupInfo'),
  savePartons      = cms.untracked.bool(True),
  ## trigger ###################################
  triggerAlias     = cms.vstring('PF','Calo','PF1','PF2','PF3','PF4','Calo1','Calo2','VBF650','VBF700','VBF750','Quad50',
       'PF40','PF80','DiPFAve40','DiPFAve80','HT200'),
  triggerSelection = cms.vstring(
    'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v* OR HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v* OR HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v* OR HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v*',
    'HLT_QuadJet75_55_35_20_BTagIP_VBF_v* OR HLT_QuadJet75_55_38_20_BTagIP_VBF_v*',
    'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v*',
    'HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v*',
    'HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v*',
    'HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v*',
    'HLT_QuadJet75_55_35_20_BTagIP_VBF_v*',
    'HLT_QuadJet75_55_38_20_BTagIP_VBF_v*',
    'HLT_DiJet35_MJJ650_AllJets_DEta3p5_VBF_v*',
    'HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF_v*',
    'HLT_DiJet35_MJJ750_AllJets_DEta3p5_VBF_v*',
    'HLT_QuadJet50_v*',
    'HLT_PFJet40_v*',
    'HLT_PFJet80_v*',
    'HLT_DiPFJetAve40_v*',
    'HLT_DiPFJetAve80_v*',
    'HLT_HT200_v*',
  ),
  triggerConfiguration = cms.PSet(
    hltResults            = cms.InputTag('TriggerResults','','HLT'),
    l1tResults            = cms.InputTag(''),
    daqPartitions         = cms.uint32(1),
    l1tIgnoreMask         = cms.bool(False),
    l1techIgnorePrescales = cms.bool(False),
    throw                 = cms.bool(False)
  )
)
process.p = cms.Path(process.Hbb)