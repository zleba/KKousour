import FWCore.ParameterSet.Config as cms 
import FWCore.ParameterSet.VarParsing as VarParsing

SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
        ignoreTotal = cms.untracked.int32(1) )



options = VarParsing.VarParsing ('analysis')

options.register ('listFile',
        #"/afs/desy.de/user/z/zlebcr/cms/CMSSW_9_3_0/src/KKousour/TopAnalysis/test/farm/fileLists/Aug17/runG.txt", # default value
        "/afs/desy.de/user/z/zlebcr/cms/CMSSW_9_3_0/src/KKousour/TopAnalysis/test/farm/fileLists/QCD_TuneCUETP8M1_13TeV_pythia8/1000to1400.txt", # default value
        VarParsing.VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.VarParsing.varType.string,         # string, int, or float
        "Name of file with list of files")

options.register ('startFile',
        1, # default value
        VarParsing.VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.VarParsing.varType.int,         # string, int, or float
        "Starting file")
options.register ('nFiles',
        1, # default value
        VarParsing.VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.VarParsing.varType.int,         # string, int, or float
        "n Files")



options.parseArguments()

file = open(options.listFile, "r")
#runList = ["root://cms-xrd-global.cern.ch/"+r.rstrip() for r in file.readlines()]
runList = [r.rstrip() for r in file.readlines()]

ff = options.startFile
lf = ff + options.nFiles
curFiles =  runList[ff:lf]

process = cms.Process('myprocess')
process.TFileService=cms.Service("TFileService",fileName=cms.string(options.outputFile))
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '80X_dataRun2_ICHEP16_repro_v0'

print 'current file is', curFiles

if 'mc' in curFiles[0]:
    process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6' #forMC
    print "MC tag used"
else:
    process.GlobalTag.globaltag = '80X_dataRun2_2016LegacyRepro_v4' #Curr for legacy Data
    print "DATA tag used"

##-------------------- Define the source  ----------------------------
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
print "MaxEvents="+str(process.maxEvents)
process.source = cms.Source("PoolSource",
  fileNames = cms.untracked.vstring(
  #      "file:/nfs/dust/cms/user/zlebcr/D2102E03-E415-E611-A4AD-02163E01395E.root"),
  #"root://cms-xrd-global.cern.ch//store/data/Run2017C/JetHT/MINIAOD/PromptReco-v1/000/299/368/00000/189F9B4C-876D-E711-9B34-02163E019BA4.root"),
  #'root://cms-xrd-global.cern.ch//store/data/Run2017C/JetHT/MINIAOD/PromptReco-v1/000/299/368/00000/189F9B4C-876D-E711-9B34-02163E019BA4.root'),
  #"/store/data/Run2016B/JetHT/MINIAOD/PromptReco-v1/000/273/017/00000/8CDA052B-BC19-E611-9CC4-02163E014298.root"),
  #"root://cms-xrd-global.cern.ch//store/data/Run2016H/JetHT/MINIAOD/07Aug17-v1/110000/002C947E-8881-E711-9B13-E0071B7A6890.root"),
  #"root://cms-xrd-global.cern.ch//store/data/Run2016H/JetHT/MINIAOD/07Aug17-v1/110001/F86FCA70-5E7E-E711-9258-0CC47A4C8EEA.root"),
  #"/store/data/Run2016G/JetHT/MINIAOD/23Sep2016-v1/100000/0085E379-F887-E611-AF46-047D7B881D72.root"),
    curFiles),
        #"file:/nfs/dust/cms/user/zlebcr/D2102E03-E415-E611-A4AD-02163E01395E.root"),
#"root://cms-xrd-global.cern.ch//store/data/Run2017C/JetHT/MINIAOD/PromptReco-v1/000/299/368/00000/189F9B4C-876D-E711-9B34-02163E019BA4.root"),
)

task = cms.Task()

############ MET for CHS

def clean_met_(met):
    del met.t01Variation
    del met.t1Uncertainties
    del met.t1SmearedVarsAndUncs
    del met.tXYUncForRaw
    del met.tXYUncForT1
    del met.tXYUncForT01
    del met.tXYUncForT1Smear
    del met.tXYUncForT01Smear

from PhysicsTools.PatAlgos.tools.metTools import addMETCollection

## Raw PAT METs
process.load('RecoMET.METProducers.PFMET_cfi')
process.pfMet.src = cms.InputTag("packedPFCandidates")
task.add(process.pfMet)
addMETCollection(process, labelName='patPFMetCHS', metSource='pfMet') # RAW MET
addMETCollection(process, labelName='patPFMet', metSource='pfMet') # RAW MET
process.patPFMet.addGenMET = False
process.patPFMetCHS.addGenMET = False
## Slimmed METs
from PhysicsTools.PatAlgos.slimming.slimmedMETs_cfi import slimmedMETs
#### CaloMET is not available in MiniAOD
del slimmedMETs.caloMET
# ### CHS
process.slMETsCHS = slimmedMETs.clone()
process.slMETsCHS.src = cms.InputTag("patPFMetCHS")
process.slMETsCHS.rawUncertainties = cms.InputTag("patPFMetCHS") # only central value
task.add(process.slMETsCHS)
clean_met_(process.slMETsCHS)



#############   Format MessageLogger #################
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import selectedPatJets
process.goodJets = selectedPatJets.clone(src='slimmedJets',cut='pt>30 & abs(eta)<2.4')

process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets   = cms.InputTag('slimmedJets')
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

from RecoJets.Configuration.RecoGenJets_cff import ak4GenJets, ak8GenJets

genParticleCollection = 'prunedGenParticles'

from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
from RecoJets.JetProducers.GenJetParameters_cfi import *

process.ak8GenJetsCustom = ak4GenJets.clone(
    src = genParticleCollection,
    rParam = cms.double(0.4),
    jetAlgorithm = cms.string("AntiKt")
)

genParticleCollection = 'prunedGenParticles'
genJetCollection = 'slimmedGenJetsAK8'

from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartons
process.selectedHadronsAndPartons = selectedHadronsAndPartons.clone(
    particles = genParticleCollection
)

from PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi import ak4JetFlavourInfos
process.ak8genJetFlavourInfos = ak4JetFlavourInfos.clone(
    jets = genJetCollection,
    rParam = cms.double(0.4),
)

#only needed if information of the associated b hadrons are required
from PhysicsTools.JetMCAlgos.GenHFHadronMatcher_cff import matchGenBHadron
process.matchGenBHadron = matchGenBHadron.clone(
    genParticles = genParticleCollection,
    jetFlavourInfos = cms.InputTag("ak8genJetFlavourInfos"),
    flavour = cms.int32(5),
    onlyJetClusteredHadrons = cms.bool(True),
    noBBbarResonances = cms.bool(False),
)

##-------------------- User analyzer  --------------------------------
process.ak8 = cms.EDAnalyzer('BoostedTTbarFlatTreeProducer',
  jets             = cms.InputTag('slimmedJetsAK8'),
#  muons            = cms.InputTag('slimmedMuons'),
#  electrons        = cms.InputTag('slimmedElectrons'),
  met1              = cms.InputTag('slimmedMETs'),
  met2              = cms.InputTag('slMETsCHS'),
  met3              = cms.InputTag('slimmedMETsPuppi'),
  vertices         = cms.InputTag('offlineSlimmedPrimaryVertices'),
  rho              = cms.InputTag('fixedGridRhoFastjetAll'),
  candidates       = cms.InputTag('packedPFCandidates'),
  triggerPrescales = cms.InputTag('patTrigger'),
#  xmlFile          = cms.string('boosted_mva_Fisher.weights.xml'),
  nJetsMin         = cms.int32(0),
 # nBJetsMin        = cms.int32(0),
  ptMin            = cms.double(0), 
  minMuPt          = cms.double(0),                              
  minElPt          = cms.double(0),                              
  ptMinLeading     = cms.double(6),   
 # massMin          = cms.double(50),
  htMin            = cms.double(5),
  etaMax           = cms.double(5),
  kinfit           = cms.string('kinFitTtFullHadEvent'),
  fileNames  =       cms.string(curFiles[0]),
  btagMinThreshold = cms.double(0.814),
  btagMin          = cms.double(0.1),                     
  btagMaxThreshold = cms.double(1.1),
  btagger          = cms.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'),
  qgtagger         = cms.InputTag('QGTagger','qgLikelihood'),
  pu               = cms.untracked.string("addPileupInfo"),
  genparticles     = cms.untracked.InputTag('prunedGenParticles'),
  triggerNames     = cms.vstring('HLT_AK8PFJet40_v','HLT_AK8PFJet60_v','HLT_AK8PFJet80_v','HLT_AK8PFJet140_v','HLT_AK8PFJet200_v','HLT_AK8PFJet260_v','HLT_AK8PFJet320_v','HLT_AK8PFJet400_v','HLT_AK8PFJet450_v','HLT_AK8PFJet500_v'),
  triggerResults   = cms.InputTag('TriggerResults','','HLT'),
  triggerObjects  = cms.InputTag("selectedPatTrigger"),
  isMC             = cms.untracked.bool('mc' in curFiles[0]),                              
  genjets          = cms.untracked.InputTag('slimmedGenJetsAK8'),
  SkipEvent       = cms.untracked.vstring('ProductNotFound'),
  GenptMin        = cms.untracked.double(60),
  GenetaMax       = cms.untracked.double(2.5),
  jetFlavourInfos = cms.InputTag("ak8genJetFlavourInfos"),#ak8gen                 
  isPrint         = cms.untracked.bool(True),                           
)

process.ak4 = process.ak8.clone(
    jets            = cms.InputTag('slimmedJets'),
    triggerNames    = cms.vstring('HLT_PFJet40_v','HLT_PFJet60_v','HLT_PFJet80_v','HLT_PFJet140_v','HLT_PFJet200_v','HLT_PFJet260_v','HLT_PFJet320_v','HLT_PFJet400_v','HLT_PFJet450_v','HLT_PFJet500_v',
        'HLT_DiPFJetAve40_v', 'HLT_DiPFJetAve60_v', 'HLT_DiPFJetAve80_v', 'HLT_DiPFJetAve140_v', 'HLT_DiPFJetAve200_v', 'HLT_DiPFJetAve260_v', 'HLT_DiPFJetAve320_v', 'HLT_DiPFJetAve400_v', 'HLT_DiPFJetAve500_v',
       'HLT_DiPFJetAve60_HFJEC_v', 'HLT_DiPFJetAve80_HFJEC_v', 'HLT_DiPFJetAve100_HFJEC_v', 'HLT_DiPFJetAve160_HFJEC_v', 'HLT_DiPFJetAve220_HFJEC_v', 'HLT_DiPFJetAve300_HFJEC_v',),
    genjets         = cms.untracked.InputTag('slimmedGenJets'),
    jetFlavourInfos = cms.InputTag("genJetFlavourInfos"),
)

process.ak4PUPPI = process.ak4.clone(
    jets            = cms.InputTag('slimmedJetsPuppi'),
    genjets         = cms.untracked.InputTag('slimmedGenJetsPuppi'),
    jetFlavourInfos = cms.InputTag("genJetFlavourInfos"),
)




process.load('TopQuarkAnalysis.TopKinFitter.TtFullHadKinFitProducer_cfi')

process.kinFitTtFullHadEvent.jets                = 'goodJets'
process.kinFitTtFullHadEvent.bTagAlgo            = 'pfCombinedInclusiveSecondaryVertexV2BJetTags'
process.kinFitTtFullHadEvent.minBTagValueBJet    = 0.814
process.kinFitTtFullHadEvent.maxBTagValueNonBJet = 0.814
process.kinFitTtFullHadEvent.bTags               = 2
process.kinFitTtFullHadEvent.maxNJets            = 8

#process.ttFilter = cms.EDFilter('AllHadronicPartonFilter',
#  genparticles    = cms.InputTag('prunedGenParticles'),
#  forceTopDecay   = cms.bool(True),
#  forceHiggsDecay = cms.bool(False)
#)

process.kinFitTtFullHadEventNoW = process.kinFitTtFullHadEvent.clone(constraints = cms.vuint32(5))
process.ak8NoW = process.ak8.clone(kinfit = 'kinFitTtFullHadEventNoW')

process.p = cms.Path(
#   process.goodJets * 
#   process.kinFitTtFullHadEvent * 
#   process.ak8*process.ak4 * process.ak4PUPPI 
   process.ak4*
process.ak4PUPPI
)
process.p.associate(task)
process.p.associate(process.patAlgosToolsTask)
