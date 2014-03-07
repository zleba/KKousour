#! /usr/bin/env python
import os
##----------------------------------------------------------------------------------
SAMPLE_LIST = [
  '/QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/QCD_HT-1000ToInf_TuneZ2star_8TeV-madgraph-pythia6/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/ZJetsFullyHadronic_Ht100_Pt50_Pt30_deta22_Mqq200_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/WJetsFullyHadronic_Ht100_Pt50_Pt30_deta22_Mqq200_8TeV-madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/TTJets_MassiveBinDECAY_TuneZ2star_8TeV-madgraph-tauola/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/T_s-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/T_t-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola/Summer12_DR53X-PU_S10_START53_V7A-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/VBF_HToBB_M-115_8TeV-powheg-pythia6_ext/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/VBF_HToBB_M-120_8TeV-powheg-pythia6_ext/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/VBF_HToBB_M-125_8TeV-powheg-pythia6_ext/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/VBF_HToBB_M-130_8TeV-powheg-pythia6_ext/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/VBF_HToBB_M-135_8TeV-powheg-pythia6_ext/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/VBF_HToBB_M-125_8TeV-powheg-pythia8/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/GluGluToHToBB_M-115_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/GluGluToHToBB_M-120_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/GluGluToHToBB_M-125_8TeV-powheg-pythia6_ext/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/GluGluToHToBB_M-130_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/GluGluToHToBB_M-135_8TeV-powheg-pythia6/Summer12_DR53X-PU_S10_START53_V19-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/GluGluToHToBB_M-125_8TeV-powheg-herwigpp/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',
  '/GluGluToHToBB_M-125_8TeV-powheg-pythia8/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0',    
  '/ggH125_HTobb_8TeV_madgraph/Summer12_DR53X-PU_S10_START53_V7C-v1/AODSIM/V5_B/PAT_CMG_V5_17_0'
]
ALIAS        = ['QCD100','QCD250','QCD500','QCD1000','ZJets','WJets',
                'TTJets','T_s-channel','T_t-channel','T_tW-channel','Tbar_s-channel','Tbar_t-channel','Tbar_tW-channel',
                'VBFPowheg115','VBFPowheg120','VBFPowheg125','VBFPowheg130','VBFPowheg135','VBFPowheg125Pythia8',
                'GFPowheg115','GFPowheg120','GFPowheg125','GFPowheg130','GFPowheg135',
                'GFPowheg125Herwig','GFPowheg125Pythia8','GFMadgraph125']
MAX_ETA      = '4.5'
MIN_PT       = ['30','30','30','30','30','30','30','30','30','30','30','30','30','30','30','30',
                '30','30','30','30','30','30','30','30','30','30','30'] 
MIN_DETA     = ['2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2']
MAX_DPHI     = ['2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2','2']
SAVE_SOFT    = ['False','False','False','False','False','False','False','False','False','False','False','False','False',
                'False','False','False','False','False','False','False','False','False','False','False','False','False','False']
SAVE_PROPERTIES= ['False','False','False','False','False','False','False','False','False','False','False','False','False',
                'True','True','True','True','True','False','False','False','False','False','False','False','False','False']
SAVE_PARTONS = ['False','False','False','False','False','False','False','False','False','False','False','False','False',
                'True','True','True','True','True','False','False','False','False','False','False','False','False','False']
CORRECT_CSV = 'False'
CORRECT_QGL = 'False'
##----------------------------------------------------------------------------------
i=0

for ss in SAMPLE_LIST:
  print 'Creating file '+'flat-'+ALIAS[i]+'-cfg.py'
  file = open('flat-'+ALIAS[i]+'-cfg.py','w')
  ##################################################################
  file.write('import FWCore.ParameterSet.Config as cms \n')
  file.write('process = cms.Process(\'myprocess\')\n')
  file.write('process.load(\'FWCore.MessageService.MessageLogger_cfi\')\n')
  file.write('process.TFileService=cms.Service("TFileService",fileName=cms.string(\'flatTree.root\'))\n')
  file.write('from CMGTools.Production.datasetToSource import *\n')
  file.write('##-------------------- Define the source  ----------------------------\n')
  file.write('process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))\n')
  file.write('process.source = datasetToSource(\'cmgtools\',\''+ss+'\',\'cmgTuple_.*.root\')\n')
  file.write('#############   Format MessageLogger #################\n')
  file.write('process.MessageLogger.cerr.FwkReport.reportEvery = 1000\n')
  file.write('##-------------------- User analyzer  --------------------------------\n')
  file.write('process.Hbb = cms.EDAnalyzer(\'VbfHbbFlatTreeProducer\',\n')
  file.write('  jets             = cms.InputTag(\'cmgPFJetSel\'),\n')
  file.write('  softJets         = cms.InputTag(\'ak5SoftPFJetsForVbfHbb\'),\n')
  file.write('  met              = cms.InputTag(\'cmgPFMETRaw\'),\n')
  file.write('  rho              = cms.InputTag(\'kt6PFJets\',\'rho\'),\n')
  file.write('  shiftJES         = cms.double(0.0),\n')
  file.write('  putag            = cms.string(\'full53x\'),\n')
  file.write('  etaMax           = cms.double('+MAX_ETA+'),\n') 
  file.write('  ptMin            = cms.double('+MIN_PT[i]+'),\n')
  file.write('  dEtaMin          = cms.double('+MIN_DETA[i]+'),\n')
  file.write('  dPhiMax          = cms.double('+MAX_DPHI[i]+'),\n')
  file.write('  btagThresholds   = cms.vdouble(0.244,0.679,0.898),\n')
  file.write('  btagger          = cms.string(\'combinedSecondaryVertexBJetTags\'),\n')
  file.write('  saveSoftJets     = cms.untracked.bool('+SAVE_SOFT[i]+'),\n')
  file.write('  saveJetProperties= cms.untracked.bool('+SAVE_PROPERTIES[i]+'),\n')
  file.write('  correctCSV       = cms.untracked.bool('+CORRECT_CSV+'),\n')
  file.write('  correctQGL       = cms.untracked.bool('+CORRECT_QGL+'),\n')
  file.write('  ## MC ########################################\n')
  file.write('  genjets          = cms.untracked.InputTag(\'genJetSel\'),\n')
  file.write('  genparticles     = cms.untracked.InputTag(\'genParticlesPruned\'),\n')
  file.write('  pu               = cms.untracked.string(\'addPileupInfo\'),\n')
  file.write('  savePartons      = cms.untracked.bool('+SAVE_PARTONS[i]+'),\n')
  file.write('  ## trigger ###################################\n')
  file.write('  triggerAlias     = cms.vstring(\'PF\',\'Calo\',\'PF1\',\'PF2\',\'PF3\',\'PF4\',\'Calo1\',\'Calo2\',\'VBF650\',\'VBF700\',\'VBF750\',\'Quad50\',\n')
  file.write('       \'PF40\',\'PF80\',\'DiPFAve40\',\'DiPFAve80\',\'HT200\'),\n')
  file.write('  triggerSelection = cms.vstring(\n')
  file.write('    \'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v* OR HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v* OR HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v* OR HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v*\',\n')
  file.write('    \'HLT_QuadJet75_55_35_20_BTagIP_VBF_v* OR HLT_QuadJet75_55_38_20_BTagIP_VBF_v*\',\n')
  file.write('    \'HLT_QuadPFJet75_55_35_20_BTagCSV_VBF_v*\',\n')
  file.write('    \'HLT_QuadPFJet75_55_38_20_BTagCSV_VBF_v*\',\n')
  file.write('    \'HLT_QuadPFJet78_61_44_31_BTagCSV_VBF_v*\',\n')
  file.write('    \'HLT_QuadPFJet82_65_48_35_BTagCSV_VBF_v*\',\n')
  file.write('    \'HLT_QuadJet75_55_35_20_BTagIP_VBF_v*\',\n')
  file.write('    \'HLT_QuadJet75_55_38_20_BTagIP_VBF_v*\',\n')
  file.write('    \'HLT_DiJet35_MJJ650_AllJets_DEta3p5_VBF_v*\',\n')
  file.write('    \'HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF_v*\',\n')
  file.write('    \'HLT_DiJet35_MJJ750_AllJets_DEta3p5_VBF_v*\',\n')
  file.write('    \'HLT_QuadJet50_v*\',\n')  
  file.write('    \'HLT_PFJet40_v*\',\n')
  file.write('    \'HLT_PFJet80_v*\',\n')
  file.write('    \'HLT_DiPFJetAve40_v*\',\n') 
  file.write('    \'HLT_DiPFJetAve80_v*\',\n')
  file.write('    \'HLT_HT200_v*\',\n') 
  file.write('  ),\n')
  file.write('  triggerConfiguration = cms.PSet(\n')
  file.write('    hltResults            = cms.InputTag(\'TriggerResults\',\'\',\'HLT\'),\n')
  file.write('    l1tResults            = cms.InputTag(\'\'),\n')
  file.write('    daqPartitions         = cms.uint32(1),\n')
  file.write('    l1tIgnoreMask         = cms.bool(False),\n')
  file.write('    l1techIgnorePrescales = cms.bool(False),\n')
  file.write('    throw                 = cms.bool(False)\n')
  file.write('  )\n')
  file.write(')\n')
  file.write('process.p = cms.Path(process.Hbb)')
  file.close()
  i+=1
