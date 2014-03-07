#! /usr/bin/env python
import os
##----------------------------------------------------------------------------------
SAMPLE_LIST = [
  '/MultiJet/Run2012A-22Jan2013-v1/AOD/PAT_CMG_V5_17_0',
  '/BJetPlusX/Run2012B-22Jan2013-v1/AOD/PAT_CMG_V5_17_0',
  '/BJetPlusX/Run2012C-22Jan2013-v1/AOD/PAT_CMG_V5_17_0',
  '/BJetPlusX/Run2012D-22Jan2013-v1/AOD/PAT_CMG_V5_17_0',
  '/VBF1Parked/Run2012B/22Jan2013/PAT_CMG_V5_17_0',
  '/VBF1Parked/Run2012C/22Jan2013/PAT_CMG_V5_17_0',
  '/VBF1Parked/Run2012D/22Jan2013/PAT_CMG_V5_17_0'
       
]
ALIAS        = ['MultiJetA','BJetPlusXB','BJetPlusXC','BJetPlusXD','VBF1ParkedB','VBF1ParkedC','VBF1ParkedD']
MAX_ETA      = '4.5'
MIN_PT       = '30' 
MIN_DETA     = '2'
MAX_DPHI     = '2'
SAVE_SOFT    = 'False'
SAVE_PROPERTIES = 'False'
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
  file.write('from CMGTools.Common.Tools.applyJSON_cff import *\n')
  file.write('from CMGTools.Production.eostools import *\n')
  file.write('##-------------------- Define the source  ----------------------------\n')
  file.write('process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))\n')
  if ss[:5] == '/VBF1':
    file.write('process.source = cms.Source(\'PoolSource\',\n')
    file.write('    fileNames = cms.untracked.vstring(ls_EOS(\'/store/cmst3/group/vbfhbb/CMG/'+ss+'\'))\n')
    file.write(')\n')
  else:
    file.write('process.source = datasetToSource(\'cmgtools\',\''+ss+'\',\'cmgTuple_.*.root\')\n')
  file.write('applyJSON(process,\'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt\')\n')
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
  file.write('  ptMin            = cms.double('+MIN_PT+'),\n')
  file.write('  dEtaMin          = cms.double('+MIN_DETA+'),\n')
  file.write('  dPhiMax          = cms.double('+MAX_DPHI+'),\n')
  if ss[:5] == '/VBF1':
    file.write('  forceVBF         = cms.untracked.bool(True),\n')
  else:
    file.write('  forceNOM         = cms.untracked.bool(True),\n')
  file.write('  btagThresholds   = cms.vdouble(0.244,0.679,0.898),\n')
  file.write('  btagger          = cms.string(\'combinedSecondaryVertexBJetTags\'),\n')
  file.write('  saveSoftJets     = cms.untracked.bool('+SAVE_SOFT+'),\n')
  file.write('  saveJetProperties= cms.untracked.bool('+SAVE_PROPERTIES+'),\n')
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
  file.write('\n')
  file.write('process.json = cms.EDAnalyzer(\'JSONProducer\',\n')
  file.write('   filename = cms.string(\'json.txt\')\n')
  file.write(')\n')
  file.write('process.p = cms.Path(process.json * process.Hbb)')
  file.close()
  i+=1
