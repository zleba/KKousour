from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()
config.General.requestName = 'JetHT_Run2016C-PromptReco-v2'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'flat-data-cfg.py'
config.JobType.maxJobRuntimeMin = 2750
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset = '/JetHT/Run2016C-PromptReco-v2/MINIAOD'
config.Data.inputDBS = 'global'
config.JobType.inputFiles = ['Spring16_25nsV3_DATA.db']
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.lumiMask = 'Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt'
config.Data.outLFNDirBase = '/store/group/cmst3/user/kkousour/ttH/'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
