from CRABClient.UserUtilities import config, getUsernameFromSiteDB

config = config()
config.General.requestName = 'QCD_Pt_170to300'
config.General.transferOutputs = True
config.General.transferLogs = False
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'flat-QCD-cfg.py'
config.JobType.maxJobRuntimeMin = 5000
config.Data.inputDataset = '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISpring15DR74-Asympt50ns_MCRUN2_74_V9A-v2/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/group/cmst3/user/kkousour/ttbar/'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
