from WMCore.Configuration import Configuration
config = Configuration()

config.section_("User")
config.User.voGroup = 'dcms'

config.section_("General")
config.General.requestName = 'JetData-March2017-Run2016B_TopAnalysis-v3'
config.General.workArea = 'JetData-March2017-Run2016B_TopAnalysis-v3'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'flatData-TTJets-cfg.py'
#config.JobType.psetName = 'flatData-TTJets-FullTriggers_cfg.py'
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB  = 3000

config.section_("Data")
#config.Data.inputDataset = '/JetHT/Run2016G-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016C-03Feb2017-v1/MINIAOD'
config.Data.inputDataset = '/JetHT/Run2016B-03Feb2017_ver2-v2/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016D-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016H-03Feb2017_ver2-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016E-03Feb2017-v1/MINIAOD'
#config.Data.inputDataset = '/JetHT/Run2016F-03Feb2017-v1/MINIAOD'
config.Data.splitting = 'LumiBased' #LumiBased' 
config.Data.unitsPerJob = 10
config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.runRange = '271036-284044' 
config.Data.outputDatasetTag = 'CRAB3_JetData-October2017-Run2016B_TopAnalysis-v3'

config.section_("Site")
config.Site.storageSite = "T2_DE_DESY"

