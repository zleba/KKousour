#! /usr/bin/env python

crabSubmitFile = open("SubmitCrabJobsData.csh","w")
crabSubmitFile.write("#!/bin/tcsh\n")

SAMPLES = [
  "/JetHT/Run2016B-PromptReco-v2/MINIAOD",
  "/JetHT/Run2016C-PromptReco-v2/MINIAOD",
  "/JetHT/Run2016D-PromptReco-v2/MINIAOD",
  "/JetHT/Run2016E-PromptReco-v2/MINIAOD",
  "/JetHT/Run2016F-PromptReco-v2/MINIAOD"
]

for ss in SAMPLES:
  tag1 = (ss.split("/")[1]).replace("/","")
  tag2 = (ss.split("/")[2]).replace("/","")
  tag = tag1+"_"+tag2
  pset = "flat-data-cfg.py"
  name = "crab_"+tag+".py"
  print "Creating file: "+name+", cfg file: "+pset

  crabSubmitFile.write("rm -rf crab_"+tag+"\n")
  crabSubmitFile.write("crab submit "+name+"\n")

  file = open(name,"w")
  file.write("from CRABClient.UserUtilities import config, getUsernameFromSiteDB\n")
  file.write("\n")
  file.write("config = config()\n")
  file.write("config.General.requestName = \'"+tag+"\'\n")
  file.write("config.General.transferOutputs = True\n")
  file.write("config.General.transferLogs = False\n")
  file.write("config.JobType.pluginName = \'Analysis\'\n")
  file.write("config.JobType.psetName = \'"+pset+"\'\n")
  file.write("config.JobType.maxJobRuntimeMin = 2750\n")
  file.write("config.JobType.allowUndistributedCMSSW = True\n")
  file.write("config.Data.inputDataset = \'"+ss+"\'\n")
  file.write("config.Data.inputDBS = \'global\'\n")
  file.write("config.Data.splitting = \'LumiBased\'\n")
  file.write("config.Data.unitsPerJob = 50\n")
  file.write("config.Data.lumiMask = \'Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON.txt\'\n")
  file.write("config.Data.outLFNDirBase = \'/store/group/cmst3/user/kkousour/ttbar/\'\n")
  file.write("config.Data.publication = False\n")
  file.write("config.Site.storageSite = \'T2_CH_CERN\'\n")
  file.close()

crabSubmitFile.close()
