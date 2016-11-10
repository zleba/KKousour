#! /usr/bin/env python
import os
from KKousour.TopAnalysis.eostools import *

path = '/store/cmst3/user/kkousour/ttV/'
 
sample = [
  "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8",
  "TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8",
  "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8",
  "TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8",
  "DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", 
  "WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "TT_TuneCUETP8M1_13TeV-powheg-pythia8",
  "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
  "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
  "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
  "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
  "WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8",
  "WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8",
  "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",
  "WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
  "WZ_TuneCUETP8M1_13TeV-pythia8",
  "WWTo2L2Nu_13TeV-powheg",
  "WWToLNuQQ_13TeV-powheg",
  "WWTo4Q_13TeV-powheg", 
  "ZZTo4L_13TeV_powheg_pythia8",
  "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8",
  "ZZTo2L2Nu_13TeV_powheg_pythia8",
  "ZZTo2Q2Nu_13TeV_amcatnloFXFX_madspin_pythia8",
  "ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8"
]

for ss in sample:
  dirs = ls_EOS(path+ss+"/crab_"+ss)
  timestamps = []
  if len(dirs) == 0: continue
  for dd in dirs:
    timestamps.append(dd.split('/')[-1])     
  timestamps.sort(reverse=True)
  #--- delete older folders -------------
  print "Removing old directories: "
  print timestamps[1:] 
  for ii in range(1,len(timestamps)):
    location = path+ss+"/crab_"+ss+"/"+timestamps[ii]
    command = "/afs/cern.ch/project/eos/installation/cms/bin/eos.select rm -r "+location
    print command
    os.system(command)
  #--- find the folders ------------------------------------
  dirs = ls_EOS(path+ss+"/crab_"+ss+"/"+timestamps[0])
  print dirs
  if len(dirs) == 0: continue
  for dd in dirs:
    files = ls_EOS(dd)
    nfiles = len(files)
    print 'Adding ROOT files from location: '+dd
    tag = str(dd.split('/')[-1])
    command1 = 'hadd -f ./flatTree_'+ss+'_'+tag+'_1.root '
    command2 = 'hadd -f ./flatTree_'+ss+'_'+tag+'_2.root '
    for ff in files:
      if (ff.find('.root') > 0):
        nn = int((ff.split('_')[-1]).split('.')[0])
        pfn = lfnToPFN(ff)
        print '\"'+pfn+'\",'
        if ((nn % 1000) < 500):
          command1 += pfn + ' '
        else:
          command2 += pfn + ' '
    #if (tag != "0000"):
    os.system(command1)
    #os.system(command2)




