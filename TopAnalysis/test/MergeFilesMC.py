#! /usr/bin/env python
import os
from eostools import *

path = '/store/cmst3/user/kkousour/ttH/'
 
sample = ['TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola',
          'TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola',
          'QCD_HT_250To500_13TeV-madgraph',
          'QCD_HT-500To1000_13TeV-madgraph',
          'QCD_HT_1000ToInf_13TeV-madgraph' 
         ]
alias = ['TTH','TTJets','QCD250','QCD500','QCD1000']

#sample = ['TTbarH_M-125_13TeV_amcatnlo-pythia8-tauola']
#alias = ['TTH']

k = 0
for ss in sample:
  print 'Sample = '+ss
  dirs = ls_EOS(path+ss+"/crab_"+alias[k])
  timestamps = []
  for dd in dirs:
    timestamps.append(dd.split('/')[-1])     
  timestamps.sort(reverse=True)
  #--- delete older folders -------------
  print "Removing old directories: "
  print timestamps[1:] 
  for ii in range(1,len(timestamps)):
    location = path+ss+"/crab_"+alias[k]+"/"+timestamps[ii]
    command = "cmsRm -r "+location
    print command
    os.system(command)
  #--- merge ROOT files from the newest folder ------------- 
  location = path+ss+"/crab_"+alias[k]+"/"+timestamps[0]+"/0000/"
  print 'Adding ROOT files from location: '+location
  files = ls_EOS(location)
  command = 'hadd -f ./flatTree_'+alias[k]+'.root '
  for ff in files:
    if (ff.find('.root') > 0):
      nn = int((ff.split('_')[-1]).split('.')[0])
      pfn = lfnToPFN(ff)
      print '\"'+pfn+'\",'
      command += pfn + ' '
  os.system(command)
  k += 1
