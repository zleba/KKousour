#! /usr/bin/env python
import os, csv, sys
from CMGTools.Production.eostools import *

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("--submit",action="store_true",default=False,dest="submit")
(options, args) = parser.parse_args()
SUBMIT = options.submit
    
path = '/store/cmst3/group/vbfhbb/flat/'
#DIR = ['MultiJetA','BJetPlusXB','BJetPlusXC','BJetPlusXD','VBF1ParkedB','VBF1ParkedC','VBF1ParkedD']
DIR = [ 'MultiJetA','BJetPlusXB','BJetPlusXC','BJetPlusXD','VBF1ParkedB','VBF1ParkedC','VBF1ParkedD', 
        'QCD100','QCD250','QCD500','QCD1000','ZJets','WJets',
        'TTJets','T_s-channel','T_t-channel','T_tW-channel','Tbar_s-channel','Tbar_t-channel','Tbar_tW-channel',
        'VBFPowheg115','VBFPowheg120','VBFPowheg125','VBFPowheg130','VBFPowheg135',
        'GFPowheg115','GFPowheg120','GFPowheg125','GFPowheg130','GFPowheg135',
        'VBFPowheg125Pythia8','GFPowheg125Herwig','GFPowheg125Pythia8','GFMadgraph125'
      ]
TYPE = 'log'

PWD = os.path.abspath(".")+'/'
for d in DIR:
  print 'Checking directory: '+TYPE+d
  jobs = os.listdir(TYPE+d)
  jobNum = []
  for ff in jobs:
    if (ff.find('Job_') == 0):
      jj = str(ff.split('_')[1])
      jobNum.append(jj)
  #print "Requested jobs:"
  #print jobNum

  ntuples = ls_EOS(path+d)
  missing = []
  found = []
  
  for nn in ntuples:
    if (nn.find('flatTree') > 0):
      jj = str((nn.split('flatTree_')[1]).split('.')[0])
      found.append(jj)

  for jj in jobNum:
    if jj not in found:
      missing.append(jj)    

  print "Missing jobs:"
  print missing
  if (SUBMIT): 
    print "ReSubmitting failed jobs"
    for mm in missing:
      command = 'cd '+PWD+TYPE+d+'/Job_'+mm
      print command
      os.chdir(PWD+TYPE+d+'/Job_'+mm)
      command = 'bsub -q 8nh < batchScript.sh'
      print command
      os.system(command)
      os.chdir(PWD) 
