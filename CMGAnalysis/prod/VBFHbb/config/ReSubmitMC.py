#! /usr/bin/env python
import os, csv, sys

DIR = ['QCD100','QCD250','QCD500','QCD1000','ZJets','WJets','TTJets','T_s-channel','T_t-channel','T_tW-channel','Tbar_s-channel','Tbar_t-channel','Tbar_tW-channel','VBFPowheg125']

TYPE = 'log'
SUBMIT = False

for d in DIR:
  print 'Re-submitting failed jobs in directory: '+TYPE+d
  files = os.listdir(TYPE+d)
  for ff in files:
    if (ff.find('Job_') == 0):
      local_files = os.listdir(TYPE+d+'/'+ff)
      if ('json.txt' not in local_files):
        command = 'bsub -q 8nh < '+TYPE+d+'/'+ff+'/batchScript.sh'
        print command
        if (SUBMIT):
          os.system(command)
