#! /usr/bin/env python
import os, csv, sys

DIR = ['MultiJetA','BJetPlusXB','BJetPlusXC','BJetPlusXD','VBF1ParkedB','VBF1ParkedC']

TYPE = 'log'
SUBMIT = False
PWD = os.path.abspath(".")+'/'
for d in DIR:
  print 'Re-submitting failed jobs in directory: '+TYPE+d
  files = os.listdir(TYPE+d)
  for ff in files:
    if (ff.find('Job_') == 0):
      local_files = os.listdir(TYPE+d+'/'+ff)
      if ('json.txt' not in local_files):
        command = 'cd '+PWD+TYPE+d+'/'+ff
        print command
        os.system(command)
        command = 'bsub -q 8nh < batchScript.sh'
        print command
        if (SUBMIT):
          os.system(command)
        os.system('cd '+PWD) 
        
