#! /usr/bin/env python
import os
from CMGTools.Production.eostools import *

path = '/store/cmst3/group/vbfhbb/flat/'

datasets = ['MultiJetA','BJetPlusXB','BJetPlusXC','BJetPlusXD','VBF1ParkedB','VBF1ParkedC','VBF1ParkedD']
#datasets = ['BJetPlusXD','VBF1ParkedB','VBF1ParkedC','VBF1ParkedD']
for ds in datasets:
  print 'Adding ROOT files from dataset: /'+ds
  files = ls_EOS(path+'/'+ds)
  command = 'hadd -f ./flatTree_'+ds+'.root '

  for ff in files:
    if (ff.find('.root') > 0) :
      pfn = lfnToPFN(ff)
      print '\"'+pfn+'\",'
      command += pfn + ' '

  os.system(command)
