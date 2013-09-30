#! /usr/bin/env python
import os
from CMGTools.Production.eostools import *

path = '/store/cmst3/group/vbfhbb/flat/'

datasets = ['QCD100','QCD250','QCD500','QCD1000','ZJets','WJets',
            'TTJets','T_s-channel','T_t-channel','T_tW-channel','Tbar_s-channel','Tbar_t-channel','Tbar_tW-channel',
            'VBFPowheg115','VBFPowheg120','VBFPowheg125','VBFPowheg130','VBFPowheg135',
            'GFPowheg115','GFPowheg120','GFPowheg125','GFPowheg130','GFPowheg135']

for ds in datasets:
  print 'Adding ROOT files from dataset: '+ds
  files = ls_EOS(path+ds)
  command = 'hadd -f ./flatTree_'+ds+'.root '

  for ff in files:
    if (ff.find('.root') > 0) :
      pfn = lfnToPFN(ff)
      print '\"'+pfn+'\",'
      command += pfn + ' '

  os.system(command)
