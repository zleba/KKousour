#! /usr/bin/env python
import os, csv, sys

DIR = ['MultiJetA','BJetPlusXB','BJetPlusXC','BJetPlusXD','VBF1ParkedB','VBF1ParkedC']

counter = 0
for d in DIR:
  print 'Creating csv files from directory: log'+d
  files = os.listdir('log'+d)
  for ff in files:
    if (ff.find('Job_') == 0):
      csv_name = 'json.csv'
      command = 'pixelLumiCalc.py -i log'+d+'/'+ff+'/json.txt -o log'+d+'/'+ff+'/'+csv_name+' --nowarning overview'
      os.system(command)
  counter += 1
