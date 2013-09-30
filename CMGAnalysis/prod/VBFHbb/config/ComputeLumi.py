#! /usr/bin/env python
import os, csv, sys

DIR = ['MultiJetA','BJetPlusXB','BJetPlusXC','BJetPlusXD','VBF1ParkedB','VBF1ParkedC']

counter = 0
for d in DIR:
  print 'Reading csv files from directory: log'+d
  total_lumi = 0
  subdir = os.listdir('log'+d)
  for ss in subdir:
    if (ss.find('Job_') == 0):
      files = os.listdir('log'+d+'/'+ss)
      csv_name = 'json.csv'
      if (csv_name in files):  
        lumi_list_csv = open('log'+d+'/'+ss+'/'+csv_name,'rb')
        for i in range(1):
          lumi_list_csv.next()
        lumi_dict = csv.DictReader(lumi_list_csv,delimiter=',',fieldnames=['Run','DeliveredLS','Delivered','SelectedLS','Recorded'])
        lumi = {}
        for l in lumi_dict:
          if (l['Recorded'] != 'n/a'):
            total_lumi+=float(l['Recorded'])      
  print 'Total Lumi = %1.3f /pb' % (total_lumi/1000000)
  counter += 1
