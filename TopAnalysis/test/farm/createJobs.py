#!/usr/local/bin/python2.7

import os
curDir = os.path.dirname(os.path.realpath(__file__))
#listADDRESS = curDir + '/fileLists/Aug17'
listADDRESS = curDir + '/fileLists/QCD_TuneCUETP8M1_13TeV_pythia8'

DUST = '/nfs/dust/cms/user/zlebcr/JEC/ntuplesJCunc'

step=6

Tags = filter(lambda x: x.endswith('.txt'), os.listdir(listADDRESS))
Tags = [s.replace('.txt','') for s in Tags]
print Tags

#import sys
#sys.exit(0)

for run in Tags:
    File= run + '.txt'
    nFiles = sum(1 for line in open(listADDRESS +'/' + File ))

    for i in range( (nFiles-1) / step + 1):
        #print nFiles, i+1, step*i+1
        name='jets{}_{:0>3}'.format(run, i+1)
        subFile='sub/' + name + '.sub'

        fSub = open(subFile, "w")
        fSub.writelines('\n'.join([       
        "#!/bin/bash",
        "#$ -V",
        "#$ -e {}/err/{}.err".format(curDir, name),
        "#$ -o {}/out/{}.out".format(curDir, name),
        "#$ -l h_vmem=2G",
        "pwd",
        "cmsRun "+curDir+"/flatData-new.py  listFile="+listADDRESS+"/"+File+
        " startFile="+str(i*step+1)+"  nFiles="+str(step)+
        " outputFile="+ DUST+"/"+name+".root",
        ]))
        fSub.close()

        import stat
        os.chmod(subFile, stat.S_IRWXU + stat.S_IRGRP)


