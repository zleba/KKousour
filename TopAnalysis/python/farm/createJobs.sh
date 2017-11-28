ADDRESS=/afs/desy.de/user/z/zlebcr/cms/CMSSW_8_0_20/src/KKousour/TopAnalysis/python/farm
DUST=/nfs/dust/cms/user/zlebcr/JEC/histos

step=3
nFiles=`cat runsG | wc -l`

rm sub/*.sub

for i in `seq  1 $step $nFiles`
do
    name=jetCorr$i
    subFile=sub/${name}.sub
    echo "#!/bin/bash" > $subFile
    echo "#$ -V" >> $subFile
    echo "#$ -e $ADDRESS/err/${name}.err" >> $subFile
    echo "#$ -o $ADDRESS/out/${name}.out" >> $subFile

    echo 'cd $TMPDIR' >> $subFile
    echo 'pwd' >> $subFile
    echo "cp $ADDRESS/flatData-TTJets-cfg.py ." >> $subFile
    echo "cmsRun flatData-TTJets-cfg.py startFile=$i  nFiles=$step outputFile=$DUST/${name}.root" >> $subFile 
    #echo "cp *.root $DUST/${name}.root" >> $subFile 

done
