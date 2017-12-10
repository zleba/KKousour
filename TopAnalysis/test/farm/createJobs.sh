ADDRESS=/afs/desy.de/user/z/zlebcr/cms/CMSSW_8_0_20/src/KKousour/TopAnalysis/test/farm
listADDRESS=$ADDRESS/fileLists
DUST=/nfs/dust/cms/user/zlebcr/JEC/histos

step=6

#rm sub/*.sub

for run in B C D E F G H
do
    file=run${run}.txt
    nFiles=`cat $listADDRESS/$file | wc -l`

    for i in `seq  1 $step $nFiles`
    do
        name=jets$run$i
        subFile=sub/${name}.sub
        echo "#!/bin/bash" > $subFile
        echo "#$ -V" >> $subFile
        echo "#$ -e $ADDRESS/err/${name}.err" >> $subFile
        echo "#$ -o $ADDRESS/out/${name}.out" >> $subFile
        echo "#$ -l h_vmem=2G" >> $subFile

        #echo 'cd $TMPDIR' >> $subFile
        echo 'pwd' >> $subFile
        #echo "cp $ADDRESS/flatData-TTJets-cfg.py ." >> $subFile
        echo "cmsRun $ADDRESS/flatData-TTJets-cfg.py  listFile=$listADDRESS/$file  startFile=$i  nFiles=$step outputFile=$DUST/${name}.root" >> $subFile 
        #echo "cp *.root $DUST/${name}.root" >> $subFile 

        chmod u+x $subFile

    done

done
