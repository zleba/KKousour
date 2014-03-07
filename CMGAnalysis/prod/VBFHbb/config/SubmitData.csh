#! /bin/csh

foreach SAMPLE ('MultiJetA' 'BJetPlusXB' 'BJetPlusXC' 'BJetPlusXD' 'VBF1ParkedB' 'VBF1ParkedC' 'VBF1ParkedD')
#foreach SAMPLE ('JetA' 'JetMonB' 'JetMonC' 'JetMonD')
	echo $SAMPLE
        echo 'cleaning' 
        eos rm -r /eos/cms/store/cmst3/group/vbfhbb/flat/$SAMPLE
        eos rm -r /eos/cms/store/cmst3/group/vbfhbb/flat/$SAMPLE
        rm -r log$SAMPLE 
        echo 'submitting'
        cmsBatch.py 20 flat-$SAMPLE-cfg.py -o log$SAMPLE -r /store/cmst3/group/vbfhbb/flat/$SAMPLE -b 'bsub -q 8nh < ./batchScript.sh'
end
