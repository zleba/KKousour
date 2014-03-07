#! /bin/csh

#foreach SAMPLE ('JetA' 'JetMonB' 'JetMonC' 'JetMonD' 'MultiJetA' 'BJetPlusXB' 'BJetPlusXC' 'BJetPlusXD' 'VBF1ParkedB' 'VBF1ParkedC' 'VBF1ParkedD')
foreach SAMPLE ('MultiJetA' 'BJetPlusXB' 'BJetPlusXC' 'BJetPlusXD' 'VBF1ParkedB' 'VBF1ParkedC' 'VBF1ParkedD')
	echo 'stage in' $SAMPLE 
        cmsStage -f flatTree_$SAMPLE.root /store/cmst3/group/vbfhbb/flat/flatTree_$SAMPLE.root
end
