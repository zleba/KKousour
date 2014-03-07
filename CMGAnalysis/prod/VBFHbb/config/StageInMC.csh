#! /bin/csh

foreach SAMPLE ('QCD100' 'QCD250' 'QCD500' 'QCD1000' 'ZJets' 'WJets' 'TTJets' 'T_s-channel' 'T_t-channel' 'T_tW-channel' 'Tbar_s-channel' 'Tbar_t-channel' 'Tbar_tW-channel' 'VBFPowheg115' 'VBFPowheg120' 'VBFPowheg125' 'VBFPowheg130' 'VBFPowheg135' 'GFPowheg115' 'GFPowheg120' 'GFPowheg125' 'GFPowheg130' 'GFPowheg135' 'VBFPowheg125Pythia8' 'GFPowheg125Herwig' 'GFPowheg125Pythia8' 'GFMadgraph125')
	echo 'staging in' $SAMPLE 
        cmsStage -f flatTree_$SAMPLE.root /store/cmst3/group/vbfhbb/flat/flatTree_$SAMPLE.root
end
