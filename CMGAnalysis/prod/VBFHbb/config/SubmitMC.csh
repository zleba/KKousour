#! /bin/csh

foreach SAMPLE ('QCD100' 'QCD250' 'QCD500' 'QCD1000' 'ZJets' 'WJets' 'TTJets' 'T_s-channel' 'T_t-channel' 'T_tW-channel' 'Tbar_s-channel' 'Tbar_t-channel' 'Tbar_tW-channel' 'VBFPowheg115' 'VBFPowheg120' 'VBFPowheg125' 'VBFPowheg130' 'VBFPowheg135' 'GFPowheg115' 'GFPowheg120' 'GFPowheg125' 'GFPowheg130' 'GFPowheg135' 'VBFPowheg125Pythia8' 'GFPowheg125Herwig' 'GFPowheg125Pythia8' 'GFMadgraph125')
#foreach SAMPLE ('VBFPowheg125Pythia8' 'GFPowheg125Herwig' 'GFPowheg125Pythia8' 'GFMadgraph125')
#foreach SAMPLE ('VBFPowheg115-NewGenjets' 'VBFPowheg120-NewGenjets' 'VBFPowheg125-NewGenjets' 'VBFPowheg130-NewGenjets' 'VBFPowheg135-NewGenjets')
	echo $SAMPLE
        echo 'cleaning' 
        eos rm -r /eos/cms/store/cmst3/group/vbfhbb/flat/$SAMPLE
        eos rm -r /eos/cms/store/cmst3/group/vbfhbb/flat/$SAMPLE
        rm -r log$SAMPLE 
        echo 'submitting'
        cmsBatch.py 20 flat-$SAMPLE-cfg.py -o log$SAMPLE -r /store/cmst3/group/vbfhbb/flat/$SAMPLE -b 'bsub -q 8nh < ./batchScript.sh'
end
