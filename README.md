1) Clone the repository into CMSSW src directory
2) Compile using scram b -j8
3) run test run using cmsRun flatData-TTJets-cfg.py

To do here:
1) Add all interesting variables into this nTuple producer, i.e. AK4, AK8 CHS jets and PUPI jets + our trigger elements, one can compare to the older ntuplizer of DESY group dedicated to jet analyzes
https://github.com/pgunnell/SMPJ/tree/ntuples93X
2) Check the selection when the event is stored (if possible, i.e not heavy, even having jets tagged with HLT50 tagger)
2) Run the producer on whole 2016 MINIAOD data

To do from scratch:
1) Write the analyzer of these hopefully small nTuples, fill interesting variables relevant for JEC, try matching of various flavour types to AK4CHS jets
2) Write plotting macros to plot relevant distributions for all jet types

Good luck
