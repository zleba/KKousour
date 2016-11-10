#include "FillHistograms.C"
void FillHistogramsAll()
{
  FillHistograms("ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8",true);
  FillHistograms("ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix",true);

  FillHistograms("TT_TuneCUETP8M1_13TeV-powheg-pythia8",true);
  
  FillHistograms("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true);
  FillHistograms("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true);
  FillHistograms("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true);
  FillHistograms("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true);
  FillHistograms("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true);
  FillHistograms("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true);
  
  FillHistograms("DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",true);
  FillHistograms("WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",true);
  FillHistograms("ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",true);
  FillHistograms("ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",true);
  
  FillHistograms("WWTo4Q_13TeV-powheg",true);
  FillHistograms("ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8",true);
  
  FillHistograms("WZ_TuneCUETP8M1_13TeV-pythia8",true);
  FillHistograms("ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",true);
  FillHistograms("TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8",true);
  FillHistograms("TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8",true);
  
  FillHistograms("JetHT",false);
}
