#include "FillHistograms.C"

void FillHistogramsAll()
{
  
  FillHistograms("TT_TuneCUETP8M1_13TeV-powheg-pythia8",true,true,true);
 
  FillHistograms("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true,false,false);
  FillHistograms("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true,false,false);
  FillHistograms("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true,false,false);
  FillHistograms("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true,false,false);
  FillHistograms("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true,false,false);
  FillHistograms("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true,false,false);
  
  FillHistograms("JetHT",false,false,false);
  
  FillHistograms("DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",true,false,false);
  FillHistograms("WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",true,false,false);
  
  FillHistograms("ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",true,false,false);
  FillHistograms("ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",true,false,false);
  FillHistograms("ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",true,false,false);
  FillHistograms("ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",true,false,false);
 
  /*
  FillHistograms("TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",true,false); 
  FillHistograms("TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true,false);
  FillHistograms("TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8",true,false);
  FillHistograms("TT_TuneCUETP8M1_13TeV-powheg-scaleup-pythia8",true,false);
  FillHistograms("TT_TuneCUETP8M1_mtop1665_13TeV-powheg-pythia8",true,false);
  FillHistograms("TT_TuneCUETP8M1_mtop1695_13TeV-powheg-pythia8",true,false);
  FillHistograms("TT_TuneCUETP8M1_mtop1715_13TeV-powheg-pythia8",true,false);
  FillHistograms("TT_TuneCUETP8M1_mtop1735_13TeV-powheg-pythia8",true,false);
  FillHistograms("TT_TuneCUETP8M1_mtop1755_13TeV-powheg-pythia8",true,false);
  FillHistograms("TT_TuneCUETP8M1_mtop1785_13TeV-powheg-pythia8",true,false);
  FillHistograms("TT_TuneCUETP8M1mpiOFF_13TeV-powheg-pythia8",true,false);
  FillHistograms("TT_TuneEE5C_13TeV-powheg-herwigpp",true,false);
  
  //FillHistograms("QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",true,false);
  
  */
 
  
}
