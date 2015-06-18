#include "DrawVariable.C"

void DrawVariablesAll()
{
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","mva"   ,true,25,-1,1,"BDT Output");

  //DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","nJets"   ,true,6,6,12,"Number of jets");
  //DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","nBJets"  ,true,8,2,10,"Number of b jets");
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","ht"      ,true,45,500,1400,"H_{T} (GeV)");
  //DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","mbbMin"  ,true,40,0,400,"m_{bb} (min #DeltaR pair) (GeV)");
  //DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","dRbbMin" ,true,60,0,6,"Minimum #DeltaR_{bb}");
  //DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","sphericity",true,50,0,1,"Sphericity");
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","aplanarity",true,50,0,0.5,"Aplanarity");
  //DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","foxWolfram[0]",true,50,0.2,0.45,"H0");
  //DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","foxWolfram[1]",true,50,-0.2,0.15,"H1");
  //DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","foxWolfram[2]",true,50,-0.05,0.4,"H2");
  //DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","foxWolfram[3]",true,50,-0.2,0.2,"H3");
  /*
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","qglAve",true,50,0,1.001,"qglAve");
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","qglMin",true,50,0,1.001,"qglMin");
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","qglMedian",true,50,0,1.001,"qglMedian");

  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","jetPt[0]",true,80,0,800,"JetPt0 (GeV)");
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","jetPt[1]",true,100,0,500,"JetPt1 (GeV)");
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","jetPt[2]",true,60,0,300,"JetPt2 (GeV)");
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","jetPt[3]",true,40,0,200,"JetPt3 (GeV)");
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","jetPt[4]",true,32,0,160,"JetPt4 (GeV)");
  DrawVariable("nBJets>1 && ht>500 && jetPt[5]>40","jetPt[5]",true,24,0,120,"JetPt5 (GeV)");
  */
}
