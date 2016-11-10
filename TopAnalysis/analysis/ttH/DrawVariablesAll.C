#include "DrawVariable.C"

void DrawVariablesAll(TString CAT,TString TRG,TString REG,bool PRINT)
{
  TString ss(TRG+"_"+CAT);
  if (REG != "") {
    ss += "_"+REG;
  }
  /*
  DrawVariable("nvtx_"+ss      ,false,1,0,50,"Number of vertices",PRINT);
  */  
  DrawVariable("mvaQCD_"+ss       ,false,4,-1,1,"QCD BDT Output",PRINT,false,0.5,1.0);
  DrawVariable("mvaTTbar_"+ss     ,false,4,-1,1,"TTbar BDT Output",PRINT);
  DrawVariable("nBJets_"+ss       ,false,1,2,7,"Number of b jets",PRINT);
  DrawVariable("nJets_"+ss        ,false,1,6,15,"Number of jets",PRINT);
  
  DrawVariable("mbbMin_"+ss       ,false,2,0,600,"m_{bb} (min #DeltaR pair) (GeV)",PRINT);
  DrawVariable("dRbbTop_"+ss      ,false,2,0,6,"#DeltaR_{bb} (kin. fit)",PRINT);
  DrawVariable("dRbbMin_"+ss      ,false,1,0,6,"Minimum #DeltaR_{bb}",PRINT);
  
  
  
  DrawVariable("chi2_"+ss         ,false,5,0,200,"#chi^{2}",PRINT);
  DrawVariable("mTTbar_"+ss       ,false,50,0,3000,"m_{t#bar{t}} (GeV)",PRINT);
  DrawVariable("ptTTbar_"+ss      ,false,10,0,600,"p_{T,t#bar{t}} (GeV)",PRINT);
  DrawVariable("mTop[0]_"+ss      ,false,4,100,400,"m_{t} (GeV)",PRINT);
  
  
  DrawVariable("ht_"+ss           ,false,2,400,2000,"H_{T} (GeV)",PRINT);
  
  DrawVariable("sphericity_"+ss   ,false,2,0,1,"Sphericity",PRINT);
  DrawVariable("aplanarity_"+ss   ,false,2,0,0.5,"Aplanarity",PRINT);
  DrawVariable("centrality_"+ss   ,false,2,0,1,"Centrality",PRINT);
  DrawVariable("cosThetaStar1_"+ss,false,4,-1,1,"cosThetaStar1",PRINT);
  DrawVariable("cosThetaStar2_"+ss,false,4,-1,1,"cosThetaStar2",PRINT);
  DrawVariable("EtStar1_"+ss      ,false,5,0,500,"EtStar1",PRINT);
  DrawVariable("EtStar2_"+ss      ,false,5,0,500,"EtStar2",PRINT);
  DrawVariable("FW0_"+ss          ,false,1,0.2,0.45,"H0",PRINT);
  DrawVariable("FW1_"+ss          ,false,1,-0.2,0.15,"H1",PRINT);
  DrawVariable("FW2_"+ss          ,false,1,-0.05,0.4,"H2",PRINT);
  DrawVariable("FW3_"+ss          ,false,1,-0.2,0.2,"H3",PRINT);
  
  /*  
  DrawVariable("jetPt0_"+ss       ,false,5,40,1000,"JetPt0 (GeV)",PRINT);
  DrawVariable("jetPt1_"+ss       ,false,5,40,600,"JetPt1 (GeV)",PRINT);
  DrawVariable("jetPt2_"+ss       ,false,2,40,400,"JetPt2 (GeV)",PRINT);
  DrawVariable("jetPt3_"+ss       ,false,2,40,300,"JetPt3 (GeV)",PRINT);
  DrawVariable("jetPt4_"+ss       ,false,1,40,200,"JetPt4 (GeV)",PRINT);
  DrawVariable("jetPt5_"+ss       ,false,1,40,160,"JetPt5 (GeV)",PRINT);
  */ 
  /*
  DrawVariable("jetEta0_"+ss      ,false,2,-3,3,"JetEta0",PRINT);
  DrawVariable("jetEta1_"+ss      ,false,2,-3,3,"JetEta0",PRINT);
  DrawVariable("jetEta2_"+ss      ,false,2,-3,3,"JetEta0",PRINT);
  DrawVariable("jetEta3_"+ss      ,false,2,-3,3,"JetEta0",PRINT);
  DrawVariable("jetEta4_"+ss      ,false,2,-3,3,"JetEta0",PRINT);
  DrawVariable("jetEta5_"+ss      ,false,2,-3,3,"JetEta0",PRINT);
  */
  
  /*
  DrawVariable("qglMin_"+ss    ,false,2,0,1,"JetQGLMin",PRINT);
  DrawVariable("qglMedian_"+ss ,false,2,0,1,"JetQGLMedian",PRINT);
  DrawVariable("qglAve_"+ss    ,false,2,0,1,"JetQGLAve",PRINT);
  DrawVariable("jetQGL0_"+ss   ,false,2,0,1,"JetQGL0",PRINT);
  DrawVariable("jetQGL1_"+ss   ,false,2,0,1,"JetQGL1",PRINT);
  DrawVariable("jetQGL2_"+ss   ,false,2,0,1,"JetQGL2",PRINT);
  DrawVariable("jetQGL3_"+ss   ,false,2,0,1,"JetQGL3",PRINT);
  DrawVariable("jetQGL4_"+ss   ,false,2,0,1,"JetQGL4",PRINT);
  DrawVariable("jetQGL5_"+ss   ,false,2,0,1,"JetQGL5",PRINT); 
  */
  /*
  DrawVariable("jetBtag0_"+ss   ,false,2,0,1,"JetBtag0",PRINT);
  DrawVariable("jetBtag1_"+ss   ,false,2,0,1,"JetBtag1",PRINT);
  DrawVariable("jetBtag2_"+ss   ,false,2,0,1,"JetBtag2",PRINT);
  DrawVariable("jetBtag3_"+ss   ,false,2,0,1,"JetBtag3",PRINT);
  DrawVariable("jetBtag4_"+ss   ,false,2,0,1,"JetBtag4",PRINT);
  DrawVariable("jetBtag5_"+ss   ,false,2,0,1,"JetBtag5",PRINT);
  */
  /*
  DrawVariable("jetChf0_"+ss      ,false,2,-3,3,"JetChf0",PRINT);
  DrawVariable("jetChf1_"+ss      ,false,2,-3,3,"JetChf1",PRINT);
  DrawVariable("jetChf2_"+ss      ,false,2,-3,3,"JetChf2",PRINT);
  DrawVariable("jetChf3_"+ss      ,false,2,-3,3,"JetChf3",PRINT);
  DrawVariable("jetChf4_"+ss      ,false,2,-3,3,"JetChf4",PRINT);
  DrawVariable("jetChf5_"+ss      ,false,2,-3,3,"JetChf5",PRINT);
  */
}
