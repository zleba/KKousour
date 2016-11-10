#include "DrawVariable.C"

void DrawVariablesAll(TString ss,bool PRINT)
{
  //float LUMI = 513;//control path: 0 btag
  float LUMI = 15941;//signal path: 1,2 btag
   
  DrawVariable("boosted","mTop_"+ss    ,LUMI,false,10,50,300,"Top mass (GeV)",false,510,PRINT);
  DrawVariable("boosted","mW_"+ss    ,LUMI,false,5,10,200,"W mass (GeV)",false,510,PRINT);
  
  //DrawVariable("boosted","nvtx_"+ss                ,LUMI,false,1,0,35,"Number of vertices",false,510,PRINT);
  //DrawVariable("boosted","ht_"+ss                  ,LUMI,true,5,600,3000,"H_{T} (GeV)",false,510,PRINT);
  //DrawVariable("boosted","nJets_"+ss               ,LUMI,true,1,2,5,"Number of jets",true,105,PRINT);
  //DrawVariable("boosted","mJJ_"+ss                 ,LUMI,true,100,400,4000,"mJJ (GeV)",false,510,PRINT);
  //DrawVariable("boosted","ptJJ_"+ss                ,LUMI,true,20,0,1000,"p_{T}JJ (GeV)",false,510,PRINT);
  //DrawVariable("boosted","yJJ_"+ss                 ,LUMI,true,10,-3,3,"yJJ",false,510,PRINT);
  //DrawVariable("boosted","dPhiJJ_"+ss              ,LUMI,true,10,0,3.142,"d#phiJJ",false,510,PRINT);
  //DrawVariable("boosted","nBJets_"+ss              ,LUMI,true,1,0,5,"Number of b jets",true,105,PRINT);
  //DrawVariable("boosted","mva_"+ss                 ,LUMI,true,5,-1.5,1.5,"Fisher discriminant",false,510,PRINT);
  //DrawVariable("boosted","jetTau32_"+ss            ,LUMI,true,2,0,1,"Leading jet #tau3/#tau2",false,510,PRINT);
  //DrawVariable("boosted","jetTau31_"+ss            ,LUMI,true,2,0,1,"Leading jet #tau3/#tau1",false,510,PRINT);
  //DrawVariable("boosted","jetMassSoftDrop1_"+ss    ,LUMI,true,10,50,300,"Second jet softDrop mass (GeV)",false,510,PRINT);
  //DrawVariable("boosted","jetMassSub0_"+ss         ,LUMI,true,5,0,200,"Leading jet leading subjet mass (GeV)",false,510,PRINT);
  //DrawVariable("boosted","jetPt_"+ss               ,LUMI,true,20,350,1200,"Leading jet p_{T} (GeV)",false,510,PRINT);
  //DrawVariable("boosted","jetEta_"+ss              ,LUMI,true,2,-2.5,2.5,"Leading jet #eta",false,510,PRINT);
 
}






