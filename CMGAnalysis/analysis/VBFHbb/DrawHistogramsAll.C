#include "DrawHistogram.C"
void DrawHistogramsAll(TString SELECTION,bool PRINT)
{ 
  //---- event variables -----------------------------------------------------
  
  DrawHistogram("h_nVtx"      ,SELECTION,"Number of vertices"    ,0   ,50  ,1,true,PRINT,true,510,"");
  DrawHistogram("h_dEtaqq"    ,SELECTION,"|#Delta#eta_{qq}|"     ,2.5 ,8.0 ,2,true,PRINT,true,510,"");
  DrawHistogram("h_softHt"    ,SELECTION,"Soft H_{T} (GeV)"      ,0   ,200 ,4,true,PRINT,true,510,"");
  
  DrawHistogram("h_nSoftJets" ,SELECTION,"Soft Multiplicity"     ,0   ,100 ,1,true,PRINT,true,510,"");
  DrawHistogram("h_nSoftJets2",SELECTION,"Soft Multiplicity"     ,0   ,100 ,1,true,PRINT,true,510,"");
  //DrawHistogram("h_nSoftJets5",SELECTION,"Soft Multiplicity"     ,0   ,100 ,1,true,PRINT,true,510,"");
  
  DrawHistogram("h_cosTheta"  ,SELECTION,"cos#theta"             ,-1  ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_mqq"       ,SELECTION,"M_{qq} (GeV)"          ,250   ,2500,10,true,PRINT,true,510,"");
  //DrawHistogram("h_dPhibb"    ,SELECTION,"#Delta#phi_{bb}"       ,0   ,1.99,2,true,PRINT,true,510,"");
  DrawHistogram("h_dPhiqq"    ,SELECTION,"#Delta#phi_{qq}"       ,0   ,3.15,5,true,PRINT,true,510,"");
  DrawHistogram("h_BDT"       ,SELECTION,"BDT Output"            ,1   ,1   ,1,true,PRINT,true,510,"");
  
  DrawHistogram("h_mbbReg"    ,SELECTION,"Regressed M_{bb} (GeV)",30  ,250 ,2,true,PRINT,true,510,"");
  //DrawHistogram("h_mbbRegCut" ,SELECTION,"Regressed M_{bb} (GeV)",30  ,250 ,5,true,PRINT,true,510,"");
  DrawHistogram("h_x1_ov_ht"  ,SELECTION,"log(x_{1}/H_{T})"      ,-12 ,-5  ,2,true,PRINT,true,510,"");
  DrawHistogram("h_x2_ov_ht"  ,SELECTION,"log(x_{2}/H_{T})"      ,-12 ,-5  ,2,true,PRINT,true,510,"");
  //DrawHistogram("h_ht"        ,SELECTION,"H_{T} (GeV)"           ,0 ,2000  ,5,true,PRINT,true,510,"");
  //---- jet variables --------------------------------------------------------
  /*
  DrawHistogram("h_jetEta0"       ,SELECTION,"Jet0 #eta"            ,-5  ,5   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetEta1"       ,SELECTION,"Jet1 #eta"            ,-5  ,5   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetEta2"       ,SELECTION,"Jet2 #eta"            ,-5  ,5   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetEta3"       ,SELECTION,"Jet3 #eta"            ,-5  ,5   ,1,true,PRINT,true,510,""); 

  DrawHistogram("h_jetEtaBtag0"       ,SELECTION,"BJet0 #eta"       ,-5  ,5   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetEtaBtag1"       ,SELECTION,"BJet1 #eta"       ,-5  ,5   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetEtaBtag2"       ,SELECTION,"BJet2 #eta"       ,-5  ,5   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetEtaBtag3"       ,SELECTION,"BJet3 #eta"       ,-5  ,5   ,1,true,PRINT,true,510,"");

  DrawHistogram("h_jetPtBtag0"        ,SELECTION,"BJet0 p_{T} (GeV)",40  ,1000,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetPtBtag1"        ,SELECTION,"BJet1 p_{T} (GeV)",40  ,1000,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetPtBtag2"        ,SELECTION,"BJet2 p_{T} (GeV)",40  ,1000,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetPtBtag3"        ,SELECTION,"BJet3 p_{T} (GeV)",40  ,1000,2,true,PRINT,true,510,"");

  DrawHistogram("h_jetPt0"            ,SELECTION,"Jet0 p_{T} (GeV)" ,85  ,1000,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetPt1"            ,SELECTION,"Jet1 p_{T} (GeV)" ,70  ,1000,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetPt2"            ,SELECTION,"Jet2 p_{T} (GeV)" ,60  ,500 ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetPt3"            ,SELECTION,"Jet3 p_{T} (GeV)" ,40  ,350 ,1,true,PRINT,true,510,"");
  
  DrawHistogram("h_jetBtag0"          ,SELECTION,"BJet0 CSV"        ,0   ,1   ,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetBtag1"          ,SELECTION,"BJet1 CSV"        ,0   ,1   ,2,true,PRINT,true,510,"");
 
  DrawHistogram("h_jetQGL0"           ,SELECTION,"BJet1 QGL"        ,0   ,1   ,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetQGL1"           ,SELECTION,"BJet2 QGL"        ,0   ,1   ,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetQGL2"           ,SELECTION,"BJet2 QGL"        ,0   ,1   ,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetQGL3"           ,SELECTION,"BJet3 QGL"        ,0   ,1   ,2,true,PRINT,true,510,"");
 
  DrawHistogram("h_jetPuMva0"         ,SELECTION,"BJet0 puMVA"      ,-1  ,1   ,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetPuMva1"         ,SELECTION,"BJet1 puMVA"      ,-1  ,1   ,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetPuMva2"         ,SELECTION,"BJet2 puMVA"      ,-1  ,1   ,2,true,PRINT,true,510,"");
  DrawHistogram("h_jetPuMva3"         ,SELECTION,"BJet3 puMVA"      ,-1  ,1   ,2,true,PRINT,true,510,"");
  
  DrawHistogram("h_jetChf0"           ,SELECTION,"Jet0 Chf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetChf1"           ,SELECTION,"Jet1 Chf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetChf2"           ,SELECTION,"Jet2 Chf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetChf3"           ,SELECTION,"Jet3 Chf"         ,0   ,1   ,1,true,PRINT,true,510,"");

  DrawHistogram("h_jetNhf0"           ,SELECTION,"Jet0 Nhf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetNhf1"           ,SELECTION,"Jet1 Nhf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetNhf2"           ,SELECTION,"Jet2 Nhf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetNhf3"           ,SELECTION,"Jet3 Nhf"         ,0   ,1   ,1,true,PRINT,true,510,"");

  DrawHistogram("h_jetPhf0"           ,SELECTION,"Jet0 Phf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetPhf1"           ,SELECTION,"Jet1 Phf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetPhf2"           ,SELECTION,"Jet2 Phf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetPhf3"           ,SELECTION,"Jet3 Phf"         ,0   ,1   ,1,true,PRINT,true,510,"");
 
  DrawHistogram("h_jetMuf0"           ,SELECTION,"Jet0 Muf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetMuf1"           ,SELECTION,"Jet1 Muf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetMuf2"           ,SELECTION,"Jet2 Muf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetMuf3"           ,SELECTION,"Jet3 Muf"         ,0   ,1   ,1,true,PRINT,true,510,"");

  DrawHistogram("h_jetElf0"           ,SELECTION,"Jet0 Elf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetElf1"           ,SELECTION,"Jet1 Elf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetElf2"           ,SELECTION,"Jet2 Elf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  DrawHistogram("h_jetElf3"           ,SELECTION,"Jet3 Elf"         ,0   ,1   ,1,true,PRINT,true,510,"");
  */
}
