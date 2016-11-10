#include "DrawQCDClosureBoosted.C"
void DrawQCDClosureBoostedAll(TString CUT)
{
  //DrawQCDClosureBoosted("mTop",CUT,"Jet SoftDrop mass (GeV)",10,70,300,"pol1",false);
  DrawQCDClosureBoosted("mW",CUT,"Subjet mass (GeV)",2,0,120,"pol1",false);
  /*
  DrawQCDClosureBoosted("jetPt",CUT,"Jet p_{T} (GeV)",25,350,1200,"pol1",false);
  DrawQCDClosureBoosted("mJJ",CUT,"m_{t#bar{t}}",50,800,2500,"pol1",false);  
  DrawQCDClosureBoosted("yJJ",CUT,"y_{t#bar{t}}",25,-2.5,2.5,"pol2",false);
  DrawQCDClosureBoosted("ptJJ",CUT,"p_{T,t#bar{t}}",10,0,500,"pol1",false);
  DrawQCDClosureBoosted("mva",CUT,"",5,-1,1,"pol1",false);
  */
}
