#include "KKousour/CMGAnalysis/plugins/QGLCalculator.h"

QGLCalculator::QGLCalculator()
{
  std::string tmp_pt_cat[7]  = {"30to50","50to80","80to120","120to170","170to300","300to470","470to600"};
  std::string tmp_eta_cat[3] = {"central","transition","forward"};
  float  tmp_pt_bnd[8]  = {0,50,80,120,170,300,470,1000};
  float  tmp_eta_bnd[4] = {0,2,3,5};
  float  tmp_cor_rho[7][3][6] = {
  {
    {0.000696, 0.000642, 0.000060, -0.012049, -0.002013 , -0.002029}, 
    {0.001469, 0.001128, 0.000099, 0.230114, -0.003204, -0.003301},
    {0.000960, 0.000958, 0.000063, 0.228834, -0.003283, -0.003322}
  },  
  {
    {0.000446, 0.000406, 0.000033, -0.002534, -0.001154 , -0.001312},
    {0.001347, 0.000964, 0.000100, 0.196515, -0.002769, -0.002907}, 
    {0.000845, 0.000784, 0.000059, 0.211113, -0.003340, -0.003351}
  },
  { 
    {0.000291, 0.000264, 0.000018, 0.003283, -0.000713 , -0.000919}, 
    {0.001045, 0.000736, 0.000076, 0.182445, -0.002352, -0.002531}, 
    {0.000630, 0.000590, 0.000043, 0.201430, -0.003097, -0.003090} 
  },
  {
    {0.000208, 0.000184, 0.000011, 0.005118, -0.000486 , -0.000667}, 
    {0.000754, 0.000540, 0.000049, 0.170863, -0.001904, -0.002085}, 
    {0.000470, 0.000451, 0.000030, 0.196330, -0.002770, -0.002766} 
  },
  {
    {0.000146, 0.000127, 0.000007, 0.004974, -0.000337 , -0.000476}, 
    {0.000541, 0.000395, 0.000030, 0.168255, -0.001517, -0.001685}, 
    {0.000384, 0.000361, 0.000020, 0.188181, -0.002517, -0.002458} 
  },
  {
    {0.000091, 0.000070, 0.000004, 0.004071, -0.000231 , -0.000275}, 
    {0.000310, 0.000243, 0.000014, 0.170718, -0.001004, -0.001141},
    {0.000410, 0.000226, -0.000097, 0.135866, -0.003550, -0.003048} 
  },
  {
    {0.000062, 0.000043, 0.000003, 0.003028, -0.000189 , -0.000139}, 
    {0.000211, 0.000163, 0.000009, 0.182038, -0.000792, -0.000895}, 
    {-0.002664, -0.001128, -0.000142, 0.114548, 0.005115, 0.000035} 
  }
  };
  for(int ipt=0;ipt<7;ipt++) {
    PT_CAT_[ipt] = tmp_pt_cat[ipt];
    PT_BND_[ipt] = tmp_pt_bnd[ipt];
  }
  PT_BND_[7] = tmp_pt_bnd[7];
  for(int ieta=0;ieta<3;ieta++) {
    ETA_CAT_[ieta] = tmp_eta_cat[ieta];
    ETA_BND_[ieta] = tmp_eta_bnd[ieta];
  }
  ETA_BND_[3] = tmp_eta_bnd[3];
  for(int ipt=0;ipt<7;ipt++) {
    for(int ieta=0;ieta<3;ieta++) {
      for(int ivar=0;ivar<6;ivar++) {
        COR_RHO_[ipt][ieta][ivar] = tmp_cor_rho[ipt][ieta][ivar];
      }
    }
  }
  std::cout<<"Booking QGL readers"<<std::endl;
  for(int ipt=0;ipt<7;ipt++) {
    for(int ieta=0;ieta<3;ieta++) {
      reader_[ipt][ieta] = new TMVA::Reader("!Color:!Silent");
      reader_[ipt][ieta]->AddVariable("axis1"  ,&var_[0]);
      reader_[ipt][ieta]->AddVariable("axis2"  ,&var_[1]);
      reader_[ipt][ieta]->AddVariable("Mult"   ,&var_[2]);
      reader_[ipt][ieta]->AddVariable("JetR"   ,&var_[3]);
      reader_[ipt][ieta]->AddVariable("JetPull",&var_[4]);
      std::string ss = PT_CAT_[ipt]+"_"+ETA_CAT_[ieta]+"_Likelihood.xml";
      edm::FileInPath f1("KKousour/CMGAnalysis/data/"+ss);
      reader_[ipt][ieta]->BookMVA("LIK_"+PT_CAT_[ipt]+"_"+ETA_CAT_[ieta],f1.fullPath());
    }
  } 
}
//-------------------------------------------------------------
int QGLCalculator::FindIndex(int N, const float BND[], float x)
{
  int index(-1);
  for(int i=0;i<N-1;i++) {
    if ((x > BND[i]) && (x < BND[i+1])) {
      index = i;
      return index;
    }
  }
  return index;
}
//-------------------------------------------------------------
float QGLCalculator::getQGL(cmg::PFJet const& jet,float rho)
{
  int ipt  = TMath::Max(0,FindIndex(8,PT_BND_,jet.pt()));
  int ieta = TMath::Max(0,FindIndex(4,ETA_BND_,fabs(jet.eta())));
  if (ieta == 0) {
    var_[0] = jet.axisMajorQC()-COR_RHO_[ipt][ieta][0]*rho;
    var_[1] = jet.axisMinorQC()-COR_RHO_[ipt][ieta][1]*rho;
    var_[2] = jet.nChargedQC()-COR_RHO_[ipt][ieta][3]*rho;
    var_[3] = jet.fmax()-COR_RHO_[ipt][ieta][4]*rho;
    var_[4] = jet.pullQC()-COR_RHO_[ipt][ieta][2]*rho;
  }
  else {
    var_[0] = jet.axisMajor()-COR_RHO_[ipt][ieta][0]*rho;
    var_[1] = jet.axisMinor()-COR_RHO_[ipt][ieta][1]*rho;
    var_[2] = jet.nChargedPtCut()+jet.nNeutralPtCut()-COR_RHO_[ipt][ieta][3]*rho;
    var_[3] = jet.fmax()-COR_RHO_[ipt][ieta][4]*rho;
    var_[4] = jet.pull()-COR_RHO_[ipt][ieta][2]*rho;
  }
  return reader_[ipt][ieta]->EvaluateMVA("LIK_"+PT_CAT_[ipt]+"_"+ETA_CAT_[ieta]);
}
//-------------------------------------------------------------
QGLCalculator::~QGLCalculator()
{
  delete reader_;
}
