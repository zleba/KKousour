#include <cmath>
#include "KKousour/CMGAnalysis/plugins/DiscriminatorMVA.h"

DiscriminatorMVA::DiscriminatorMVA(std::string weights)
{
  weights_ = weights;
  reader_ = new TMVA::Reader("!Color:!Silent");
  
  reader_->AddVariable("mqq"            ,&var_[0]);
  reader_->AddVariable("dEtaqq"         ,&var_[1]);
  reader_->AddVariable("dPhiqq"         ,&var_[2]);
  reader_->AddVariable("btag_[0]"       ,&var_[3]);
  reader_->AddVariable("btag_[1]"       ,&var_[4]);
  reader_->AddVariable("qgl_[0]"        ,&var_[5]);
  reader_->AddVariable("qgl_[1]"        ,&var_[6]);
  reader_->AddVariable("qgl_[2]"        ,&var_[7]);
  reader_->AddVariable("qgl_[3]"        ,&var_[8]);
  reader_->AddVariable("softHt"         ,&var_[9]);
  reader_->AddVariable("nSoftJets2"     ,&var_[10]);
  reader_->AddVariable("cosTheta"       ,&var_[11]);
  reader_->AddVariable("log(x1)"        ,&var_[12]);
  reader_->AddVariable("log(x2)"        ,&var_[13]);
  reader_->AddVariable("log(sphericity)",&var_[14]);
  reader_->AddVariable("log(aplanarity)",&var_[15]);
  reader_->AddVariable("etaRatio"       ,&var_[16]);
  reader_->AddVariable("ptRatio"        ,&var_[17]);
  
  edm::FileInPath f1(weights_);
  reader_->BookMVA("BDT",f1.fullPath());  
  
}
//-------------------------------------------------------------
float DiscriminatorMVA::eval(float mqq,float dEtaqq,float dPhiqq,float btag0,float btag1,float qgl0,float qgl1,float qgl2,float qgl3,float softHt,int nSoftJets,float cosTheta,float x1,float x2,float sphericity,float aplanarity,float etaRatio,float ptRatio)
{
  var_[0]  = mqq;
  var_[1]  = dEtaqq;
  var_[2]  = dPhiqq;
  var_[3]  = btag0;
  var_[4]  = btag1;
  var_[5]  = qgl0;
  var_[6]  = qgl1;
  var_[7]  = qgl2;
  var_[8]  = qgl3;
  var_[9]  = softHt;
  var_[10] = nSoftJets;
  var_[11] = cosTheta;
  var_[12] = log(x1);
  var_[13] = log(x2);  
  var_[14] = log(sphericity);
  var_[15] = log(aplanarity);
  var_[16] = etaRatio;
  var_[17] = ptRatio;

  return reader_->EvaluateMVA("BDT");
}
//-------------------------------------------------------------
DiscriminatorMVA::~DiscriminatorMVA()
{
  delete reader_;
}
