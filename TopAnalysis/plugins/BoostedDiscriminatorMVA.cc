#include <cmath>
#include "KKousour/TopAnalysis/plugins/BoostedDiscriminatorMVA.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

BoostedDiscriminatorMVA::BoostedDiscriminatorMVA(std::string weights)
{
  weights_ = weights;
  reader_  = new TMVA::Reader("!Color:!Silent");

  reader_->AddVariable("jetTau3[0]/jetTau2[0]",&var_[0]);
  reader_->AddVariable("jetTau3[0]/jetTau1[0]",&var_[1]);
  reader_->AddVariable("jetTau3[1]/jetTau2[1]",&var_[2]);
  reader_->AddVariable("jetTau3[1]/jetTau1[1]",&var_[3]);

  edm::FileInPath f1(weights_);
  reader_->BookMVA("BDT",f1.fullPath());  
}
//-------------------------------------------------------------
float BoostedDiscriminatorMVA::eval(float tau320,float tau310,float tau321,float tau311)
{
  var_[0] = tau320;
  var_[1] = tau310;
  var_[2] = tau321;
  var_[3] = tau311;

  return reader_->EvaluateMVA("BDT");
}
//-------------------------------------------------------------
BoostedDiscriminatorMVA::~BoostedDiscriminatorMVA()
{
  delete reader_;
}
