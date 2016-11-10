#include <cmath>
#include "KKousour/TopAnalysis/plugins/TTHDiscriminatorMVA.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

TTHDiscriminatorMVA::TTHDiscriminatorMVA(std::string weights)
{
  weights_ = weights;
  reader_  = new TMVA::Reader("!Color:!Silent");

  reader_->AddVariable("nJets"                ,&var_[0]);
  reader_->AddVariable("ht"                   ,&var_[1]);
  reader_->AddVariable("jetPt[5]"             ,&var_[2]);
  reader_->AddVariable("mbbMin"               ,&var_[3]);
  reader_->AddVariable("dRbbMin"              ,&var_[4]);
  reader_->AddVariable("qglMedian"            ,&var_[5]);
  reader_->AddVariable("abs(cosThetaStar1)"   ,&var_[6]); 
  reader_->AddVariable("abs(cosThetaStar2)"   ,&var_[7]);
  reader_->AddVariable("sphericity"           ,&var_[8]);
  reader_->AddVariable("aplanarity"           ,&var_[9]);
  reader_->AddVariable("centrality"           ,&var_[10]); 
  reader_->AddVariable("foxWolfram[0]"        ,&var_[11]);
  reader_->AddVariable("foxWolfram[1]"        ,&var_[12]);
  reader_->AddVariable("foxWolfram[2]"        ,&var_[13]);
  reader_->AddVariable("foxWolfram[3]"        ,&var_[14]);
  reader_->AddVariable("mTop[0]"              ,&var_[15]);
  reader_->AddVariable("ptTTbar"              ,&var_[16]);
  reader_->AddVariable("mTTbar"               ,&var_[17]);
  reader_->AddVariable("dRbbTop"              ,&var_[18]);
  reader_->AddVariable("chi2"                 ,&var_[19]); 
  
  edm::FileInPath f1(weights_);
  reader_->BookMVA("BDT",f1.fullPath());  
}
//-------------------------------------------------------------
float TTHDiscriminatorMVA::eval(int nJets,float ht,float jetPt5,float mbbMin,float dRbbMin,float qglMedian,float cosThetaStar1,float cosThetaStar2,float sphericity,float aplanarity,float centrality,float H0,float H1,float H2,float H3,float mTop,float ptTTbar,float mTTbar,float dRbbTop,float chi2)
{
  var_[0]  = nJets;
  var_[1]  = ht;
  var_[2]  = jetPt5;
  var_[3]  = mbbMin;
  var_[4]  = dRbbMin;
  var_[5]  = qglMedian;
  var_[6]  = fabs(cosThetaStar1);
  var_[7]  = fabs(cosThetaStar2);
  var_[8]  = sphericity;
  var_[9]  = aplanarity;  
  var_[10] = centrality;
  var_[11] = H0;
  var_[12] = H1;
  var_[13] = H2;
  var_[14] = H3;
  var_[15] = mTop;
  var_[16] = ptTTbar;
  var_[17] = mTTbar;
  var_[18] = dRbbTop;
  var_[19] = chi2;
  
  return reader_->EvaluateMVA("BDT");
}
//-------------------------------------------------------------
TTHDiscriminatorMVA::~TTHDiscriminatorMVA()
{
  delete reader_;
}
