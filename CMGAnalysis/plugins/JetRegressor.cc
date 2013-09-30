#include "KKousour/CMGAnalysis/plugins/JetRegressor.h"

JetRegressor::JetRegressor(std::string weights)
{
  weights_ = weights;
  reader_ = new TMVA::Reader("!Color:!Silent");
  
  reader_->AddVariable("jetBtag"   ,&var_[0]);
  reader_->AddVariable("jetPt"     ,&var_[1]);
  reader_->AddVariable("jetEta"    ,&var_[2]);
  reader_->AddVariable("jetMetPhi" ,&var_[3]);
  reader_->AddVariable("jetChf"    ,&var_[4]);
  reader_->AddVariable("jetPhf"    ,&var_[5]);
  reader_->AddVariable("jetNhf"    ,&var_[6]);
  reader_->AddVariable("jetElf"    ,&var_[7]);
  reader_->AddVariable("jetMuf"    ,&var_[8]);
  reader_->AddVariable("jetPtD"    ,&var_[9]);
  reader_->AddVariable("jetVtxPt"  ,&var_[10]);
  reader_->AddVariable("jetVtx3dL" ,&var_[11]);
  reader_->AddVariable("jetVtx3deL",&var_[12]);
  reader_->AddVariable("met"       ,&var_[13]);
  reader_->AddVariable("rho"       ,&var_[14]);
  
  edm::FileInPath f1(weights_);
  reader_->BookMVA("MLP",f1.fullPath());  
  
}
//-------------------------------------------------------------
std::vector<float> JetRegressor::getTarget(cmg::PFJet const& jet,float jetMetPhi,float rho,float met)
{
  var_[0]  = jet.bDiscriminator("combinedSecondaryVertexBJetTags");
  var_[1]  = jet.pt(); 
  var_[2]  = jet.eta();
  var_[3]  = jetMetPhi;
  var_[4]  = jet.component(1).fraction();
  var_[5]  = jet.component(4).fraction();
  var_[6]  = jet.component(5).fraction();
  var_[7]  = jet.component(2).fraction();
  var_[8]  = jet.component(3).fraction();
  var_[9]  = jet.ptd();
  var_[10] = jet.vtxPt();
  var_[11] = jet.vtx3dL();
  var_[12] = jet.vtx3deL();
  var_[13] = met;
  var_[14] = rho;
  std::vector<float> result(reader_->EvaluateRegression("MLP"));
  return result;
}
//-------------------------------------------------------------
JetRegressor::~JetRegressor()
{
  delete reader_;
}
