#include "KKousour/CMGAnalysis/plugins/JetRegressor.h"

JetRegressor::JetRegressor(std::string weights)
{
  weights_ = weights;
  reader_ = new TMVA::Reader("!Color:!Silent");
  /*
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
  */
  reader_->AddVariable("jetPt"          ,&var_[0]);
  reader_->AddVariable("jetEta"         ,&var_[1]);
  reader_->AddVariable("jetM"           ,&var_[2]);
  reader_->AddVariable("met"            ,&var_[3]);
  reader_->AddVariable("jetMetPhi"      ,&var_[4]);
  reader_->AddVariable("jetLeadTrkPt"   ,&var_[5]);
  reader_->AddVariable("rho"            ,&var_[6]);
  reader_->AddVariable("jetPart"        ,&var_[7]);
  reader_->AddVariable("jetPhf"         ,&var_[8]);
  reader_->AddVariable("jetNhf"         ,&var_[9]);
  reader_->AddVariable("jetVtx3deL"     ,&var_[10]);
  reader_->AddVariable("jetVtxMass"     ,&var_[11]);
  reader_->AddVariable("jetSoftLepPt"   ,&var_[12]);
  reader_->AddVariable("jetSoftLepPtRel",&var_[13]);
  
  edm::FileInPath f1(weights_);
  reader_->BookMVA("BDT",f1.fullPath());  
  
}
//-------------------------------------------------------------
std::vector<float> JetRegressor::getTarget(cmg::PFJet const& jet,float jetMetPhi,float rho,float met,float softLepPt,float softLepPtRel)
{
  var_[0]  = jet.pt(); 
  var_[1]  = jet.eta();
  var_[2]  = jet.mass();
  var_[3]  = met;
  var_[4]  = jetMetPhi;
  var_[5]  = jet.fmaxCharged()*jet.pt()*jet.rawFactor();
  var_[6]  = rho;
  var_[7]  = jet.nConstituents();
  var_[8]  = jet.component(4).fraction();
  var_[9]  = jet.component(5).fraction();
  var_[10] = jet.vtx3deL();
  var_[11] = jet.secvtxMass();
  var_[12] = softLepPt;
  var_[13] = softLepPtRel;

  std::vector<float> result(reader_->EvaluateRegression("BDT"));
  return result;
}
//-------------------------------------------------------------
JetRegressor::~JetRegressor()
{
  delete reader_;
}
