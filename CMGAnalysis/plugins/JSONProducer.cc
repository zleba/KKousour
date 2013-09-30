#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>

#include "KKousour/CMGAnalysis/plugins/JSONProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"


static bool pair_rule(const std::pair<int,std::vector<int> >& i, const std::pair<int,std::vector<int> >& j) 
{
  return (i.first < j.first);
}
//////////////////////////////////////////////////////////////////////////////////////////
JSONProducer::JSONProducer(edm::ParameterSet const& cfg) 
{
  filename_ = cfg.getParameter<std::string> ("filename");
}
//////////////////////////////////////////////////////////////////////////////////////////
void JSONProducer::beginJob() 
{
  //--- create the text file ----------------
  outf_.open(filename_);
  //--- clean the vector --------------------
  vrun_.clear();
}
//////////////////////////////////////////////////////////////////////////////////////////
void JSONProducer::write(const std::vector<std::pair<int,std::vector<int> > >& v)
{
  outf_<<"{";
  for(unsigned iRun=0;iRun<vrun_.size();iRun++) {
    outf_<<"\""<<vrun_[iRun].first<<"\": [";
    int l1 = vrun_[iRun].second[0];
    int prev = l1;
    for(unsigned iLumi=0;iLumi<vrun_[iRun].second.size();iLumi++) {
      int l2 = vrun_[iRun].second[iLumi]; 
      if ((iLumi == vrun_[iRun].second.size()-1) && (l2 == prev)) {
        outf_<<"["<<l1<<", "<<l1<<"]";
      }
      if ((iLumi == vrun_[iRun].second.size()-1) && (l2-prev == 1)) {
        outf_<<"["<<l1<<", "<<l2<<"]";
      }   
      if (l2-prev > 1) {
        if (iLumi == vrun_[iRun].second.size()-1) {
          outf_<<"["<<l1<<", "<<prev<<"],";
          outf_<<"["<<l2<<", "<<l2<<"]"; 
        }
        else {
          outf_<<"["<<l1<<", "<<prev<<"],"; 
        }  
        l1 = l2;
      }
      prev = l2;
    }
    outf_<<"]";
    if (iRun < vrun_.size()-1) {
      outf_<<","; 
    }
  }
  outf_<<"}";
}
//////////////////////////////////////////////////////////////////////////////////////////
void JSONProducer::endJob() 
{
  std::sort(vrun_.begin(),vrun_.end(),pair_rule);
  for(unsigned iRun=0;iRun<vrun_.size();iRun++) {
    std::cout<<"Run: "<<vrun_[iRun].first<<std::endl;
    std::sort(vrun_[iRun].second.begin(),vrun_[iRun].second.end());
    for(unsigned iLumi=0;iLumi<vrun_[iRun].second.size();iLumi++) {
      std::cout<<vrun_[iRun].second[iLumi]<<" ";
    }
    std::cout<<std::endl;
  }
  write(vrun_);
}
//////////////////////////////////////////////////////////////////////////////////////////
void JSONProducer::analyze(edm::Event const& iEvent, const edm::EventSetup& iSetup) 
{
  run_    = iEvent.id().run();
  lumi_   = iEvent.id().luminosityBlock();

  int iRun = find(run_,vrun_);
  if (iRun == -1) {
    std::vector<int> vlumi(0);
    vlumi.push_back(lumi_);
    std::pair<int,std::vector<int> > p;
    p.first = run_;
    p.second = vlumi;
    vrun_.push_back(p);
    iRun = vrun_.size()-1;
  }
  int iLumi = find(lumi_,vrun_[iRun].second);
  if (iLumi == -1) {
    vrun_[iRun].second.push_back(lumi_);
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
int JSONProducer::find(int x, const std::vector<std::pair<int,std::vector<int> > >& v)
{
  int result(-1);
  for(unsigned i=0;i<v.size();i++)
    if (x == v[i].first)
      return i;
  return result;
}
//////////////////////////////////////////////////////////////////////////////////////////
int JSONProducer::find(int x, const std::vector<int>& v)
{
  int result(-1);
  for(unsigned i=0;i<v.size();i++)
    if (x == v[i])
      return i;
  return result;
}
//////////////////////////////////////////////////////////////////////////////////////////
JSONProducer::~JSONProducer() 
{
}

DEFINE_FWK_MODULE(JSONProducer);
