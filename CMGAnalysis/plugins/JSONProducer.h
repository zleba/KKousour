#ifndef JSONProducer_h
#define JSONProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <iostream>
#include "TTree.h"

class JSONProducer : public edm::EDAnalyzer 
{
  public:
    explicit JSONProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, const edm::EventSetup& iSetup);
    virtual void endJob();
    virtual ~JSONProducer();

  private:  
    int find(int x, const std::vector<std::pair<int,std::vector<int> > >& v);
    int find(int x, const std::vector<int>& v);
    void write(const std::vector<std::pair<int,std::vector<int> > >& v);
    std::string filename_; 
    ofstream outf_;
    int run_,lumi_;
    std::vector<std::pair<int,std::vector<int> > > vrun_;   
};

#endif
