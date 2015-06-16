#ifndef AllHadronicPartonFilter_h
#define AllHadronicPartonFilter_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Particle.h"

class AllHadronicPartonFilter : public edm::EDFilter 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    AllHadronicPartonFilter(edm::ParameterSet const& cfg);
    virtual bool filter(edm::Event& iEvent, edm::EventSetup const& iSetup);
    ~AllHadronicPartonFilter();

  private:  
    //---- configurable parameters --------  
    edm::InputTag srcGenParticles_;
    bool forceTopDecay_,forceHiggsDecay_;
};

#endif
