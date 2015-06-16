#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TMath.h"
#include "TLorentzVector.h"

#include "KKousour/TopAnalysis/plugins/AllHadronicPartonFilter.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace std;
using namespace reco;

AllHadronicPartonFilter::AllHadronicPartonFilter(edm::ParameterSet const& cfg) 
{
  srcGenParticles_ = cfg.getParameter<edm::InputTag>("genparticles");
  forceTopDecay_   = cfg.getParameter<bool>("forceTopDecay");
  forceHiggsDecay_ = cfg.getParameter<bool>("forceHiggsDecay");
}
//////////////////////////////////////////////////////////////////////////////////////////
bool AllHadronicPartonFilter::filter(edm::Event& iEvent, edm::EventSetup const& iSetup) 
{
  edm::Handle<GenParticleCollection> genParticles;
  iEvent.getByLabel(edm::InputTag(srcGenParticles_),genParticles);

  bool WPlusLep(false),WMinusLep(false),HToBB(false),result(true);
  for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
    const GenParticle &p = (*genParticles)[ip];
    
    if (p.pdgId() == 24) {
      for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
        int daughterID = p.daughter(k)->pdgId();
        if (daughterID == -11 || daughterID == -13 || daughterID == -15) {
          WPlusLep = true;
        }
      }
    }
    if (p.pdgId() == -24) {
      for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
        int daughterID = p.daughter(k)->pdgId();
        if (daughterID == 11 || daughterID == 13 || daughterID == 15) {
          WMinusLep = true;
        }
      }
    }
    if (p.pdgId() == 25) {
      for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
        int daughterID = p.daughter(k)->pdgId();
        if (fabs(daughterID) == 5) {
          HToBB = true;
        }
      }
    }
  }// end of particle loop
  if (forceTopDecay_) {
    result *= (!WPlusLep && !WMinusLep);
  }
  if (forceHiggsDecay_) {
    result *= HToBB;
  }
  return result;
}  
//////////////////////////////////////////////////////////////////////////////////////////
AllHadronicPartonFilter::~AllHadronicPartonFilter() 
{
}

DEFINE_FWK_MODULE(AllHadronicPartonFilter);
