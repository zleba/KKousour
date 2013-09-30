#ifndef JetRegressor_h
#define JetRegressor_h

#include "TMVA/Reader.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"

class JetRegressor
{
  public:
    JetRegressor(std::string weights);
    ~JetRegressor();
    std::vector<float> getTarget(cmg::PFJet const& jet,float jetMetPhi,float rho,float met);
    
  private:
    std::string weights_;
    TMVA::Reader *reader_;
    float var_[15];
};
#endif
