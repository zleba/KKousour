#ifndef DiscriminatorMVA_h
#define DiscriminatorMVA_h

#include "TMVA/Reader.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

class DiscriminatorMVA
{
  public:
    DiscriminatorMVA(std::string weights);
    ~DiscriminatorMVA();
    float eval(float mqq,float dEtaqq,float dPhiqq,float btag0,float btag1,float qgl0,float qgl1,float qgl2,float qgl3,float softHt,int nSoftJets,float cosTheta,float x1,float x2,float sphericity,float aplanarity,float etaRatio,float ptRatio);
    
  private:
    std::string weights_;
    TMVA::Reader *reader_;
    float var_[18];
};
#endif
