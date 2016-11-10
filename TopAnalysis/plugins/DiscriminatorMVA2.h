#ifndef DiscriminatorMVA2_h
#define DiscriminatorMVA2_h

#include "TMVA/Reader.h"

class DiscriminatorMVA2
{
  public:
    DiscriminatorMVA2(std::string weights);
    ~DiscriminatorMVA2();
    float eval(int nBJets,float yTop0,float yTop1,float ptTop0,float ptTop1,float yTTbar,float mbbMin,float dRbbAve,float dRbbTop,float EtStar1,float EtStar2,float jetEta0,float jetEta1,float jetEta2,float jetEta3,float jetEta4,float jetEta5);
    
  private:
    std::string weights_;
    TMVA::Reader *reader_;
    float var_[16];
    float spc_[1];
};
#endif
