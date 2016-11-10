#ifndef BoostedDiscriminatorMVA_h
#define BoostedDiscriminatorMVA_h

#include "TMVA/Reader.h"

class BoostedDiscriminatorMVA
{
  public:
    BoostedDiscriminatorMVA(std::string weights);
    ~BoostedDiscriminatorMVA();
    float eval(float tau320,float tau310,float tau321,float tau311);
    
  private:
    std::string weights_;
    TMVA::Reader *reader_;
    float var_[4];
};
#endif
