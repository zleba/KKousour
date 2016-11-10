#ifndef DiscriminatorMVA_h
#define DiscriminatorMVA_h

#include "TMVA/Reader.h"

class DiscriminatorMVA
{
  public:
    DiscriminatorMVA(std::string weights);
    ~DiscriminatorMVA();
    float eval(int nBJets,float ht,float jetPt2,float jetPt3,float jetPt4,float jetPt5,float qglMin,float qglMedian,float sphericity,float aplanarity,float H0,float H1,float H2,float H3,float mTop,float ptTTbar,float mTTbar,float chi2,float centrality,float cosThetaStar1,float cosThetaStar2);
    
  private:
    std::string weights_;
    TMVA::Reader *reader_;
    float var_[21];
    float spc_[1];
};
#endif
