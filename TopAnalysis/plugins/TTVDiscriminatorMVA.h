#ifndef TTVDiscriminatorMVA_h
#define TTVDiscriminatorMVA_h

#include "TMVA/Reader.h"

class TTVDiscriminatorMVA
{
  public:
    TTVDiscriminatorMVA(std::string weights,std::string selection);
    ~TTVDiscriminatorMVA();
    float evalZ(int nJets,int nBJets,float lepPt0,float lepPt1,float lepEta0,float lepEta1,float dRlbMin,float dRlbMax,float dRbbTop,float chi2);
    float evalW(int nJets,int nBJets,float lepPt0,float lepEta0,float met,float dRlbMin,float dRlbMax,float dRbbTop,float mTop,float chi2,float ht);
    
  private:
    std::string weights_;
    std::string selection_;
    TMVA::Reader *reader_;
    float var_[20];
};
#endif
