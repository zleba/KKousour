#include "../plugins/JEC.h"

#include <iostream>

using namespace edm;

using namespace std;

void GetJEC()
{
  string jecTag = "Summer16_07Aug2017";
  int version = 4;
  char period = 'B';
  string jetType = "AK4PFchs";
  vector<string> dumy;
  bool isMC = false;

  JECs jetEcorrs;
  jetEcorrs.Init(isMC, jecTag, period, version, jetType, "", dumy);

  //Jet.eta(0.5);

  double rho = 0.2;
  double CorFactor, Unc;
  double pt = 300;
  double eta = 3;
  double area = 0.1;
  //double jet = jetEcorrs.JEC_CHScorrections(pt, eta, area, rho, dumy, CorFactor, Unc);

  double corr = jetEcorrs.GetJECL2L3Residual(pt, eta);


  for(eta = -6; eta <= 6; eta += 0.1) {
  for(pt = 100; pt < 2600; pt += 400) {
      //cout << pt <<" "<<eta << " "<< jetEcorrs.GetJECL2L3Residual(pt, eta) << endl;
      //cout << pt <<" "<<eta << " "<< jetEcorrs.GetJECL2Relative(pt, eta) << endl;
      cout << pt <<" "<<eta << " "<< jetEcorrs.GetJECL3Absolute(pt, eta) << endl;
  }
  }

}

int main()
{

    GetJEC();

    return 0;
}
