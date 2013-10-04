#define ReadForHisto_cxx
#include "ReadForHisto.h"
#include <iostream>
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include "TMath.h"
using std::cin;
using std::cout;
using std::endl;

void ReadForHisto::Loop(TFile *OUTF, TString SELECTION, bool isMC)
{
  int iINTER(0);
  if (SELECTION == "NOM") {
    iINTER = 1;
  }
  if (SELECTION == "VBF") {
    iINTER = 2;
  }
  //---- define histograms ------------------
  const int NVAR = 17;
  TString var[NVAR] = {"mbb","mbbReg","mqq","dPhibb","dPhiqq","dEtaqq","softHt","nSoftJets","nSoftJets2","nSoftJets5",
                       "cosTheta","BDT","nVtx","nJets","x1_ov_ht","x2_ov_ht","ht"}; 
  double dX[NVAR]   = {2,2,5,0.02,0.02,0.05,1,1,1,1,0.02,0.05,1,1,0.1,0.1,5};
  double XMIN[NVAR] = {0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-15,-15,0};
  double XMAX[NVAR] = {1000,1000,4000,3.14159,3.14159,10,400,30,30,30,1,1,40,40,0,0,4000};
  TH1F *hVar[NVAR],*hVarCut[NVAR];
  cout<<"Booking histograms.........."<<endl;
  for(int ivar=0;ivar<NVAR;ivar++) {
    int NBINS = (XMAX[ivar]-XMIN[ivar])/dX[ivar];
    hVar[ivar]    = new TH1F("h_"+var[ivar],"h_"+var[ivar],NBINS,XMIN[ivar],XMAX[ivar]);
    hVarCut[ivar] = new TH1F("h_"+var[ivar]+"Cut","h_"+var[ivar]+"Cut",NBINS,XMIN[ivar],XMAX[ivar]);
  }
  cout<<"Booking jet histograms.........."<<endl;
  const int NJETVAR = 13;
  TString varJet[NJETVAR] = {"jetPt","jetPtBtag","jetEta","jetEtaBtag","jetPhi","jetBtag","jetQGL","jetChf",
                             "jetNhf","jetPhf","jetMuf","jetElf","jetPuMva"}; 
  int NJETBINS[NJETVAR]   = {200,200,50,50,50,50,50,50,50,50,50,50,100};
  double XMINJET[NJETVAR] = {0,0,-5,-5,-3.14159,0,0,0,0,0,0,0,-1};
  double XMAXJET[NJETVAR] = {1000,1000,5,5,3.14159,1.0001,1,1,1,1,1,1,1};
  TH1F *hJetVar[NJETVAR][4],*hJetVarCut[NJETVAR][4];
  char name[1000];
  for(int ivar=0;ivar<NJETVAR;ivar++) {
    for(int j=0;j<4;j++) {
      sprintf(name,"h_%s%d",varJet[ivar].Data(),j);
      hJetVar[ivar][j]    = new TH1F(name,name,NJETBINS[ivar],XMINJET[ivar],XMAXJET[ivar]);
      sprintf(name,"h_%sCut%d",varJet[ivar].Data(),j);
      hJetVarCut[ivar][j] = new TH1F(name,name,NJETBINS[ivar],XMINJET[ivar],XMAXJET[ivar]);
    }
  }
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"Reading "<<nentries<<" events"<<endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry % (nentries/10) == 0) cout<<100*jentry/nentries<<"%"<<endl;
    fChain->GetEntry(jentry);
    bool cut_trigger(false);
    bool cut_sel(false);
    float mva(-1);
    vector<int> blikIdx;
    if (SELECTION == "NOM") {
      cut_trigger = ((*triggerResult)[0] || (*triggerResult)[1]);  
      cut_sel = selNOM[1];
      mva = mvaNOM;
      blikIdx = *blikNOMIdx;
    }
    if (SELECTION == "VBF") {
      cut_trigger = ((*triggerResult)[9] && !(*triggerResult)[0] && !(*triggerResult)[1]);  
      cut_sel = selVBF[2];
      mva = mvaVBF;
      blikIdx = *blikVBFIdx;
    }
    bool cut_btag = (((*jetBtag)[b1[0]] > 0) && ((*jetBtag)[b2[0]] > 0));
    bool cut_electron = ((*jetElf)[0]<0.7 && (*jetElf)[1]<0.7 && (*jetElf)[2]<0.7 && (*jetElf)[3]<0.7); 
    bool cut_muon     = ((*jetMuf)[0]<0.7 && (*jetMuf)[1]<0.7 && (*jetMuf)[2]<0.7 && (*jetMuf)[3]<0.7);
    if (!cut_trigger)  continue;
    if (!cut_sel)      continue;
    if (!cut_btag)     continue;
    if (dPhibb[iINTER] > 2.0)  continue;
    if (softHt > 200)  continue;
    if (!cut_electron) continue;
    if (!cut_muon)     continue;
    float x[NVAR] = {mbb[iINTER],mbbReg[iINTER],mqq[iINTER],dPhibb[iINTER],dPhiqq[iINTER],dEtaqq[iINTER],softHt,nSoftJets,
                     nSoftJets2,nSoftJets5,fabs(cosTheta[iINTER]),mva,nvtx,nJets,log(x1/ht),log(x2/ht),ht};
    float wt(1.0);
    if (isMC) {
      wt = puWt; 
    } 
    for(int ivar=0;ivar<NVAR;ivar++) {
      hVar[ivar]->Fill(x[ivar],wt);
      if (mva > 0.2) {
        hVarCut[ivar]->Fill(x[ivar],wt);
      }
    }
    for(int j=0;j<4;j++) {
      if ((*jetPt)[j]<0) continue;
      int ib = blikIdx[j]; 
      float xJet[NJETVAR] = {(*jetPt)[j],(*jetPt)[ib],(*jetEta)[j],(*jetEta)[ib],(*jetPhi)[j],(*jetBtag)[ib],(*jetQGL)[ib],(*jetChf)[j],(*jetNhf)[j],(*jetPhf)[j],(*jetMuf)[j],(*jetElf)[j],(*jetPuMva)[ib]}; 
      for(int ivar=0;ivar<NJETVAR;ivar++) {
        hJetVar[ivar][j]->Fill(xJet[ivar],wt);
        if (mva > 0.2) {
          hJetVarCut[ivar][j]->Fill(xJet[ivar],wt);
        }
      }
    }    
  }
  OUTF->Write();
}









