#define TreeClass_cxx
#include "TreeClass.h"
#include <iostream>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDirectory.h>
using std::cin;
using std::cout;
using std::endl;

void TreeClass::Loop(TDirectory *DIR, bool isMC)
{
  //----- get the PU weight --------
  //TFile *puf  = TFile::Open("PileupWeight.root");
  //TH1F *hPuSF = (TH1F*)puf->Get("pileupSF");
  //hPuSF->SetDirectory(0);
  //puf->Close();
  DIR->cd();
  //---- Number of categories ---------------
  const int NCAT = 3;
  //---- define histograms ------------------
  char name[1000];
  const int NVAR = 31;
  TString var[NVAR] = {"mvaQCD","mvaTTbar","ht","htBtag","nvtx","nJets","nBJets","met",
                       "sphericity","aplanarity","FW0","FW1","FW2","FW3",
                       "qglAve","qglMedian","qglMin","dRbbMin","mbbMin","dRbbAve","mbbAve",
                       "mTop[0]","dRbbTop","mTTbar","ptTTbar","chi2",
                       "centrality","cosThetaStar1","cosThetaStar2","EtStar1","EtStar2"}; 
  double dX[NVAR]   = {0.01,0.01,10,10,1,1,1,1,
                       0.01,0.005,0.005,0.005,0.005,0.005,
                       0.01,0.01,0.01,0.1,5,0.1,5,
                       1,0.05,1,1,0.5,
                       0.01,0.01,0.01,1,1};
  double XMIN[NVAR] = {-1,-1,0,0,0,0,0,0,
                       0,0,0.2,-0.3,-0.1,-0.3,
                       0,0,0,0,0,0,0,
                       0,0,0,0,0,
                       0,-1,-1,0,0};
  double XMAX[NVAR] = {1,1,3000,1000,50,15,10,150,
                       1,0.5,0.5,0.2,0.5,0.3,
                       1,1,1,6,800,6,800,
                       3000,7,5000,2000,200,
                       1,1,1,500,500};

  TH2F *hMvaVsMva[NCAT+1]; 
  for(int icat=0;icat<=NCAT;icat++) {
    sprintf(name,"h_MvaVsMva_CAT%d",icat); 
    hMvaVsMva[icat] = new TH2F(name,name,100,-1,1,100,-1,1); 
  }

  TH1F *hCat = new TH1F("hCat","hCat",NCAT+1,0,NCAT+1);
  TH1F *hVar[NVAR][NCAT+1][2][3];
  cout<<"Booking histograms.........."<<endl;
  for(int ivar=0;ivar<NVAR;ivar++) {
    int NBINS = (XMAX[ivar]-XMIN[ivar])/dX[ivar];
    for(int icat=0;icat<=NCAT;icat++) {
      sprintf(name,"h_%s_sigTrg_CAT%d",var[ivar].Data(),icat); 
      hVar[ivar][icat][0][0] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
      sprintf(name,"h_%s_ctlTrg_CAT%d",var[ivar].Data(),icat); 
      hVar[ivar][icat][1][0] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
      sprintf(name,"h_%s_sigTrg_CAT%d_sigReg",var[ivar].Data(),icat); 
      hVar[ivar][icat][0][1] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
      sprintf(name,"h_%s_ctlTrg_CAT%d_sigReg",var[ivar].Data(),icat); 
      hVar[ivar][icat][1][1] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
      sprintf(name,"h_%s_sigTrg_CAT%d_ctlReg",var[ivar].Data(),icat); 
      hVar[ivar][icat][0][2] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
      sprintf(name,"h_%s_ctlTrg_CAT%d_ctlReg",var[ivar].Data(),icat); 
      hVar[ivar][icat][1][2] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
    }
  }
  cout<<"Booking jet histograms.........."<<endl;
  const int NJETVAR = 9;
  TString varJet[NJETVAR] = {"jetPt","jetEta","jetBtag","jetQGL","jetChf","jetNhf","jetPhf","jetMuf","jetElf"}; 
  double dXJET[NJETVAR]   = {2,0.1,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
  double XMINJET[NJETVAR] = {0,-3,-10,0,0,0,0,0,0};
  double XMAXJET[NJETVAR] = {1000,3,1,1,1.0001,1.0001,1.0001,1.0001,1.0001};
  TH1F *hJetVar[NJETVAR][6][NCAT+1][2][3];
  
  for(int ivar=0;ivar<NJETVAR;ivar++) {
    int NBINS = (XMAXJET[ivar]-XMINJET[ivar])/dXJET[ivar];
    for(int j=0;j<6;j++) {
      for(int icat=0;icat<=NCAT;icat++) {
        sprintf(name,"h_%s%d_sigTrg_CAT%d",varJet[ivar].Data(),j,icat);
        hJetVar[ivar][j][icat][0][0] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"h_%s%d_ctlTrg_CAT%d",varJet[ivar].Data(),j,icat);
        hJetVar[ivar][j][icat][1][0] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"h_%s%d_sigTrg_CAT%d_sigReg",varJet[ivar].Data(),j,icat);
        hJetVar[ivar][j][icat][0][1] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"h_%s%d_ctlTrg_CAT%d_sigReg",varJet[ivar].Data(),j,icat);
        hJetVar[ivar][j][icat][1][1] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"h_%s%d_sigTrg_CAT%d_ctlReg",varJet[ivar].Data(),j,icat);
        hJetVar[ivar][j][icat][0][2] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"h_%s%d_ctlTrg_CAT%d_ctlReg",varJet[ivar].Data(),j,icat);
        hJetVar[ivar][j][icat][1][2] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
      }
    }
  }
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"Reading "<<nentries<<" events"<<endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry % (nentries/10) == 0) cout<<100*jentry/nentries<<"%"<<endl;
    fChain->GetEntry(jentry);
    bool cut_ctl_trigger(true);
    bool cut_sig_trigger(true);
    if (!isMC) {
      cut_ctl_trigger = ((*triggerBit)[1] || (*triggerBit)[3]);
      cut_sig_trigger = ((*triggerBit)[0] || (*triggerBit)[2]);
    }
    bool cut_leptons = (nLeptons == 0); 
    bool cut_ht      = (ht > 450);
    bool cut_jetPt   = ((*jetPt)[5] > 40);
    bool cut_nJets   = (nJets > 5);
    bool cut_nBJets  = (nBJets > 1);
    bool cut_status  = (status == 0);
    bool cut_sigReg(false);
    bool cut_ctlReg(false);

    if (!cut_ctl_trigger && !cut_sig_trigger)  continue;
    if (!cut_leptons)  continue;
    if (!cut_ht)       continue;
    if (!cut_jetPt)    continue;
    if (!cut_nJets)    continue;
    if (!cut_nBJets)   continue;
    if (!cut_status)   continue;
    
    //---- gen weight -----------------
    float wtGen = 1.0;

    //---- pu weight ----------------
    float wtPU = 1.0; 
    if (isMC) {
      wtGen = genEvtWeight;
      /*
      int pubin = hPuSF->FindBin(npu);
      if (pubin < 40) {
        wtPU = hPuSF->GetBinContent(pubin);
      }
      else {
        wtPU = hPuSF->GetBinContent(40);
      }
      */
    }

    int category(0);
    
    if (nBJets == 2) {
      category = 1;
      cut_sigReg = (mvaQCD > -0.5 && mvaTTbar > 0.0);
      cut_ctlReg = (mvaTTbar < 0.0); 
    } 
    if (nBJets == 3) {
      category = 2;
      cut_sigReg = (mvaQCD > -0.5 && mvaTTbar > 0.0);
      cut_ctlReg = (mvaTTbar < 0.0); 
    }
    if (nBJets > 3) {
      category = 3;
      cut_sigReg = (mvaQCD > -0.5 && mvaTTbar > 0.0);
      cut_ctlReg = (mvaTTbar < 0.0);
    }
    
    hCat->Fill(category);
    float x[NVAR] = {mvaQCD,mvaTTbar,ht,htBtag,float(nvtx),float(nJets),float(nBJets),met,
                     sphericity,aplanarity,foxWolfram[0],foxWolfram[1],foxWolfram[2],foxWolfram[3],
                     qglAve,qglMedian,qglMin,dRbbMin,mbbMin,dRbbAve,mbbAve,mTop[0],dRbbTop,mTTbar,ptTTbar,chi2,
                     centrality,cosThetaStar1,cosThetaStar2,EtStar1,EtStar2};
     
    hMvaVsMva[0]->Fill(mvaQCD,mvaTTbar,wtGen*wtPU); 
    hMvaVsMva[category]->Fill(mvaQCD,mvaTTbar,wtGen*wtPU);

    for(int ivar=0;ivar<NVAR;ivar++) {
      if (cut_sig_trigger) {
        hVar[ivar][0][0][0]->Fill(x[ivar],wtGen*wtPU);
        hVar[ivar][category][0][0]->Fill(x[ivar],wtGen*wtPU); 
        if (cut_sigReg) {
          hVar[ivar][0][0][1]->Fill(x[ivar],wtGen*wtPU);
          hVar[ivar][category][0][1]->Fill(x[ivar],wtGen*wtPU);
        }
        if (cut_ctlReg) {
          hVar[ivar][0][0][2]->Fill(x[ivar],wtGen*wtPU);
          hVar[ivar][category][0][2]->Fill(x[ivar],wtGen*wtPU);
        } 
      }
      if (cut_ctl_trigger) {
        hVar[ivar][0][1][0]->Fill(x[ivar],wtGen*wtPU);
        hVar[ivar][category][1][0]->Fill(x[ivar],wtGen*wtPU); 
        if (cut_sigReg) {
          hVar[ivar][0][1][1]->Fill(x[ivar],wtGen*wtPU);
          hVar[ivar][category][1][1]->Fill(x[ivar],wtGen*wtPU);
        }
        if (cut_ctlReg) {
          hVar[ivar][0][1][2]->Fill(x[ivar],wtGen*wtPU);
          hVar[ivar][category][1][2]->Fill(x[ivar],wtGen*wtPU);
        }
      }
    }
    for(int j=0;j<6;j++) {
      float xJet[NJETVAR] = {(*jetPt)[j],(*jetEta)[j],(*jetBtag)[j],(*jetQGL)[j],(*jetChf)[j],(*jetNhf)[j],(*jetPhf)[j],(*jetMuf)[j],(*jetElf)[j],}; 
      for(int ivar=0;ivar<NJETVAR;ivar++) {
        if (cut_sig_trigger) {
          hJetVar[ivar][j][0][0][0]->Fill(xJet[ivar],wtGen*wtPU);
          hJetVar[ivar][j][category][0][0]->Fill(xJet[ivar],wtGen*wtPU);
          if (cut_sigReg) {
            hJetVar[ivar][j][0][0][1]->Fill(xJet[ivar],wtGen*wtPU);
            hJetVar[ivar][j][category][0][1]->Fill(xJet[ivar],wtGen*wtPU);
          }
          if (cut_ctlReg) {
            hJetVar[ivar][j][0][0][2]->Fill(xJet[ivar],wtGen*wtPU);
            hJetVar[ivar][j][category][0][2]->Fill(xJet[ivar],wtGen*wtPU);
          }
        }
        if (cut_ctl_trigger) {
          hJetVar[ivar][j][0][1][0]->Fill(xJet[ivar],wtGen*wtPU);
          hJetVar[ivar][j][category][1][0]->Fill(xJet[ivar],wtGen*wtPU);
          if (cut_sigReg) {
            hJetVar[ivar][j][0][1][1]->Fill(xJet[ivar],wtGen*wtPU);
            hJetVar[ivar][j][category][1][1]->Fill(xJet[ivar],wtGen*wtPU);
          }
          if (cut_ctlReg) {
            hJetVar[ivar][j][0][1][2]->Fill(xJet[ivar],wtGen*wtPU);
            hJetVar[ivar][j][category][1][2]->Fill(xJet[ivar],wtGen*wtPU);
          }
        }
      }
    }    
  }
  DIR->Write();
}
