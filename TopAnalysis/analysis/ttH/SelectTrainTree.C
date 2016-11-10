#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TTree.h"
#include "TCut.h"
#include "TDirectoryFile.h"
#include <iostream>
void SelectTrainFile(TString FILENAME)
{
  TString PATH("root://eoscms//eos/cms/store/cmst3/user/kkousour/ttH/flat/");
  //TString PATH("/afs/cern.ch/work/k/kkousour/private/CMSSW_7_4_12/src/KKousour/TopAnalysis/prod/ttH/");
  TFile *inf = TFile::Open(PATH+"flatTree_"+FILENAME+".root");
  TH1F *hGenEventWeight = (TH1F*)inf->Get("eventCounter/GenEventWeight");
  TTree *tIN = (TTree*)inf->Get("ttH/events");

  tIN->SetBranchStatus("*",0);
  tIN->SetBranchStatus("nJets",1);
  tIN->SetBranchStatus("nBJets",1);
  tIN->SetBranchStatus("nLeptons",1);
  tIN->SetBranchStatus("status",1);
  tIN->SetBranchStatus("ht",1);
  tIN->SetBranchStatus("mbbMin",1);
  tIN->SetBranchStatus("mbbAve",1);
  tIN->SetBranchStatus("dRbbMin",1);
  tIN->SetBranchStatus("dRbbAve",1);
  tIN->SetBranchStatus("qglAve",1);
  tIN->SetBranchStatus("qglMin",1);
  tIN->SetBranchStatus("qglMedian",1); 
  tIN->SetBranchStatus("btagMax",1);
  tIN->SetBranchStatus("btagAve",1);
  tIN->SetBranchStatus("sphericity",1);
  tIN->SetBranchStatus("aplanarity",1);
  tIN->SetBranchStatus("centrality",1);
  tIN->SetBranchStatus("cosThetaStar1",1);
  tIN->SetBranchStatus("cosThetaStar2",1);
  tIN->SetBranchStatus("EtStar1",1);
  tIN->SetBranchStatus("EtStar2",1);
  tIN->SetBranchStatus("foxWolfram",1);
  tIN->SetBranchStatus("hcMoments",1);
  tIN->SetBranchStatus("yTop",1);
  tIN->SetBranchStatus("ptTop",1);
  tIN->SetBranchStatus("mTop",1);
  tIN->SetBranchStatus("ptTTbar",1);
  tIN->SetBranchStatus("mTTbar",1);
  tIN->SetBranchStatus("yTTbar",1);
  tIN->SetBranchStatus("dRbbTop",1);
  tIN->SetBranchStatus("chi2",1);
  tIN->SetBranchStatus("prob",1);
  tIN->SetBranchStatus("genEvtWeight",1);
  tIN->SetBranchStatus("jetPt",1); 
  tIN->SetBranchStatus("jetEta",1);
  tIN->SetBranchStatus("triggerBit",1);

  TFile *outf = TFile::Open("flatTree_"+FILENAME+"_train.root","RECREATE");
  TCut SEL = "ht>450 && jetPt[5]>40 && nBJets>1 && nLeptons==0 && status==0"; 
  
  TDirectoryFile *dir = (TDirectoryFile*)outf->mkdir("ttH");
  dir->cd();
  TTree *tOUT = (TTree*)tIN->CopyTree(SEL);
  cout<<"Events: "<<tOUT->GetEntries()<<endl;
  tOUT->Write("events");
  dir->Close();

  dir = (TDirectoryFile*)outf->mkdir("eventCounter");
  dir->cd();
  hGenEventWeight->Write("GenEventWeight");
  dir->Close(); 
  
  outf->Close();
  inf->Close();
}
