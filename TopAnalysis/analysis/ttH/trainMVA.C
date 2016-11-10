#include "TMVA/Factory.h"
#include "TMVA/MethodCategory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h" 
#include "TH1F.h"
using namespace TMVA;
using namespace TMath;

void trainMVA(TString TYPE_BKG, TString CAT);

void trainMVA(TString TYPE_BKG, TString CAT)
{
  TFile *sigSrc    = TFile::Open("flatTree_ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8_train.root");
  TTree *sigTree   = (TTree*)sigSrc->Get("ttH/events"); 
   
  char name[1000];
  
  TFile *outf = new TFile("mva_"+TYPE_BKG+"_"+CAT+".root","RECREATE");
  TMVA::Factory* factory = new TMVA::Factory("TTH_mva_"+TYPE_BKG+"_"+CAT,outf,"!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification");  

  factory->AddSignalTree(sigTree);
  //factory->SetWeightExpression("genEvtWeight");

  if (TYPE_BKG == "TTbar") {
    TFile *bkgSrc  = TFile::Open("flatTree_TT_TuneCUETP8M1_13TeV-powheg-pythia8_train.root");
    TTree *bkgTree = (TTree*)bkgSrc->Get("ttH/events");
    factory->AddBackgroundTree(bkgTree);
  }
  else {
    float XSEC[6] = {347700,32100,6831,1207,119.9,25.24}; 
    TFile *bkgSrc[6];
    TTree *bkgTree[6];
    bkgSrc[0] = TFile::Open("flatTree_QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_train.root");
    bkgSrc[1] = TFile::Open("flatTree_QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_train.root");
    bkgSrc[2] = TFile::Open("flatTree_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_train.root");
    bkgSrc[3] = TFile::Open("flatTree_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_train.root");
    bkgSrc[4] = TFile::Open("flatTree_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_train.root");
    bkgSrc[5] = TFile::Open("flatTree_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_train.root");
    for(int k=0;k<6;k++) {
      float norm = ((TH1F*)bkgSrc[k]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
      bkgTree[k] = (TTree*)bkgSrc[k]->Get("ttH/events");
      factory->AddBackgroundTree(bkgTree[k],XSEC[k]/norm);
    }
  }  
  
  const int NVAR = 20;
  TString VAR[NVAR] = {
    "nJets",
    "ht",
    "jetPt[5]",
    "mbbMin","dRbbMin",
    "qglMedian",
    "abs(cosThetaStar1)","abs(cosThetaStar2)",
    "sphericity","aplanarity","centrality","foxWolfram[0]","foxWolfram[1]","foxWolfram[2]","foxWolfram[3]",
    "mTop[0]","ptTTbar","mTTbar","dRbbTop","chi2"
  };
  char TYPE[NVAR] = {
    'I',
    'F', 
    'F','F',
    'F','F','F','F','F','F','F','F','F','F','F',
    'F','F','F','F','F'
  };

  for(int i=0;i<NVAR;i++) {
    factory->AddVariable(VAR[i],TYPE[i]);
  }
  if (CAT == "CAT0") {
    factory->PrepareTrainingAndTestTree("nBJets==2","nTrain_Signal=200000:nTrain_Background=200000:nTest_Signal=200000:nTest_Background=200000");
  } 
  else {
    if (TYPE_BKG == "QCD") {
      factory->PrepareTrainingAndTestTree("nBJets>2","nTrain_Signal=30000:nTrain_Background=30000:nTest_Signal=27000:nTest_Background=27000"); 
    }
    else {
      factory->PrepareTrainingAndTestTree("nBJets>2","nTrain_Signal=200000:nTrain_Background=200000:nTest_Signal=200000:nTest_Background=200000");
    } 
  }
  factory->BookMethod(TMVA::Types::kFisher,"Fisher");
  factory->BookMethod(TMVA::Types::kBDT,"BDT","NTrees=200:BoostType=Grad:Shrinkage=0.1:DoBoostMonitor=True");
  //factory->BookMethod(TMVA::Types::kBDT,"BDT");
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods(); 
  outf->Close();
}
