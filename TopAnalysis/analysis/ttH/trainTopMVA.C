#include "TMVA/Factory.h"
#include "TMVA/MethodCategory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h" 
#include "TH1F.h"
using namespace TMVA;
using namespace TMath;

void trainTopMVA();

void trainTopMVA()
{
  TFile *sigSrc  = TFile::Open("flatTree_ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8_train.root");
  TTree *sigTree = (TTree*)sigSrc->Get("ttH/events"); 
  TFile *bkgSrc  = TFile::Open("flatTree_TT_TuneCUETP8M1_13TeV-powheg-pythia8_train.root");
  TTree *bkgTree = (TTree*)bkgSrc->Get("ttH/events"); 

  char name[1000];
  
  TFile *outf = new TFile("mva_Top_CAT2.root","RECREATE");
  TMVA::Factory* factory = new TMVA::Factory("factory_mva_Top_CAT2_",outf,"!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification");
  
  factory->AddSignalTree(sigTree);
  factory->AddBackgroundTree(bkgTree);
  
  //const int NVAR = 16;
  //TString VAR[NVAR] = {"yTop[0]","yTop[1]","ptTop[0]","ptTop[1]","yTTbar","mbbMin","dRbbAve","dRbbTop","EtStar1","EtStar2","jetEta[0]","jetEta[1]","jetEta[2]","jetEta[3]","jetEta[4]","jetEta[5]"};
  //char TYPE[NVAR] = {'F','F','F','F','F','F','F','F','F','F','F','F','F','F','F','F'};

  const int NVAR = 25;
  TString VAR[NVAR] = {
    "nJets",
    "nBJets",
    "ht",
    "jetPt[2]","jetPt[3]","jetPt[4]","jetPt[5]",
    "mbbMin","dRbbMin",
    "qglMin","qglMedian",
    "cosThetaStar1","cosThetaStar2",
    "sphericity","aplanarity","centrality","foxWolfram[0]","foxWolfram[1]","foxWolfram[2]","foxWolfram[3]",
    "mTop[0]","ptTTbar","mTTbar","dRbbTop","chi2"
  };
  char TYPE[NVAR] = {
    'I',
    'I',
    'F',
    'F','F','F','F', 
    'F','F',
    'F','F','F','F','F','F','F','F','F','F','F',
    'F','F','F','F','F'
  };

  for(int i=0;i<NVAR;i++) {
    factory->AddVariable(VAR[i],TYPE[i]);
  }

  //factory->AddSpectator("nBJets",'I');

  factory->PrepareTrainingAndTestTree("","nTrain_Signal=40000:nTrain_Background=40000:nTest_Signal=40000:nTest_Background=40000");
  
  /*
  TMVA::IMethod* Fisher_Category = factory->BookMethod( TMVA::Types::kCategory,"Fisher_Category");
  TMVA::MethodCategory* mcategory_Fisher = dynamic_cast<TMVA::MethodCategory*>(Fisher_Category);
  
  mcategory_Fisher->AddMethod("nBJets == 2",
                      "mbbMin:dRbbTop:EtStar1:EtStar2:jetEta[0]:jetEta[1]:jetEta[2]:jetEta[3]:jetEta[4]:jetEta[5]",
                      TMVA::Types::kFisher,
                      "Fisher_Cat1","H:!V:Fisher");

  mcategory_Fisher->AddMethod("nBJets > 2",
                      "mbbMin:dRbbTop:EtStar1:EtStar2:jetEta[0]:jetEta[1]:jetEta[2]:jetEta[3]:jetEta[4]:jetEta[5]",
                      TMVA::Types::kFisher,
                      "Fisher_Cat2","H:!V:Fisher");
  
  //----------- BDT ----------------------------------------------------
  TMVA::IMethod* BDT_Category = factory->BookMethod( TMVA::Types::kCategory,"BDT_Category");
  TMVA::MethodCategory* mcategory_BDT = dynamic_cast<TMVA::MethodCategory*>(BDT_Category);
  
  mcategory_BDT->AddMethod("nBJets == 2",
                      "mbbMin:dRbbTop:EtStar1:EtStar2:jetEta[0]:jetEta[1]:jetEta[2]:jetEta[3]:jetEta[4]:jetEta[5]",
                      TMVA::Types::kBDT,
                      "BDT_Cat1","NTrees=1000:BoostType=Grad:Shrinkage=0.1:DoBoostMonitor=True");

  mcategory_BDT->AddMethod("nBJets > 2",
                      "mbbMin:dRbbTop:EtStar1:EtStar2:jetEta[0]:jetEta[1]:jetEta[2]:jetEta[3]:jetEta[4]:jetEta[5]",
                      TMVA::Types::kBDT,
                      "BDT_Cat2","NTrees=1000:BoostType=Grad:Shrinkage=0.1:DoBoostMonitor=True");
  */ 
  factory->BookMethod(TMVA::Types::kFisher,"Fisher");
  //factory->BookMethod(TMVA::Types::kBDT,"BDT","NTrees=500:BoostType=Grad:Shrinkage=0.1:DoBoostMonitor=True");
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods(); 
  outf->Close();
}
