#include "TMVA/Factory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h" 
using namespace TMVA;
using namespace TMath;

void trainMVA();

void trainMVA()
{
  char name[1000];
  float XSEC[4] = {670500,26740,769.7,0.46*424.5};
  float NORM[4];
  TCut preselectionCut = "nBJets>1 && qglAve>0 && ht>500 && jetPt[5]>40";
  TFile *bkgSrc[4];
  bkgSrc[0] = TFile::Open("flatTree_QCD250.root");
  bkgSrc[1] = TFile::Open("flatTree_QCD500.root");
  bkgSrc[2] = TFile::Open("flatTree_QCD1000.root");
  //bkgSrc[3] = TFile::Open("flatTree_TTJets.root");

  TFile *sigSrc = TFile::Open("flatTree_TTH.root");
  TTree *sigTree = (TTree*)sigSrc->Get("hadtop/events"); 
  TTree *bkgTree[4];
  
  
  TFile *outf = new TFile("mva_QCD.root","RECREATE");
  TMVA::Factory* factory = new TMVA::Factory("factory_mva_QCD_",outf,"!V:!Silent:Color:DrawProgressBar:Transformations=I;G:AnalysisType=Classification");
  factory->AddSignalTree(sigTree);

  for(int k=0;k<3;k++) {
    NORM[k] = ((TH1F*)bkgSrc[k]->Get("hadtop/pileup"))->GetEntries();
    bkgTree[k] = (TTree*)bkgSrc[k]->Get("hadtop/events");
    factory->AddBackgroundTree(bkgTree[k],XSEC[k]/NORM[k]);
  }
  
  int N_SIG(sigTree->GetEntries(preselectionCut));
  
  int N_BKG0(bkgTree[0]->GetEntries(preselectionCut));
  int N_BKG1(bkgTree[1]->GetEntries(preselectionCut));
  int N_BKG2(bkgTree[2]->GetEntries(preselectionCut));
  //int N_BKG3(bkgTree[3]->GetEntries(preselectionCut));

  //float N_BKG_EFF = N_BKG0*XSEC[0]/NORM[0]+N_BKG1*XSEC[1]/NORM[1]+N_BKG2*XSEC[2]/NORM[2]+N_BKG3*XSEC[3]/NORM[3];
  
  //int N = TMath::Min((float)N_SIG,N_BKG_EFF);

  //cout<<N_SIG<<" "<<N_BKG_EFF<<endl;
  
  const int NVAR = 20;
  TString VAR[NVAR] = {
    "nJets",
    //"nBJets",
    "ht","htBtag",
    "jetPt[0]","jetPt[1]","jetPt[2]","jetPt[3]","jetPt[4]","jetPt[5]",
    "mbbMin","dRbbMin",
    //"dRbbAve","mbbAve",
    //"btagAve","btagMax","btagMin",
    "qglAve","qglMin","qglMedian",
    "sphericity","aplanarity","foxWolfram[0]","foxWolfram[1]","foxWolfram[2]","foxWolfram[3]"
    //"mTop[0]","(ptTop[0]+ptTop[1])/2","ptTTbar","mTTbar","dRbbTop","chi2"
  };
  char TYPE[NVAR] = {
    'I',
    //'I',
    'F','F',
    'F','F','F','F','F','F', 
    'F','F',
    //'F','F',
    //'F','F','F',
    'F','F','F',
    'F','F','F','F','F','F' 
    //'F','F','F','F','F','F'
  };
   
  for(int i=0;i<NVAR;i++) {
    factory->AddVariable(VAR[i],TYPE[i]);
  }
  
  // spectator variables: not used for the training but recorded
  //factory->AddSpectator("nJets",'I');
  //factory->AddSpectator("nBJets",'I');

  sprintf(name,"nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d",-1,-1,-1,-1);
  factory->PrepareTrainingAndTestTree(preselectionCut,name);
  // specify the training methods
  factory->BookMethod(TMVA::Types::kFisher,"Fisher");
  //factory->BookMethod(TMVA::Types::kLikelihood,"Likelihood");
  factory->BookMethod(TMVA::Types::kBDT,"BDT_GRAD_500","NTrees=500:BoostType=Grad:Shrinkage=0.1");
  factory->BookMethod(TMVA::Types::kBDT,"BDT_GRAD_1000","NTrees=1000:BoostType=Grad:Shrinkage=0.1");
  //factory->BookMethod(TMVA::Types::kBDT,"BDT_GRAD_2000","NTrees=2000:BoostType=Grad:Shrinkage=0.1");
  //factory->BookMethod(TMVA::Types::kMLP,"MLP_ANN","NCycles=500:HiddenLayers=N,N-1:TrainingMethod=BFGS");
  
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods(); 
  outf->Close();
}
