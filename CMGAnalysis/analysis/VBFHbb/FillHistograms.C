#include "ReadForHisto.C"
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TString.h> 
void FillHistograms(TString FILENAME,TString SELECTION,bool isMC)
{
  TString PATH("root://eoscms//eos/cms/store/cmst3/group/vbfhbb/flat/flatTree_");
  TFile *inf  = TFile::Open(PATH+FILENAME+".root");
  TTree *tr   = (TTree*)inf->Get("Hbb/events");
  TH1F *hPass = (TH1F*)inf->Get("Hbb/TriggerPass");
  TFile *outf = TFile::Open("Histo_sel"+SELECTION+"_"+FILENAME+".root","RECREATE");
  outf->cd();
  hPass->Write("TriggerPass");
  ReadForHisto myTree(tr);
  myTree.Loop(outf,SELECTION,isMC);
  
  outf->Close();
  inf->Close();
}
