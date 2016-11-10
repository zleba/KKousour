//#include "TreeClassResolved.C"
#include "TreeClassBoosted.C"
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TString.h> 
#include <TCollection.h>
#include <TKey.h>
void FillHistograms(TString SAMPLE,bool isMC,bool doBtagSF,bool applyNLOWeight)
{
  cout<<"Processing sample "<<SAMPLE<<endl;
  TString PATH("../../prod/ttbar/");
  TFile *inf  = TFile::Open(PATH+"flatTree_"+SAMPLE+".root");
  TFile *outf = TFile::Open(TString::Format("Histo_%s.root",SAMPLE.Data()),"RECREATE");

  TIter nextKey(inf->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)nextKey())) {
    TString dirName(key->GetName());
    cout<<"Found directory "<<dirName<<endl;
    
    outf->mkdir(dirName);  
    TDirectory *dir = (TDirectory*)outf->Get(dirName); 
    TTree *tr   = (TTree*)inf->Get(dirName+"/events");
    if (dirName.Contains("eventCounter")) {
      dir->cd();
      TH1F *hGenEventWeight = (TH1F*)inf->Get(dirName+"/GenEventWeight");
      hGenEventWeight->Write("GenEventWeight");  
    } 
    if (dirName.Contains("boosted")) {
      TreeClassBoosted myTree(tr);
      dir->cd();
      myTree.Loop(dir,isMC,doBtagSF,applyNLOWeight);
      cout<<"Loop finished"<<endl;
      dir->Close();
      cout<<"directory closed"<<endl;
      delete tr;
      cout<<"Tree deleted"<<endl;
    }
  }
  outf->Close();
  inf->Close();
}
