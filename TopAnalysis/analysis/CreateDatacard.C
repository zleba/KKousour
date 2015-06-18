void CreateDatacard()
{
  gROOT->ForceStyle();
  TString FileName[5] = {"TTH","TTJets","QCD250","QCD500","QCD1000"};
  float XSEC[5]       = {0.46*0.577*0.5085,0.46*424.5,670500,26740,769.7};
  float LUMI(10000);
  const int NCAT = 2;
  TFile *inf[5];
  TTree *tr[5];
  TH1F  *h[5][NCAT];
  int COLOR[5] = {kRed,kBlack,kBlue,kBlue-5,kBlue-8};
  char name[1000];
  TCut CUT[NCAT];
  CUT[0] = "ht>500 && jetPt[5]>40 && nBJets==2";
  CUT[1] = "ht>500 && jetPt[5]>40 && nBJets>2";

  TCanvas *can = new TCanvas("can","can",900,600);
  can->Divide(2,1);
  for(int i=0;i<5;i++) {
    inf[i] = TFile::Open("flatTree_"+FileName[i]+".root");
    tr[i]  = (TTree*)inf[i]->Get("hadtop/events");
    TH1F *hpu = (TH1F*)inf[i]->Get("hadtop/pileup");
    for(int k=0;k<NCAT;k++) {
      can->cd(k+1);
      cout<<i<<" "<<k<<endl;
      h[i][k] = new TH1F(TString::Format("h%d%d",i,k),TString::Format("h%d%d",i,k),25,-1,1);
      h[i][k]->Sumw2();
      tr[i]->Draw("mva>>"+TString::Format("h%d%d",i,k),CUT[k]);
      h[i][k]->SetLineWidth(2);
      h[i][k]->SetLineColor(COLOR[i]);
      h[i][k]->Scale(LUMI*XSEC[i]/hpu->GetEntries());
    } 
  }
  TFile *outf = new TFile("ttH-shapes.root","RECREATE");
  TH1F *hQCD[NCAT];
  for(int k=0;k<NCAT;k++) {
    can->cd(k+1); 
    hQCD[k] = (TH1F*)h[2][k]->Clone(TString::Format("hQCD%d",k));
    hQCD[k]->Add(h[3][k]);
    hQCD[k]->Add(h[4][k]); 
    h[0][k]->SetFillColor(kRed-10);
    h[1][k]->SetFillColor(kGreen-10);
    hQCD[k]->SetFillColor(kBlue-10);
    h[0][k]->SetFillStyle(3001);
    h[1][k]->SetFillStyle(3001);
    cout<<"QCD events:   "<<hQCD[k]->Integral()<<endl;
    cout<<"TTbar events: "<<h[1][k]->Integral()<<endl;
    cout<<"TTH events:   "<<h[0][k]->Integral()<<endl;
    double max1 = TMath::Max(h[1][k]->GetBinContent(h[1][k]->GetMaximumBin()),hQCD[k]->GetBinContent(hQCD[k]->GetMaximumBin()));
    float max = TMath::Max(h[0][k]->GetBinContent(h[0][k]->GetMaximumBin()),max1);
    hQCD[k]->SetMaximum(1.1*max);
    hQCD[k]->Draw("hist");
    h[0][k]->Draw("same hist"); 
    h[1][k]->Draw("same hist"); 
  
    outf->cd();
    hQCD[k]->Write(TString::Format("qcdCAT%d",k));
    h[0][k]->Write(TString::Format("signalCAT%d",k));
    h[1][k]->Write(TString::Format("ttjetsCAT%d",k));
    hQCD[k]->Write(TString::Format("data_obsCAT%d",k));
  }

  ofstream datacard;
  datacard.open("datacard_ttH.txt");
  datacard<<"imax "<<NCAT<<"\n";
  datacard<<"jmax *"<<"\n";
  datacard<<"kmax *"<<"\n";
  datacard<<"----------------"<<"\n";
  datacard<<"shapes *  * "<<outf->GetName()<<" $PROCESS$CHANNEL"<<"\n";
  datacard<<"----------------"<<"\n";
  datacard<<"bin         ";
  for(int k=0;k<NCAT;k++) { 
    sprintf(name,"CAT%d ",k);
    datacard<<name;
  }
  datacard<<"\n";
  datacard<<"observation ";
  for(int k=0;k<NCAT;k++) {
    datacard<<"-1 ";
  }  
  datacard<<"\n";
  datacard<<"----------------"<<"\n";
  datacard<<"bin  ";
  for(int k=0;k<NCAT;k++) {
    sprintf(name,"CAT%d CAT%d CAT%d ",k,k,k);
    datacard<<name;
  }  
  datacard<<"\n";
  datacard<<"process ";
  for(int k=0;k<NCAT;k++) {
    datacard<<"signal ttjets qcd ";
  }  
  datacard<<"\n";
  datacard<<"process ";
  for(int k=0;k<NCAT;k++) {
    datacard<<"0 1 2 ";
  }  
  datacard<<"\n";
  datacard<<"rate       ";
  for(int k=0;k<NCAT;k++) {
    datacard<<h[0][k]->Integral()<<" "<<h[1][k]->Integral()<<" "<<hQCD[k]->Integral()<<" ";
  }
  datacard<<"\n";
  datacard.close();
  //outf->Close();
}















