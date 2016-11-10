void DrawVariable(TString DIR,TString VAR,float LUMI,bool LOG,int REBIN,float XMIN,float XMAX,TString XTITLE,bool isINT,int XNDIV,bool PRINT)
{
  gROOT->ForceStyle();
  const int N = 14;

  TString SAMPLE[N] = {
    "JetHT",
    "TT_TuneCUETP8M1_13TeV-powheg-pythia8",
    "WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",
    "DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",
    "ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
    "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
    "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
  };
  float XSEC[N] = {1.0,0.5*832,3539,1460.,136.02,80.95,35.6,35.6,3.67e+5,2.94e+4,6.524e+03,1.064e+03,121.5,2.542e+01};

  TFile *inf[N];
  TH1F  *h[N];

  TCanvas *can = new TCanvas("DataVsMC_"+DIR+"_"+VAR,"DataVsMC_"+DIR+"_"+VAR,900,600);
  can->SetRightMargin(0.15);

  for(int i=0;i<N;i++) {
    inf[i] = TFile::Open("Histo_"+SAMPLE[i]+".root");
    h[i] = (TH1F*)inf[i]->Get(DIR+"/hWt_"+VAR);
    if (!h[i]) {
      cout<<"Histogram "<<"hWt_"+VAR<<" does not exist !!!"<<endl;
      break;
    } 
    h[i]->SetDirectory(0);
    h[i]->Sumw2();
    h[i]->Rebin(REBIN);
    h[i]->SetLineWidth(1);
    h[i]->SetLineColor(kBlack);
    if (i>0) {
      float norm = ((TH1F*)inf[i]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();   
      //cout<<SAMPLE[i]<<" "<<norm<<endl;
      h[i]->Scale(LUMI*XSEC[i]/norm);
    }
    inf[i]->Close();
  }

  TH1F *hQCD = (TH1F*)h[8]->Clone("hQCD");
  for(int i=9;i<N;i++) {
    hQCD->Add(h[i]);
  }
 
  TH1F *hST = (TH1F*)h[4]->Clone("hST");
  hST->Add(h[5]);
  hST->Add(h[6]);
  hST->Add(h[7]);
  h[0]->SetLineWidth(2);//data
  hQCD->SetFillColor(kBlue-10);//QCD
  h[1]->SetFillColor(kOrange);//ttbar
  h[2]->SetFillColor(kGreen-10);//WJets
  h[3]->SetFillColor(kGreen-8);//ZJets
  hST->SetFillColor(kOrange-1);//ST

  float kfactor = 1.0;
  if (hQCD->Integral() > 0) { 
    kfactor = (h[0]->Integral()-h[1]->Integral()-h[2]->Integral()-h[3]->Integral()-hST->Integral())/hQCD->Integral();
  }  
  hQCD->Scale(kfactor);

  TH1F *hBkg = (TH1F*)hQCD->Clone("hBkg");
  hBkg->Add(h[1]);
  hBkg->Add(h[2]);
  hBkg->Add(h[3]);
  hBkg->Add(h[4]);
  hBkg->Add(hST);
  //hBkg->SetFillColor(kGray);

  cout<<"======== "<<VAR<<"====================="<<endl;
  cout<<"Data events:  "<<h[0]->Integral()<<endl;
  cout<<"QCD events:   "<<hQCD->Integral()<<endl;
  cout<<"WJets events: "<<h[2]->Integral()<<endl;
  cout<<"ZJets events: "<<h[3]->Integral()<<endl;
  cout<<"ST events:    "<<hST->Integral()<<endl;
  cout<<"TTbar events: "<<h[1]->Integral()<<endl;
  cout<<"kfactor:      "<<kfactor<<endl;

  THStack *hs = new THStack("hs","hs");
  if (LOG) {   
    hs->Add(h[2]);
    hs->Add(h[3]);
    hs->Add(hST); 
    hs->Add(hQCD);
    hs->Add(h[1]);
  }
  else {
    hs->Add(h[3]);
    hs->Add(hST);
    hs->Add(h[2]);
    hs->Add(hQCD);
    hs->Add(h[1]);
  }
  TH1F *hRatio = (TH1F*)h[0]->Clone("Ratio");
  hRatio->SetLineWidth(2);
  hRatio->Divide(hBkg);

  TLegend *leg = new TLegend(0.86,0.7,0.99,0.9);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(h[1],"TTbar","F"); 
  leg->AddEntry(hST,"ST","F"); 
  leg->AddEntry(h[2],"WJets","F");
  leg->AddEntry(h[3],"ZJets","F"); 
  
  can->SetBottomMargin(0.25);
  TH1F *hAux = (TH1F*)h[0]->Clone("aux");
  hAux->Reset();
  hAux->GetXaxis()->SetNdivisions(XNDIV);
  if (isINT) {
    hAux->GetXaxis()->CenterLabels();
  } 
  hAux->GetYaxis()->SetRangeUser(0.5,1.1*TMath::Max(hBkg->GetBinContent(hBkg->GetMaximumBin()),h[0]->GetBinContent(h[0]->GetMaximumBin())));  
  if (LOG) {
    gPad->SetLogy(); 
    hAux->GetYaxis()->SetRangeUser(0.5,2*TMath::Max(hBkg->GetBinContent(hBkg->GetMaximumBin()),h[0]->GetBinContent(h[0]->GetMaximumBin())));
  }
  hAux->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hAux->GetYaxis()->SetTitle(TString::Format("Number of events / %1.2f fb^{-1}",LUMI/1000));
  hAux->GetXaxis()->SetTitle("");
  hAux->GetXaxis()->SetLabelSize(0.0);
  hAux->Draw();
  hs->Draw("hist same");
  //hBkg->Draw("sames hist"); 
  h[0]->Draw("sames E");
  
  leg->Draw();
  gPad->RedrawAxis();

  TPad *pad = new TPad("pad","pad",0.,0.,1.,1.);
  pad->SetTopMargin(0.77);
  pad->SetRightMargin(0.15);
  pad->SetFillColor(0);
  pad->SetFillStyle(0);
  pad->Draw();
  pad->cd(0);
  pad->SetGridy();
  hRatio->GetXaxis()->SetTitleOffset(0.95);
  hRatio->GetYaxis()->SetTitleOffset(1.5);
  hRatio->GetYaxis()->SetTickLength(0.06);
  hRatio->GetYaxis()->SetTitleSize(0.03);
  hRatio->GetYaxis()->SetLabelSize(0.03);
  hRatio->GetYaxis()->SetTitle("Data/MC");
  hRatio->GetXaxis()->SetTitle(XTITLE);
  hRatio->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hRatio->GetYaxis()->SetRangeUser(0.5,1.5);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetXaxis()->SetNdivisions(XNDIV);
  if (isINT) {
    hRatio->GetXaxis()->CenterLabels();
  }
  hRatio->Draw();
  if (PRINT) {
    can->Print("plots/"+TString(can->GetName())+".pdf"); 
    can->Print("plots/"+TString(can->GetName())+".png");
  }
}















