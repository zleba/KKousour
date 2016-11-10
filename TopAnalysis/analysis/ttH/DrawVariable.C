void DrawVariable(TString VAR,bool SHAPE,int REBIN,float XMIN,float XMAX,TString XTITLE,bool PRINT,bool BLIND=false,float BLIND_MIN=0,float BLIND_MAX=0)
{
  gROOT->ForceStyle();
  const int N = 20;

  TString SAMPLE[N] = {
    "JetHT",
    "ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8",
    "ttHJetToNonbb_M125_13TeV_amcatnloFXFX_madspin_pythia8_mWCutfix",
    "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",
    "WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8",
    "TT_TuneCUETP8M1_13TeV-powheg-pythia8",
    "ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1",
    "ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
    "ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1",
    "TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8",
    "TTZToQQ_TuneCUETP8M1_13TeV-amcatnlo-pythia8",
    "WZ_TuneCUETP8M1_13TeV-pythia8",
    "WWTo4Q_13TeV-powheg",
    "ZZTo4Q_13TeV_amcatnloFXFX_madspin_pythia8"
  };

  float XSEC[N] = {
    1.0,
    0.2934,
    0.2151,
    347700, 
    32100,
    6831,
    1207,
    119.9,
    25.24,
    1460,
    3539,
    832,
    217,
    35.6,
    35.6,
    0.4062,
    0.5297,
    47.13,
    51.723,
    22.29
  }; 
 
  float LUMI(5760);
  TFile *inf[N];
  TH1F  *h[N];

  TCanvas *can = new TCanvas("can_"+VAR,"can_"+VAR,900,600);
  can->SetRightMargin(0.15);

  for(int i=0;i<N;i++) {
    inf[i] = TFile::Open("Histo_"+SAMPLE[i]+".root");
    //cout<<inf[i]->GetName()<<" "<<VAR<<endl;
    h[i]   = (TH1F*)inf[i]->Get("ttH/h_"+VAR);
    h[i]->SetDirectory(0);
    h[i]->Rebin(REBIN);
    h[i]->Sumw2();
    h[i]->SetLineWidth(1);
    h[i]->SetLineColor(kBlack);
    if (i>0) {
      float norm = ((TH1F*)inf[i]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
      //cout<<inf[i]->GetName()<<" "<<norm<<endl; 
      h[i]->Scale(LUMI*XSEC[i]/norm);
    }
    else {
      h[i]->SetLineWidth(2);
    }
    inf[i]->Close();
  }

  TH1F *hQCD = (TH1F*)h[3]->Clone("hQCD");
  hQCD->Add(h[4]);
  hQCD->Add(h[5]);
  hQCD->Add(h[6]); 
  hQCD->Add(h[7]);
  hQCD->Add(h[8]);

  TH1F *hVV = (TH1F*)h[17]->Clone("hVV");
  hVV->Add(h[18]);
  hVV->Add(h[19]);

  TH1F *hST = (TH1F*)h[12]->Clone("hST");
  hST->Add(h[13]);
  hST->Add(h[14]);

  h[1]->SetLineColor(kRed);//ttHbb
  h[1]->SetFillColor(kRed-9);
  h[1]->SetLineWidth(2);
  h[2]->SetLineColor(kRed);//ttHNobb
  h[2]->SetFillColor(kRed-10);
  h[2]->SetLineWidth(2);

  TH1F *hTT     = (TH1F*)h[11]->Clone("TTbar");
  TH1F *hTTZ    = (TH1F*)h[16]->Clone("TTZ");
  TH1F *hTTW    = (TH1F*)h[15]->Clone("TTW");
  TH1F *hWJets  = (TH1F*)h[10]->Clone("WJets");
  TH1F *hDYJets = (TH1F*)h[9]->Clone("DYJets");
 
  hQCD->SetFillColor(kBlue-10);//QCD
  hTT->SetFillColor(kOrange);
  hST->SetFillColor(kOrange-1);//ST
  hVV->SetFillColor(kMagenta-10);//VV
  hTTZ->SetFillColor(kYellow-10);//ttZ
  hTTW->SetFillColor(kYellow-9);//ttW
  hWJets->SetFillColor(kGreen-10);//WJets
  hDYJets->SetFillColor(kGreen-8);//ZJets
 
 

  float kfactor = 1.0;
  if (hQCD->Integral() > 0) { 
    kfactor = (h[0]->Integral()-hTT->Integral()-hWJets->Integral()-hDYJets->Integral()-hST->Integral())/hQCD->Integral();
  }  
  hQCD->Scale(kfactor);

  TH1F *hBkg = (TH1F*)hQCD->Clone("hBkg");
  hBkg->Add(hTT);
  hBkg->Add(hST);
  hBkg->Add(hTTW);
  hBkg->Add(hTTZ);
  hBkg->Add(hWJets);
  hBkg->Add(hDYJets);
  hBkg->Add(hVV);
  hBkg->SetFillColor(kGray);

  cout<<"Data events:  "<<h[0]->Integral()<<endl;
  cout<<"ttHJetTobb:   "<<h[1]->Integral()<<endl;
  cout<<"ttHJetToNonbb:"<<h[2]->Integral()<<endl;
  cout<<"QCD events:   "<<hQCD->Integral()<<endl;
  cout<<"TTbar events: "<<hTT->Integral()<<endl;
  cout<<"ST events:    "<<hST->Integral()<<endl;
  cout<<"TTZ events:   "<<hTTZ->Integral()<<endl;  
  cout<<"TTW events:   "<<hTTW->Integral()<<endl;
  cout<<"WJets events: "<<hWJets->Integral()<<endl;
  cout<<"DYJets events:"<<hDYJets->Integral()<<endl;
  cout<<"VV events:    "<<hVV->Integral()<<endl;
  cout<<"kfactor:      "<<kfactor<<endl;

  THStack *hs = new THStack("hs","hs");
  hs->Add(hVV);
  hs->Add(hTTW);
  hs->Add(hTTZ);
  hs->Add(hWJets);
  hs->Add(hDYJets);
  hs->Add(hST);
  hs->Add(hTT);
  hs->Add(hQCD);

  if (BLIND) {
    for(int i=0;i<h[0]->GetNbinsX();i++) {
      if (h[0]->GetBinCenter(i+1) > BLIND_MIN && h[0]->GetBinCenter(i+1) < BLIND_MAX) {
        h[0]->SetBinContent(i+1,0);
        h[0]->SetBinError(i+1,0);
      } 
    }
  }

  TH1F *hRatio = (TH1F*)h[0]->Clone("Ratio");
  hRatio->SetLineWidth(2);
  hRatio->Divide(hBkg);

  TLegend *leg = new TLegend(0.86,0.55,0.99,0.9);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(hTT,"TTbar","F");
  leg->AddEntry(hST,"ST","F");
  leg->AddEntry(hTTZ,"ttZ","F");
  leg->AddEntry(hTTW,"ttW","F");
  leg->AddEntry(hWJets,"WJets","F");
  leg->AddEntry(hDYJets,"ZJets","F");
  leg->AddEntry(hVV,"Diboson","F");
  leg->AddEntry(h[1],"ttHbb","L");
  leg->AddEntry(h[2],"ttHNonbb","L");
  
  if (SHAPE) {
    h[1]->SetFillStyle(3001);
    h[2]->SetFillStyle(3001);
    h[0]->Scale(1./h[0]->Integral());
    h[1]->Scale(1./h[1]->Integral());
    h[2]->Scale(1./h[2]->Integral());
    hQCD->Scale(1./hQCD->Integral());
    hBkg->Scale(1./hBkg->Integral());
    double max = TMath::Max(h[0]->GetBinContent(h[0]->GetMaximumBin()),hBkg->GetBinContent(hBkg->GetMaximumBin()));
    hBkg->SetMaximum(1.1*max);
    hBkg->GetXaxis()->SetTitle(XTITLE);
    hBkg->GetXaxis()->SetRangeUser(XMIN,XMAX);
    hBkg->Draw("hist");
    h[0]->Draw("same E"); 
    //h[1]->Draw("same hist"); 
    leg->Draw();
    gPad->RedrawAxis();
    can->Print("can_"+VAR+"_norm.pdf"); 
  }
  else {
    can->SetBottomMargin(0.25);
    gPad->SetLogy(); 
    //h[1]->Scale(h[0]->Integral()/h[1]->Integral());
    h[1]->SetFillColor(0);
    TH1F *hAux = (TH1F*)h[0]->Clone("aux");
    hAux->Reset();
    hAux->GetYaxis()->SetRangeUser(0.5,1.1*TMath::Max(h[1]->GetBinContent(h[1]->GetMaximumBin()),h[0]->GetBinContent(h[0]->GetMaximumBin()))); 
    
    hAux->GetXaxis()->SetRangeUser(XMIN,XMAX);
    hAux->GetYaxis()->SetTitle(TString::Format("Number of events / %1.2f fb^{-1}",LUMI/1000));
    hAux->GetXaxis()->SetTitle("");
    hAux->GetXaxis()->SetLabelSize(0.0);
    hAux->Draw();
    hs->Draw("hist same"); 
    h[0]->Draw("same E");
    h[1]->Draw("same hist");
    h[2]->Draw("same hist");
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
    hRatio->Draw();
    if (PRINT) {
      can->Print("plots/can_"+VAR+"_abs.pdf"); 
      can->Print("plots/can_"+VAR+"_abs.png");
    }
  }
}















