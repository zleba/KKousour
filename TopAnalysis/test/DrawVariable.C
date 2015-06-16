void DrawVariable(TString SEL,TString VAR,bool SHAPE,int NBINS,float XMIN,float XMAX,TString XTITLE)
{
  gROOT->ForceStyle();
  TString FileName[5] = {"TTH","TTJets","QCD250","QCD500","QCD1000"};
  float XSEC[5]       = {0.46*0.577*0.5085,0.46*424.5,670500,26740,769.7};
  float LUMI(10000);
  TFile *inf[5];
  TTree *tr[5];
  TH1F  *h[5];
  int COLOR[5] = {kRed,kBlack,kBlue,kBlue-5,kBlue-8};
  TCut CUT(SEL);

  TCanvas *can = new TCanvas("can_"+VAR,"can_"+VAR,900,600);
  can->SetRightMargin(0.15);

  for(int i=0;i<5;i++) {
    inf[i] = TFile::Open("flatTree_"+FileName[i]+".root");
    tr[i]  = (TTree*)inf[i]->Get("hadtop/events");
    TH1F *hpu = (TH1F*)inf[i]->Get("hadtop/pileup");
    h[i]   = new TH1F("h_"+VAR+"_"+FileName[i],"h_"+VAR+"_"+FileName[i],NBINS,XMIN,XMAX);
    h[i]->Sumw2();
    tr[i]->Draw(VAR+">>"+"h_"+VAR+"_"+FileName[i],CUT);
    h[i]->SetLineWidth(2);
    h[i]->SetLineColor(COLOR[i]);
    h[i]->Scale(LUMI*XSEC[i]/hpu->GetEntries());
  }

  TH1F *hQCD = (TH1F*)h[2]->Clone("hQCD");
  hQCD->Add(h[3]);
  hQCD->Add(h[4]); 

  h[0]->SetFillColor(kRed-10);
  h[1]->SetFillColor(kGreen-10);
  hQCD->SetFillColor(kBlue-10);

  THStack *hs = new THStack("hs","hs");
  hs->Add(h[1]);
  hs->Add(hQCD);

  cout<<"QCD events:   "<<hQCD->Integral()<<endl;
  cout<<"TTbar events: "<<h[1]->Integral()<<endl;
  cout<<"TTH events:   "<<h[0]->Integral()<<endl;

  TLegend *leg = new TLegend(0.86,0.65,0.99,0.9);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(hQCD,"QCD","F");
  leg->AddEntry(h[1],"TTbar","F");
  leg->AddEntry(h[0],"ttH","F");
  
  if (SHAPE) {
    h[0]->SetFillStyle(3001);
    h[1]->SetFillStyle(3001);
    h[0]->Scale(1./h[0]->Integral());
    h[1]->Scale(1./h[1]->Integral());
    hQCD->Scale(1./hQCD->Integral());
    double max1 = TMath::Max(h[1]->GetBinContent(h[1]->GetMaximumBin()),hQCD->GetBinContent(hQCD->GetMaximumBin()));
    float max = TMath::Max(h[0]->GetBinContent(h[0]->GetMaximumBin()),max1);
    hQCD->SetMaximum(1.1*max);
    hQCD->GetXaxis()->SetTitle(XTITLE);
    hQCD->Draw("hist");
    h[0]->Draw("same hist"); 
    h[1]->Draw("same hist"); 
    leg->Draw();
    gPad->RedrawAxis();
    can->Print("can_"+VAR+"_norm.pdf"); 
  }
  else {
    gPad->SetLogy(); 
    TH1F *hAux = (TH1F*)h[0]->Clone("aux");
    hAux->Reset();
    hAux->GetYaxis()->SetRangeUser(0.5,1.1*(hQCD->GetBinContent(hQCD->GetMaximumBin())+h[1]->GetBinContent(h[1]->GetMaximumBin()))); 
    hAux->GetXaxis()->SetTitle(XTITLE);
    hAux->GetYaxis()->SetTitle(TString::Format("Number of events / %d fb^{-1}",int(LUMI)/1000));
    hAux->Draw();
    hs->Draw("hist same"); 
    h[0]->Draw("same hist");
    leg->Draw();
    gPad->RedrawAxis();
    can->Print("can_"+VAR+"_abs.pdf"); 
  }
}















