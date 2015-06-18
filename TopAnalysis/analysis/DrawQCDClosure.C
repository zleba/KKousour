void DrawQCDClosure(TString SEL,TString VAR,int NBINS,float XMIN,float XMAX,TString XTITLE)
{
  gROOT->ForceStyle();
  TString FileName[3] = {"QCD250","QCD500","QCD1000"};
  float XSEC[3]       = {670500,26740,769.7};
  TFile *inf[3];
  TTree *tr[3],*tr1[3];
  TH1F  *h[3],*h1[3];
  TCut CUT(SEL);

  TCanvas *can = new TCanvas("can_QCDClosure_"+VAR,"can_QCDClosure_"+VAR,900,600);
  can->cd(1);
  can->SetBottomMargin(0.3);
  can->SetRightMargin(0.15);

  for(int i=0;i<2;i++) {
    inf[i] = TFile::Open("flatTree_"+FileName[i]+".root");
    tr[i]  = (TTree*)inf[i]->Get("hadtopNoBtag/events");
    tr1[i] = (TTree*)inf[i]->Get("hadtop/events");
    TH1F *hpu = (TH1F*)inf[i]->Get("hadtop/pileup");
    h[i]   = new TH1F("h_"+FileName[i],"h_"+FileName[i],NBINS,XMIN,XMAX);
    h1[i]  = new TH1F("h1_"+FileName[i],"h1_"+FileName[i],NBINS,XMIN,XMAX);
    h[i]->Sumw2();
    h1[i]->Sumw2();
    //tr[i]->Draw(VAR+">>"+"h_"+FileName[i],"(1.11-0.000162845*ht)*(nBJets>1 && mva>-0.5)");
    tr[i]->Draw(VAR+">>"+"h_"+FileName[i],"nBJets>1 && mva>-1 && ht>400 && jetPt[5]>20 && chi2<10");
    tr1[i]->Draw(VAR+">>"+"h1_"+FileName[i],"nBJets>1 && mva>-1 && ht>400 && jetPt[5]>20 && chi2<10");
    h[i]->Scale(XSEC[i]/hpu->GetEntries());
    h1[i]->Scale(XSEC[i]/hpu->GetEntries());
    cout<<hpu->GetEntries()<<endl;
  }

  TH1F *hQCD = (TH1F*)h[0]->Clone("hQCD");
  hQCD->Add(h[1]);
  //hQCD->Add(h[2]); 
  hQCD->SetFillColor(kGray);
  TH1F *hQCD1 = (TH1F*)h1[0]->Clone("hQCD1");
  hQCD1->Add(h1[1]);
  //hQCD1->Add(h1[2]);

  hQCD->Scale(1./hQCD->Integral());
  hQCD1->Scale(1./hQCD1->Integral());
  hQCD->GetXaxis()->SetLabelSize(0.0);
  double max = 1.1*TMath::Max(hQCD->GetBinContent(hQCD->GetMaximumBin()),hQCD1->GetBinContent(hQCD1->GetMaximumBin()));
  hQCD->SetMinimum(1e-5);
  hQCD->SetMaximum(max);
  hQCD->Draw("hist");
  hQCD1->Draw("sameE");
  gPad->RedrawAxis();

  TLegend *leg = new TLegend(0.86,0.65,0.99,0.9);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(hQCD,"Control","F");
  leg->AddEntry(hQCD1,"Signal","LP");
  leg->Draw();

  TH1F *hRatio = (TH1F*)hQCD1->Clone("Ratio");
  hRatio->Divide(hQCD);
  TF1 *fit = new TF1("fit","[0]+[1]*TMath::Erf([2]*x+[3])",200,2000);
  fit->SetParameters(0.5,0.5,0.01,1);

  TPad* pad = new TPad("pad", "pad", 0., 0., 1., 1.);
  pad->SetTopMargin(0.7);
  pad->SetRightMargin(0.15);
  pad->SetFillColor(0);
  pad->SetFillStyle(0);
  pad->Draw();
  pad->cd(0);
  gPad->SetGridy();
  hRatio->GetXaxis()->SetTitle(XTITLE);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetRangeUser(0,2);
  hRatio->GetYaxis()->SetLabelSize(0.04);
 // hRatio->Fit(fit,"R");
  hRatio->Draw();
}















