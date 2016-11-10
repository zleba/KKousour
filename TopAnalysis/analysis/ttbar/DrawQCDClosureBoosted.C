TGraphErrors *getUnc(TMatrixDSym COV, TF1 *fit, TH1F *h);
void FilterHisto(TH1F *h);
void DrawQCDClosureBoosted(TString VAR,TString CUT,TString XTITLE,int REBIN,float XMIN,float XMAX,TString FUNC,bool LOG)
{
  gROOT->ForceStyle();

  const int N = 6;
  float XSEC[N] = {3.513e+5,3.1630e+4,6.802e+03,1.206e+03,120.4,2.524e+01};

  TString SAMPLE[N] = {
    "QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8", 
    "QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8" 
  }; 

  TFile *inf[N];
  TH1F  *h0[N],*h1[N],*h2[N];
  TH1F  *hQCD0,*hQCD1,*hQCD2,*hRatio1,*hRatio2; 
  
  for(int i=0;i<N;i++) {
    inf[i] = TFile::Open("Histo_"+SAMPLE[i]+".root");
    h0[i] = (TH1F*)inf[i]->Get("boosted/hWt_"+VAR+"_"+CUT+"_0btag");
    h1[i] = (TH1F*)inf[i]->Get("boosted/hWt_"+VAR+"_"+CUT+"_1btag");
    h2[i] = (TH1F*)inf[i]->Get("boosted/hWt_"+VAR+"_"+CUT+"_2btag");
    h0[i]->Rebin(REBIN);
    h1[i]->Rebin(REBIN);
    h2[i]->Rebin(REBIN);
    h0[i]->Sumw2();
    h1[i]->Sumw2();
    h2[i]->Sumw2();
    float norm = ((TH1F*)inf[i]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights(); 
    h0[i]->Scale(XSEC[i]/norm);
    h1[i]->Scale(XSEC[i]/norm);
    h2[i]->Scale(XSEC[i]/norm); 
    //cout<<i<<" "<<h2[i]->GetEntries()<<" "<<h1[i]->Integral()<<endl;
  }
  hQCD0 = (TH1F*)h0[0]->Clone("hQCD0");
  hQCD1 = (TH1F*)h1[0]->Clone("hQCD1");
  hQCD2 = (TH1F*)h2[0]->Clone("hQCD2");
  for(int i=1;i<N;i++) {
    hQCD0->Add(h0[i]);
    hQCD1->Add(h1[i]);
    hQCD2->Add(h2[i]);
  }  
  hQCD0->SetLineColor(kBlack);
  hQCD0->SetMarkerColor(kBlack);
  hQCD0->SetFillColor(kGray);
  hQCD1->SetLineColor(kRed);
  hQCD1->SetMarkerColor(kRed);
  hQCD1->SetLineWidth(2);
  hQCD1->SetMarkerStyle(20);
  hQCD2->SetLineColor(kBlue);
  hQCD2->SetMarkerColor(kBlue);
  hQCD2->SetLineWidth(2);
  hQCD2->SetMarkerStyle(25);
  hQCD0->Scale(1./hQCD0->Integral());
  hQCD1->Scale(1./hQCD1->Integral());
  hQCD2->Scale(1./hQCD2->Integral());

  TString ss1 = "QCDClosure_Boosted_"+CUT+"_1btag_"+VAR;
  TCanvas *can1 = new TCanvas(ss1,ss1,900,600); 
  can1->cd(1);
  can1->SetBottomMargin(0.3);
  can1->SetRightMargin(0.15);
  if (LOG) {
    gPad->SetLogy();
  }
  hQCD0->GetXaxis()->SetLabelSize(0.0);
  double max1 = 1.3*TMath::Max(hQCD0->GetBinContent(hQCD0->GetMaximumBin()),hQCD1->GetBinContent(hQCD1->GetMaximumBin()));
  hQCD0->SetMinimum(1e-5);
  hQCD0->SetMaximum(max1);
  hQCD0->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hQCD0->Draw("hist");
  hQCD1->Draw("sameE");
  gPad->RedrawAxis();

  TLegend *leg1 = new TLegend(0.86,0.75,0.99,0.9);
  leg1->AddEntry(hQCD0,"0 btag","F");
  leg1->AddEntry(hQCD1,"1 btag","LP");
  leg1->SetFillColor(0);
  leg1->SetTextFont(42);
  leg1->SetTextSize(0.03);
  leg1->Draw();

  hRatio1 = (TH1F*)hQCD1->Clone("hRatio1");
  hRatio1->Divide(hQCD0);
  FilterHisto(hRatio1);

  TPad *pad1 = new TPad("pad1","pad1", 0., 0., 1., 1.);
  pad1->SetTopMargin(0.7);
  pad1->SetRightMargin(0.15);
  pad1->SetFillColor(0);
  pad1->SetFillStyle(0);
  pad1->Draw();
  pad1->cd(0);
  gPad->SetGridy();
  hRatio1->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hRatio1->GetXaxis()->SetTitle(XTITLE);
  hRatio1->GetXaxis()->SetTitleOffset(1.0);
  hRatio1->GetYaxis()->SetNdivisions(505);
  hRatio1->GetYaxis()->SetRangeUser(0,2);
  hRatio1->GetYaxis()->SetLabelSize(0.04);
  hRatio1->Draw();

  TVirtualFitter *fitter1;
  TMatrixDSym COV1;
  
  TF1 *fit1 = new TF1("fit_"+CUT+"_1btag",FUNC,XMIN,XMAX);
  fit1->SetLineColor(kRed); 
  hRatio1->Fit(fit1,"RQ+"); 
  
  fitter1 = TVirtualFitter::GetFitter();
  COV1.Use(fitter1->GetNumberTotalParameters(),fitter1->GetCovarianceMatrix());
  TGraphErrors *gUnc1 = getUnc(COV1,fit1,hRatio1); 
  gUnc1->SetFillColor(kRed);
  gUnc1->SetFillStyle(3001);
  gUnc1->Draw("sameE3");

  can1->Print("plots/"+TString(can1->GetName())+".pdf");
  can1->Print("plots/"+TString(can1->GetName())+".png");

  TString ss2 = "QCDClosure_Boosted_"+CUT+"_2btag_"+VAR;
  TCanvas *can2 = new TCanvas(ss2,ss2,900,600); 
  can2->cd(1);
  can2->SetBottomMargin(0.3);
  can2->SetRightMargin(0.15);
  if (LOG) {
    gPad->SetLogy();
  }
  hQCD0->GetXaxis()->SetLabelSize(0.0);
  double max2 = 1.3*TMath::Max(hQCD0->GetBinContent(hQCD0->GetMaximumBin()),hQCD2->GetBinContent(hQCD2->GetMaximumBin()));
  hQCD0->SetMinimum(1e-5);
  hQCD0->SetMaximum(max2);
  hQCD0->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hQCD0->Draw("hist");
  hQCD2->Draw("sameE");
  gPad->RedrawAxis();

  TLegend *leg2 = new TLegend(0.86,0.75,0.99,0.9);
  leg2->AddEntry(hQCD0,"0 btag","F");
  leg2->AddEntry(hQCD2,"2 btag","LP");
  leg2->SetFillColor(0);
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.03);
  leg2->Draw();

  hRatio2 = (TH1F*)hQCD2->Clone("hRatio2");
  hRatio2->Divide(hQCD0);
  FilterHisto(hRatio2);

  TPad *pad2 = new TPad("pad2","pad2", 0., 0., 1., 1.);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.15);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);
  pad2->Draw();
  pad2->cd(0);
  gPad->SetGridy();
  hRatio2->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hRatio2->GetXaxis()->SetTitle(XTITLE);
  hRatio2->GetXaxis()->SetTitleOffset(1.0);
  hRatio2->GetYaxis()->SetNdivisions(505);
  hRatio2->GetYaxis()->SetRangeUser(0,2);
  hRatio2->GetYaxis()->SetLabelSize(0.04);
  hRatio2->Draw();

  TVirtualFitter *fitter2;
  TMatrixDSym COV2;
  
  TF1 *fit2 = new TF1("fit_"+CUT+"_2btag",FUNC,XMIN,XMAX);
  fit2->SetLineColor(kBlue); 
  hRatio2->Fit(fit2,"R+"); 
  
  fitter2 = TVirtualFitter::GetFitter();
  COV2.Use(fitter2->GetNumberTotalParameters(),fitter2->GetCovarianceMatrix());
  TGraphErrors *gUnc2 = getUnc(COV2,fit2,hRatio2); 
  gUnc2->SetFillColor(kBlue);
  gUnc2->SetFillStyle(3001);
  gUnc2->Draw("sameE3");

  can2->Print("plots/"+TString(can2->GetName())+".pdf");
  can2->Print("plots/"+TString(can2->GetName())+".png");

  TFile *outf = TFile::Open("ScaleFactor_"+CUT+"_"+VAR+"_boosted.root","RECREATE");
  fit1->Write();
  fit2->Write();
  outf->Close();
}

TGraphErrors *getUnc(TMatrixDSym COV, TF1 *fit, TH1F *h)
{
  int N = COV.GetNrows();
  float vx[200],vy[200],vex[200],vey[200];
  float dx = (h->GetBinLowEdge(h->GetNbinsX())+h->GetBinWidth(1)-h->GetBinLowEdge(1))/200;
  for(int b=0;b<200;b++) {
    vx[b]  = h->GetBinLowEdge(1)+b*dx;
    vy[b]  = fit->Eval(vx[b]);
    vex[b] = 0.0;
    float sum(0.0);
    for(int i=0;i<N;i++) {
      for(int j=0;j<N;j++) {
        sum += pow(vx[b],i)*pow(vx[b],j)*COV(i,j);
      }
    }
    vey[b] = sqrt(sum);
  }  
  TGraphErrors *g = new TGraphErrors(200,vx,vy,vex,vey);
  return g;
}

void FilterHisto(TH1F *h)
{
  for(int i=0;i<h->GetNbinsX();i++) {
    float y = h->GetBinContent(i+1);
    float e = h->GetBinError(i+1);
    if (y > 0) {
      if ((e/y < 0.01) || (y < 0.4)) {
        h->SetBinContent(i+1,0); 
        h->SetBinError(i+1,1); 
      }
    }  
  }
}












