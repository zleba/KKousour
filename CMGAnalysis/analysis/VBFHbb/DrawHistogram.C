void DrawHistogram(TString HISTO,TString SELECTION,TString XTITLE,float XMIN,float XMAX,int REBIN,bool LOGY,bool PRINT,bool PRINT_PAS,int XDIV,TString UNITS)
{
  gROOT->ForceStyle();
  const int N = 16;
  TFile *inf[N];
  TH1F *h[N],*hPass;
  TString PATH("");
  float LUMI = 0;
  if (SELECTION == "NOM") {
    LUMI = 19785;
  }
  if (SELECTION == "VBF") {
    LUMI = 11640;
  }
  TString fileName[N] = {
    "Histo_sel"+SELECTION+"_data.root",
    "Histo_sel"+SELECTION+"_VBFPowheg125.root",
    "Histo_sel"+SELECTION+"_GFPowheg125.root",
    "Histo_sel"+SELECTION+"_ZJets.root",
    "Histo_sel"+SELECTION+"_WJets.root", 
    "Histo_sel"+SELECTION+"_TTJets.root",
    "Histo_sel"+SELECTION+"_T_t-channel.root",
    "Histo_sel"+SELECTION+"_T_tW-channel.root",
    "Histo_sel"+SELECTION+"_T_s-channel.root",
    "Histo_sel"+SELECTION+"_Tbar_t-channel.root",
    "Histo_sel"+SELECTION+"_Tbar_tW-channel.root", 
    "Histo_sel"+SELECTION+"_Tbar_s-channel.root",
    "Histo_sel"+SELECTION+"_QCD100.root",
    "Histo_sel"+SELECTION+"_QCD250.root",
    "Histo_sel"+SELECTION+"_QCD500.root",
    "Histo_sel"+SELECTION+"_QCD1000.root"
  };
  char name[1000];
  
  float XSEC[N]    = {1,0.911,11.26,650,1.2*1205,245.8,56.4,11.1,3.79,30.7,11.1,1.76,1.036e+7,2.76e+5,8426,204.};
  int FILLCOLOR[N] = {kBlack,kRed,kBlue,kMagenta-10,kYellow-10,kGreen+2,kGreen-8,kGreen-8,kGreen-8,kGreen-8,kGreen-8,kGreen-8,kBlue-6,kBlue-8,kBlue-9,kBlue-10};
  int LINECOLOR[N] = {kBlack,kRed,kBlue,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack,kBlack};
  int LINESIZE[N]  = {1,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1};
  int LINESTYLE[N] = {1,1,2,1,1,1,1,1,1,1,1,1,1,1,1,1};
  THStack *hs = new THStack("bkg","bkg");
  TH1F *hBkg;
  TH1F *hQCD;
  TH1F *hTop;
  
  for(int iFile=0;iFile<N;iFile++) {
    inf[iFile] = TFile::Open(PATH+fileName[iFile]);
    hPass = (TH1F*)inf[iFile]->Get("TriggerPass");
    h[iFile] = (TH1F*)inf[iFile]->Get(HISTO);
    h[iFile]->Sumw2();
    h[iFile]->Rebin(REBIN);
    h[iFile]->SetLineColor(LINECOLOR[iFile]);
    h[iFile]->SetLineWidth(LINESIZE[iFile]);
    h[iFile]->SetLineStyle(LINESTYLE[iFile]);
    if (iFile == 0) {
      hBkg = (TH1F*)h[iFile]->Clone("hBkg");
      hBkg->Reset();
      hQCD = (TH1F*)h[iFile]->Clone("hQCD");
      hQCD->Reset();
      hTop = (TH1F*)h[iFile]->Clone("hTop");
      hTop->Reset(); 
    }
    else {
      float wt = LUMI*XSEC[iFile]/hPass->GetBinContent(1);
      h[iFile]->Scale(wt);
      if (iFile > 2) {
        h[iFile]->SetFillColor(FILLCOLOR[iFile]);
        hBkg->Add(h[iFile]);
      } 
       
      if (iFile > 5 && iFile < 12) {
        hTop->Add(h[iFile]);
      } 
      
      if (iFile > 11) {
        hQCD->Add(h[iFile]);
      }
    }
  }  
  
  float NQCD = hQCD->Integral();
  float NTOT = hBkg->Integral();
  float NDAT = h[0]->Integral();
  float kfactor = (NDAT - NTOT)/NQCD + 1;
  cout<<"k-factor: "<<kfactor<<endl;
  /*
  for(int iFile=6;iFile<N;iFile++) {
    if (iFile > 6) {
      h[iFile]->Scale(kfactor);
    }  
  }
  */
  hQCD->Scale(kfactor);
  hBkg->Scale(1+(kfactor-1)*NQCD/NTOT);
  TH1F *hRatio = (TH1F*)h[0]->Clone("Ratio");
  double KS = h[0]->KolmogorovTest(hBkg);
  //cout<<"KS = "<<KS<<endl;
  
  hRatio->Divide(hBkg);
  
  hTop->SetFillColor(kGreen-8);
  hQCD->SetFillColor(kBlue-8);

  hs->Add(h[4]);
  hs->Add(hTop);
  hs->Add(h[5]);
  hs->Add(h[3]);
  hs->Add(hQCD);
  
  TCanvas *can = new TCanvas("can_"+HISTO,"can_"+HISTO,900,600);
  can->cd(1);
  can->SetBottomMargin(0.3);
  can->SetRightMargin(0.2);
  if (LOGY) gPad->SetLogy();
  TH1F *haux = (TH1F*)hBkg->Clone("aux");
  haux->Reset();
  if (haux->GetBinWidth(1) >= 1) {
    sprintf(name,"Events / %1.0f %s",haux->GetBinWidth(1),UNITS.Data());
  }
  else {
    sprintf(name,"Events / %1.2f %s",haux->GetBinWidth(1),UNITS.Data());
  }  
  haux->GetYaxis()->SetTitle(name);
  
  //haux->GetXaxis()->SetTitle(XTITLE);
  haux->GetXaxis()->SetLabelSize(0.0);
  haux->GetYaxis()->SetNdivisions(505);
  haux->GetXaxis()->SetNdivisions(XDIV);
  haux->GetXaxis()->SetRangeUser(XMIN,XMAX);
  haux->SetMinimum(0.5);
  //haux->SetMaximum(1.2*h[0]->GetBinContent(h[0]->GetMaximumBin()));
  float ymax = 2*h[0]->GetBinContent(h[0]->GetMaximumBin());
  haux->SetMaximum(ymax);
  haux->Draw();
  hs->Draw("same hist");
  h[0]->Draw("same E");
  h[1]->Draw("same hist");
  h[2]->Draw("same hist");
  TLegend *leg = new TLegend(0.81,0.3,0.92,0.7);
  leg->SetHeader("Presel. & Trigger");
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->AddEntry(h[0],"Data","P");
  sprintf(name,"QCD (#times %1.2f)",kfactor);
  leg->AddEntry(hQCD,name,"F");
  leg->AddEntry(h[3],"Z+jets","F");
  leg->AddEntry(h[5],"t#bar{t}","F");
  leg->AddEntry(hTop,"single top","F");
  leg->AddEntry(h[4],"W+jets","F"); 
  leg->AddEntry(h[1],"VBF H(125)#rightarrow b#bar{b}","L");
  leg->AddEntry(h[2],"GF H(125)#rightarrow b#bar{b}","L"); 
  leg->Draw();
  gPad->RedrawAxis();
  TPad* pad = new TPad("pad", "pad", 0., 0., 1., 1.);
  pad->SetTopMargin(0.7);
  pad->SetRightMargin(0.2);
  pad->SetFillColor(0);
  pad->SetFillStyle(0);
  pad->Draw();
  pad->cd(0);
  gPad->SetGridy();
  hRatio->SetMarkerSize(0.7);
  hRatio->SetMinimum(0);
  hRatio->SetMaximum(2);
  hRatio->GetYaxis()->SetTitle("Data / MC");
  if (UNITS == "") {
    hRatio->GetXaxis()->SetTitle(XTITLE);
  }
  else {
    hRatio->GetXaxis()->SetTitle(XTITLE+" ("+UNITS+")");
  }
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTickLength(0.06);
  hRatio->GetYaxis()->SetTitleSize(0.04);
  hRatio->GetXaxis()->SetTitleSize(0.05);
  hRatio->GetYaxis()->SetLabelSize(0.03);
  hRatio->GetXaxis()->SetLabelSize(0.04);
  hRatio->GetXaxis()->SetNdivisions(XDIV);
  hRatio->GetYaxis()->CenterTitle(kTRUE);
  hRatio->GetXaxis()->SetRangeUser(XMIN,XMAX);
  hRatio->Draw();
  
  TPaveText *paveKS = new TPaveText(0.81,0.4,0.92,0.5,"NDC");
  sprintf(name,"KS = %1.1e",KS);
  paveKS->AddText(name);
  paveKS->SetFillColor(0);
  paveKS->SetLineColor(0);
  paveKS->SetBorderSize(0);
  paveKS->SetTextFont(42);
  paveKS->SetTextSize(0.03);
  //paveKS->Draw();
  
  TPaveText *pave = new TPaveText(0.81,0.75,0.92,0.92,"NDC");
  pave->AddText("CMS");
  pave->AddText("#sqrt{s} = 8 TeV");
  sprintf(name,"L = %1.1f fb^{-1}",LUMI/1000);
  pave->AddText(name); 
  pave->SetFillColor(0);
  pave->SetLineColor(0);
  pave->SetBorderSize(0);
  pave->SetTextFont(62);
  pave->SetTextSize(0.04);
  
  if (PRINT_PAS) {
    pave->Draw(); 
  }  
  if (PRINT) {
    can->Print(TString(can->GetTitle())+".pdf");
  }
}
