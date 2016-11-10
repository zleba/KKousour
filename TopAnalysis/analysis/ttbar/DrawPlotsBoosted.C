#include "CMS_lumi.C"
using namespace RooFit;
void DrawPlotsBoosted(TString CAT, TString CUT)
{
  gROOT->ForceStyle();
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  const int NVAR  = 12;
  
  TString VAR[NVAR]    = {"mTop","jetPt","jetEta","ptJJ","mJJ","mva",
                          "mW","jetMassSub1","jetPtSub0","jetPtSub1","jetTau32","jetTau31"};
  TString UNITS[NVAR]  = {"GeV","GeV","","GeV","GeV","","GeV","GeV","GeV","GeV","",""}; 
  TString XTITLE[NVAR] = {"Jet mass","Jet p_{T}","Jet #eta",
                          "p_{T,tt}","m_{tt}","Fisher discriminant","Jet lead. subjet mass","Jet second subjet mass",
                          "Jet lead. subjet p_{T}","Jet second subjet p_{T}","Jet #tau_{3}/#tau_{2}","Jet #tau_{3}/#tau_{1}"};
  float XMIN[NVAR]     = {70,350,-2.5,0,800,-1,0,0,200,0,0,0};
  float XMAX[NVAR]     = {300,1300,2.5,800,3000,1.5,200,100,1000,500,1,1};
  float NORM[NVAR]     = {1.0,2.0,2.0,1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,2.0};
  int REBIN[NVAR]      = {5,25,2,20,50,5,5,2,20,20,2,2};
  int NDIV[NVAR]       = {505,505,510,505,505,510,510,510,510,510,510,510};
  bool LOG[NVAR]       = {false,false,false,true,true,false,false,false,false,false,false};

  TH1F *h[NVAR][3];
  TH1F *hRatio[NVAR],*hAll[NVAR];
  TFile *inf[2];
  inf[0] = TFile::Open("Histo_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  //inf[0] = TFile::Open("Histo_TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root");
  inf[1] = TFile::Open("Histo_JetHT.root");
  RooWorkspace *w; 
  TCanvas *can[NVAR];
  THStack *hs[NVAR];

  TFile *fTemplates = TFile::Open("templates_"+CUT+"_"+CAT+"_workspace.root");

  TFile *fYields = TFile::Open("FittedYields_"+CUT+"_"+CAT+".root");
  TH1F  *hSig    = (TH1F*)fYields->Get("YieldSigData");
  TH1F  *hBkg    = (TH1F*)fYields->Get("YieldBkgData");
  TH1F  *hCor    = (TH1F*)fYields->Get("YieldCorrelation");

  float NSig = hSig->GetBinContent(1);
  float NBkg = hBkg->GetBinContent(1);
  float ESig = hSig->GetBinError(1);
  float EBkg = hBkg->GetBinError(1);

  w = (RooWorkspace*)fTemplates->Get("w");
  RooRealVar *x = (RooRealVar*)w->var("mTop");
  RooAbsPdf *pdf_bkg = (RooAbsPdf*)w->pdf("qcdCor_pdf");
  RooAbsPdf *pdf_sig = (RooAbsPdf*)w->pdf("ttbar_pdf_Nominal");

  RooAbsReal *fbkg;
  RooAbsReal *fsig;
  ///*
  x->setRange("signal",140,200);
  fbkg = pdf_bkg->createIntegral(*x,RooFit::NormSet(*x),Range("signal"));
  fsig = pdf_sig->createIntegral(*x,RooFit::NormSet(*x),Range("signal"));
  cout<<"Signal: "<<fsig->getVal()*NSig<<", Bkg: "<<fbkg->getVal()*NBkg<<endl;
  cout<<"Signal fraction: "<<fsig->getVal()<<", Bkg fraction: "<<fbkg->getVal()<<endl;
  //*/
  /*
  for(int i=0;i<6;i++) {
    for(int j=0;j<6;j++) {
      float m1 = 130+i*5;
      float m2 = 180+j*5;
      x->setRange("signal",m1,m2);
      fbkg = pdf_bkg->createIntegral(*x,RooFit::NormSet(*x),Range("signal"));
      fsig = pdf_sig->createIntegral(*x,RooFit::NormSet(*x),Range("signal"));
      float sovb = fsig->getVal()*NSig/(fbkg->getVal()*NBkg);
      cout<<m1<<" "<<m2<<" "<<fsig->getVal()<<" "<<fbkg->getVal()<<" "<<sovb<<endl;
    }
  } 
  */   
  for(int k=0;k<NVAR;k++) {
    h[k][0] = (TH1F*)inf[0]->Get("boosted/hWt_"+VAR[k]+"_"+CUT+"_"+CAT); //MC
    h[k][1] = (TH1F*)inf[1]->Get("boosted/h_"+VAR[k]+"_"+CUT+"_0btag");// data
    h[k][2] = (TH1F*)inf[1]->Get("boosted/h_"+VAR[k]+"_"+CUT+"_"+CAT);// data
    for(int j=0;j<3;j++) {
      h[k][j]->Rebin(REBIN[k]);
      if (j == 2) {
        cout<<VAR[k]<<", Data events: "<<h[k][j]->GetEntries()<<endl;
      }
    }
    if (VAR[k].CompareTo("mTop") == 0) {
      h[k][0]->Scale(NORM[k]*NSig/h[k][0]->Integral());
      h[k][1]->Scale(NORM[k]*NBkg/h[k][1]->Integral());
    }
    else {
      h[k][0]->Scale(NORM[k]*fsig->getVal()*NSig/h[k][0]->Integral());
      h[k][1]->Scale(NORM[k]*fbkg->getVal()*NBkg/h[k][1]->Integral()); 
    }
    h[k][0]->SetFillColor(kRed-10);  
    h[k][1]->SetFillColor(kGray);
    hs[k] = new THStack("hs_"+VAR[k],"hs_"+VAR[k]);
    hs[k]->Add(h[k][1]);
    hs[k]->Add(h[k][0]);
    
    hAll[k] = (TH1F*)h[k][0]->Clone("All_"+VAR[k]);
    hAll[k]->Add(h[k][1]);
    hAll[k]->SetFillColor(kBlack);
    hAll[k]->SetFillStyle(3004); 
    hAll[k]->SetMarkerSize(0);

    //--- set the errors of hAll including the errors on the normalization
    for(int i=0;i<hAll[k]->GetNbinsX();i++) {
      float eStat = hAll[k]->GetBinError(i+1);
      float rho   = hCor->GetBinContent(1);
      float y     = hAll[k]->GetBinContent(i+1);
      float S     = h[k][0]->GetBinContent(i+1); 
      float B     = h[k][1]->GetBinContent(i+1); 
      float eS    = ESig/NSig;
      float eB    = EBkg/NBkg;
      if (VAR[k].CompareTo("mTop") != 0) {
        eS = eS/fsig->getVal();
        eB = eB/fbkg->getVal();
      }
      float e = sqrt(pow(S*eS,2)+pow(B*eB,2)+2*rho*S*B*eS*eB);
      hAll[k]->SetBinError(i+1,sqrt(eStat*eStat+e*e));
    }

    TH1F *hAux = (TH1F*)h[k][2]->Clone("aux_"+VAR[k]);
    hAux->Reset();
    can[k] = new TCanvas("PostFit_"+CUT+"_"+CAT+"_"+VAR[k],"PostFit_"+CUT+"_"+CAT+"_"+VAR[k],600,600);
    can[k]->cd(1)->SetBottomMargin(0.25);
    if (LOG[k]) {
      gPad->SetLogy();
    }
    hAux->Draw();
    hAux->SetMinimum(0.5);
    hAux->SetMaximum(1.3*TMath::Max(h[k][2]->GetBinContent(h[k][2]->GetMaximumBin()),h[k][0]->GetBinContent(h[k][0]->GetMaximumBin())));
    hAux->GetYaxis()->SetTitleOffset(1.4);
    hAux->GetYaxis()->SetTitle("Events / "+TString::Format("%1.2g",h[k][2]->GetBinWidth(2))+" "+UNITS[k]);
    hAux->GetXaxis()->SetTitle("");
    hAux->GetXaxis()->SetLabelSize(0.0);
    hAux->GetXaxis()->SetRangeUser(XMIN[k],XMAX[k]);
    hAux->GetXaxis()->SetNdivisions(NDIV[k]);
    hAux->GetYaxis()->SetNdivisions(510);
    hs[k]->Draw("same hist");
    hAll[k]->Draw("same E2");
    h[k][2]->Draw("sameE");
    gPad->RedrawAxis();

    if (VAR[k] == "mTop") {
      gPad->Update();

      TLine *ln1 = new TLine(140,gPad->GetFrame()->GetY1(),140,gPad->GetFrame()->GetY2());
      ln1->SetLineColor(kBlack);
      ln1->SetLineWidth(2);
      ln1->SetLineStyle(2);
      ln1->Draw();

      TLine *ln2 = new TLine(200,gPad->GetFrame()->GetY1(),200,gPad->GetFrame()->GetY2());
      ln2->SetLineColor(kBlack);
      ln2->SetLineWidth(2);
      ln2->SetLineStyle(2);
      ln2->Draw();
    }

    TLegend *leg = new TLegend(0.7,0.73,0.9,0.92);
    leg->AddEntry(hAux,"Data","P");
    leg->AddEntry(h[k][0],"Signal","F");
    leg->AddEntry(h[k][1],"QCD","F");
    leg->AddEntry(hAll[k],"Fit Unc.","F");
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.045);
    leg->Draw();
      
    TPaveText *paveBIN = new TPaveText(0.14,0.98,0.6,0.99,"NDC");
    paveBIN->SetFillColor(0);
    paveBIN->SetBorderSize(0);
    paveBIN->SetTextFont(42);
    paveBIN->SetTextSize(0.04);
    paveBIN->SetTextAlign(13);

    TPaveText *paveBINin = new TPaveText(0.165,0.91,0.62,0.92,"NDC");
    paveBINin->SetFillColor(0);
    paveBINin->SetBorderSize(0);
    paveBINin->SetTextFont(42);
    paveBINin->SetTextSize(0.042);
    paveBINin->SetTextAlign(13);
       
    TH1F *hAllUncUp = (TH1F*)hAll[k]->Clone("uncUp");
    TH1F *hAllUncLo = (TH1F*)hAll[k]->Clone("uncLo");
    TH1F *hAllAux   = (TH1F*)hAll[k]->Clone("aux");
    hRatio[k] = (TH1F*)h[k][2]->Clone("Ratio_"+VAR[k]);
    for(int bin=0;bin<hAllUncUp->GetNbinsX();bin++) {
      float e = 0.0;
      float r = 0.0;
      float er = 0.0;
      if (hAll[k]->GetBinContent(bin+1) != 0) {
        e = hAll[k]->GetBinError(bin+1)/hAll[k]->GetBinContent(bin+1);
      }
      hAllUncUp->SetBinContent(bin+1,e);
      hAllUncLo->SetBinContent(bin+1,-e);
      hAllAux->SetBinError(bin+1,0);
    }
      
    hRatio[k]->Add(hAllAux,-1);
    hRatio[k]->Divide(hAllAux);
    
    TPad *pad = new TPad("pad","pad",0.,0.,1.,1.);
    pad->SetTopMargin(0.77);
    pad->SetFillColor(0);
    pad->SetFillStyle(0);
    pad->Draw();
    pad->cd(0);
    pad->SetGridy();
    hAllUncUp->SetMinimum(-0.9);
    hAllUncUp->SetMaximum(0.9);
    hAllUncUp->GetYaxis()->SetNdivisions(505);
    hAllUncUp->GetXaxis()->SetTitleOffset(0.95);
    hAllUncUp->GetYaxis()->SetTitleOffset(1.5);
    hAllUncUp->GetYaxis()->SetTickLength(0.06);
    hAllUncUp->GetYaxis()->SetTitleSize(0.03);
    hAllUncUp->GetYaxis()->SetLabelSize(0.03);
    hAllUncUp->GetYaxis()->SetTitle("Data/(S+B)-1");
    hAllUncUp->GetXaxis()->SetNdivisions(NDIV[k]);
    if (UNITS[k] == "") {
      hAllUncUp->GetXaxis()->SetTitle(XTITLE[k]);
    }
    else {
      hAllUncUp->GetXaxis()->SetTitle(XTITLE[k]+" ("+UNITS[k]+")");
    }
    hAllUncUp->GetXaxis()->SetRangeUser(XMIN[k],XMAX[k]);
    hAllUncUp->Draw("HIST"); 
    hAllUncLo->Draw("HIST SAME"); 
    hRatio[k]->Draw("SAME");

    //CMS_lumi(can[k],4,0);

    //can[k]->Print("plots/"+TString(can[k]->GetName())+".pdf");     
    //can[k]->Print("plots/"+TString(can[k]->GetName())+".png");
  }
}
