TH1F* NormHisto(TH1F *h) 
{
  TH1F *hNew = (TH1F*)h->Clone(TString(h->GetName())+"_Norm");
  float norm = hNew->GetBinContent(1);
  for(int k=0;k<hNew->GetNbinsX();k++) {
    hNew->SetBinContent(k+1,hNew->GetBinContent(k+1)/norm);
    hNew->SetBinError(k+1,0);
  }  
  return hNew;
}
void FitData(TString CAT,TString CUT,int REBIN)
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  
  TFile *inf = TFile::Open("Histo_JetHT.root");
  TH1F *h = (TH1F*)inf->Get("boosted/h_mTop_"+CUT+"_"+CAT);
  h->Rebin(REBIN);
  // -----------------------------------------
  const float LUMI = 15900;
  const int NMC = 1;
  TString ALIAS[NMC] = {"Nominal"
                        /*
                        "BkgSFUp","BkgSFDown","NoBkgSF","NoWeight",
                        "BtagUp","BtagDown","TrigUp","TrigDown","JER","JERUp","JERDown","JESUp","JESDown",
                        "aMC@NLO","Madgraph","ScaleDown","ScaleUp","mtop1665","mtop1695","mtop1715","mtop1735","mtop1755","mtop1785",
                        "mpiOFF","herwigpp"
                        */
  }; 
  TH1F *hYieldBkg    = new TH1F("YieldBkgData","YieldBkgData",NMC,0,NMC);
  TH1F *hYieldSig    = new TH1F("YieldSigData","YieldSigData",NMC,0,NMC);
  TH1F *hYieldCor    = new TH1F("YieldCorrelation","YieldCorrelation",NMC,0,NMC);
  TH1F *hYieldSigExp = new TH1F("YieldSigExp","YieldSigExp",NMC,0,NMC);
  TH1F *hAccSigExp   = new TH1F("AccSigExp","AccSigExp",NMC,0,NMC); 
  TH1F *hSigFrac     = new TH1F("SignalFraction","SignalFraction",NMC,0,NMC);
  TH1F *hBkgFrac     = new TH1F("BkgFraction","BkgFraction",NMC,0,NMC);  
  TH1F *hChi2        = new TH1F("Chi2","Chi2",NMC,0,NMC);
  TH1F *hXsecFid     = new TH1F("XsecFiducial","XsecFiducial",NMC,0,NMC);
  TH1F *hXsecFidExp  = new TH1F("XsecFiducialExp","XsecFiducialExp",NMC,0,NMC);
  TH1F *hXsecExtr    = new TH1F("XsecExtrapolated","XsecExtrapolated",NMC,0,NMC);

  TFile *fTemplates = TFile::Open("templates_"+CUT+"_"+CAT+"_workspace.root");
  RooWorkspace *wTemplates = (RooWorkspace*)fTemplates->Get("w");
  RooRealVar *x = (RooRealVar*)wTemplates->var("mTop");
  x->setRange("signal",140,200);
  RooRealVar *kMassScale = (RooRealVar*)wTemplates->var("kMassScale");
  RooRealVar *kMassResol = (RooRealVar*)wTemplates->var("kMassResol");
  kMassScale->setConstant(true);
  kMassResol->setConstant(true);
  //kMassScale->setVal(0.99);
  //kMassResol->setVal(1.0);
    
  RooDataHist *roohist_data = new RooDataHist("roohist_data","roohist_data",RooArgList(*x),h);
  RooAbsPdf *pdf_bkg = (RooAbsPdf*)wTemplates->pdf("bkg_pdf");
  RooAbsPdf *pdf_qcd = (RooAbsPdf*)wTemplates->pdf("qcd_pdf");
  cout<<pdf_bkg->GetName()<<endl;
  cout<<pdf_qcd->GetName()<<endl;
  //---- QCD correction factor ---------------------------
  RooRealVar kQCD("kQCD","kQCD",0,-1,5);
  kQCD.setConstant(false);
  RooFormulaVar qcdCor("qcdCor","1+@0*@1",RooArgList(*x,kQCD));
  //---- corrected QCD -----------------------------------
  RooEffProd pdf_qcdCor("qcdCor_pdf","qcdCor_pdf",*pdf_qcd,qcdCor);
  RooAbsReal *fqcd = pdf_qcdCor.createIntegral(*x,RooFit::NormSet(*x),RooFit::Range("signal"));
  RooRealVar *nFitBkg = new RooRealVar("nFitBkg","nFitBkg",200,0,20000);
  RooRealVar *nFitQCD = new RooRealVar("nFitQCD","nFitQCD",2000,0,200000);  
  
  TCanvas *can[NMC];
  
  for(int k=0;k<NMC;k++) {
    cout<<"========= "+ALIAS[k]<<" ============="<<endl;
    TString TTBARNAME = ALIAS[k];

    RooRealVar *exp_sig = (RooRealVar*)wTemplates->var("YieldTT_"+TTBARNAME);
    RooRealVar *exp_acc = (RooRealVar*)wTemplates->var("AccTT_"+TTBARNAME);

    RooPlot *frame = x->frame();
    roohist_data->plotOn(frame);

    RooRealVar *nFitSig = new RooRealVar("nFitSig_"+ALIAS[k],"nFitSig_"+ALIAS[k],2000,500,10000);
    RooAbsPdf *pdf_signal = (RooAbsPdf*)wTemplates->pdf("ttbar_pdf_"+TTBARNAME);

    RooAbsReal *fsig = pdf_signal->createIntegral(*x,RooFit::NormSet(*x),RooFit::Range("signal"));
    
    RooAddPdf *model = new RooAddPdf("model","model",RooArgList(*pdf_signal,pdf_qcdCor,*pdf_bkg),RooArgList(*nFitSig,*nFitQCD,*nFitBkg)); 
    model->Print();
    RooFitResult *res = model->fitTo(*roohist_data,RooFit::Save(),RooFit::Extended(kTRUE));
    res->Print();
    
    model->plotOn(frame);
    
    float chi2 = frame->chiSquare(2);
    float xsec = fsig->getVal()*nFitSig->getVal()/LUMI;
    float xsec_err = fsig->getVal()*nFitSig->getError()/LUMI;
    float xsec_exp = fsig->getVal()*exp_sig->getVal()/LUMI;
    float xsec_exp_err = fsig->getVal()*exp_sig->getError()/LUMI;

    cout<<"Signal          = "<<nFitSig->getVal()<<" +/- "<<nFitSig->getError()<<endl;
    cout<<"QCD             = "<<nFitQCD->getVal()<<" +/- "<<nFitQCD->getError()<<endl;
    cout<<"Bkg             = "<<nFitBkg->getVal()<<" +/- "<<nFitBkg->getError()<<endl;
    cout<<"Exp. Signal     = "<<exp_sig->getVal()<<" +/- "<<exp_sig->getError()<<endl;
    cout<<"r               = "<<nFitSig->getVal()/exp_sig->getVal()<<endl;
    cout<<"chi2/ndof       = "<<chi2<<endl;
    cout<<"Signal fraction = "<<fsig->getVal()<<endl;
    cout<<"QCD fraction    = "<<fqcd->getVal()<<endl;
    cout<<"correlation     = "<<res->correlation(*nFitSig,*nFitQCD)<<endl;
    
    hYieldSig->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hYieldSig->SetBinContent(k+1,nFitSig->getVal());
    hYieldSig->SetBinError(k+1,nFitSig->getError());
    hYieldBkg->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hYieldBkg->SetBinContent(k+1,nFitQCD->getVal());
    hYieldBkg->SetBinError(k+1,nFitQCD->getError());
    hYieldCor->SetBinContent(k+1,res->correlation(*nFitSig,*nFitQCD));
    
    hYieldSigExp->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hYieldSigExp->SetBinContent(k+1,exp_sig->getVal());
    hYieldSigExp->SetBinError(k+1,exp_sig->getError());

    hAccSigExp->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hAccSigExp->SetBinContent(k+1,exp_acc->getVal());
    hAccSigExp->SetBinError(k+1,exp_acc->getError());

    hXsecFid->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hXsecFid->SetBinContent(k+1,xsec);
    hXsecFid->SetBinError(k+1,xsec_err);

    hXsecFidExp->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hXsecFidExp->SetBinContent(k+1,xsec_exp);
    hXsecFidExp->SetBinError(k+1,xsec_exp_err);
    
    hXsecExtr->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hXsecExtr->SetBinContent(k+1,xsec/(fsig->getVal()*exp_acc->getVal()));
    hXsecExtr->SetBinError(k+1,(xsec/(fsig->getVal()*exp_acc->getVal()))*sqrt(pow(xsec_err/xsec,2)+pow(exp_acc->getError()/exp_acc->getVal(),2)));
    
    hChi2->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hChi2->SetBinContent(k+1,chi2);
    
    hSigFrac->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hSigFrac->SetBinContent(k+1,fsig->getVal());
    
    hBkgFrac->GetXaxis()->SetBinLabel(k+1,ALIAS[k]);
    hBkgFrac->SetBinContent(k+1,fqcd->getVal());
    
    RooHist *hpull = frame->pullHist();
    //model->plotOn(frame,RooFit::VisualizeError(*res,1,kFALSE),RooFit::FillColor(kGray),RooFit::MoveToBack());
    model->plotOn(frame,RooFit::Components("qcdCor_pdf"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
    model->plotOn(frame,RooFit::Components("bkg_pdf"),RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::LineStyle(4));
    model->plotOn(frame,RooFit::Components("ttbar_pdf_"+TTBARNAME),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(1));

    can[k] = new TCanvas("Fitted_"+CUT+"_"+CAT+"_mTop_"+ALIAS[k],"Fitted_"+CUT+"_"+CAT+"_mTop_"+ALIAS[k],900,600);
    can[k]->cd(1)->SetBottomMargin(0.3);
    frame->SetMinimum(0.1);
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetXaxis()->SetTitle("");
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->Draw(); 

    RooHist  *hist_data = (RooHist*)frame->findObject("h_roohist_data",RooHist::Class());
    RooCurve *curve_all = (RooCurve*)frame->findObject("model_Norm[mTop]",RooCurve::Class());
    RooCurve *curve_sig = (RooCurve*)frame->findObject("model_Norm[mTop]_Comp[ttbar_pdf_"+TTBARNAME+"]",RooCurve::Class());
    RooCurve *curve_qcd = (RooCurve*)frame->findObject("model_Norm[mTop]_Comp[qcdCor_pdf]",RooCurve::Class());
    RooCurve *curve_bkg = (RooCurve*)frame->findObject("model_Norm[mTop]_Comp[bkg_pdf]",RooCurve::Class());

    TLegend *leg = new TLegend(0.75,0.7,0.92,0.9);
    leg->SetHeader(ALIAS[k]);
    leg->AddEntry(hist_data,"data","P");
    leg->AddEntry(curve_all,"total","L");
    leg->AddEntry(curve_sig,"ttbar","L");
    leg->AddEntry(curve_qcd,"qcd","L");
    leg->AddEntry(curve_bkg,"other bkg","L"); 
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.04);
    leg->Draw();
        
    RooPlot *frame2 = x->frame();
    frame2->addPlotable(hpull,"p");
  
    TPad *pad = new TPad("pad","pad",0.,0.,1.,1.);
    pad->SetTopMargin(0.7);
    //pad->SetGridy();
    pad->SetFillColor(0);
    pad->SetFillStyle(0);
    pad->Draw();
    pad->cd(0);
    pad->SetGridy();
    frame2->SetMinimum(-5);
    frame2->SetMaximum(5);
    frame2->GetYaxis()->SetNdivisions(505);
    frame2->GetXaxis()->SetTitleOffset(0.9);
    frame2->GetYaxis()->SetTitleOffset(0.8);
    frame2->GetYaxis()->SetTickLength(0.06);
    frame2->GetYaxis()->SetTitleSize(0.05);
    frame2->GetYaxis()->SetTitleSize(0.03);
    frame2->GetYaxis()->SetLabelSize(0.03);
    frame2->GetYaxis()->SetTitle("(Data-Fit)/Error");
    frame2->GetXaxis()->SetTitle("m_{t} (GeV)");
    frame2->Draw();   

    can[k]->Print("plots/"+TString(can[k]->GetName())+".pdf");
    can[k]->Print("plots/"+TString(can[k]->GetName())+".png");
    /*
    RooAbsReal *nll = model->createNLL(*roohist_data,RooFit::NumCPU(2));
    RooMinuit(*nll).migrad();
    RooAbsReal *pll_sig = nll->createProfile(*nFitSig);
    //RooAbsReal *pll_bkg = nll->createProfile(RooArgSet(*nFitSig,kBkg));
    RooPlot *frame_nll = nFitSig->frame(RooFit::Bins(10),RooFit::Range(nFitSig->getVal()-2*nFitSig->getError(),nFitSig->getVal()+2*nFitSig->getError()));
    pll_sig->plotOn(frame_nll,RooFit::ShiftToZero());
    //pll_bkg->plotOn(frame_nll,RooFit::LineStyle(3),RooFit::ShiftToZero(),RooFit::ShiftToZero());
    nll->plotOn(frame_nll,RooFit::LineStyle(kDashed),RooFit::ShiftToZero());
    TCanvas *c = new TCanvas("c","c",900,600);
    frame_nll->SetMinimum(0.0);
    frame_nll->SetMaximum(3.0); 
    frame_nll->Draw();
    */
  }
  TFile *outf = TFile::Open("FittedYields_"+CUT+"_"+CAT+".root","RECREATE");
  outf->cd();
  hYieldSig->Write();
  hYieldBkg->Write(); 
  hYieldCor->Write();
  hYieldSigExp->Write();
  hAccSigExp->Write();
  hChi2->Write();
  hXsecFid->Write();
  hXsecFidExp->Write();
  hXsecExtr->Write();
  hSigFrac->Write();
  hBkgFrac->Write();

  TH1F *hYieldSigNorm    = (TH1F*)NormHisto(hYieldSig);
  TH1F *hYieldBkgNorm    = (TH1F*)NormHisto(hYieldBkg);
  TH1F *hYieldSigExpNorm = (TH1F*)NormHisto(hYieldSigExp);
  TH1F *hAccSigExpNorm   = (TH1F*)NormHisto(hAccSigExp);
  TH1F *hChi2Norm        = (TH1F*)NormHisto(hChi2);
  TH1F *hXsecFidNorm     = (TH1F*)NormHisto(hXsecFid);
  TH1F *hXsecFidExpNorm  = (TH1F*)NormHisto(hXsecFidExp); 
  TH1F *hXsecExtrNorm    = (TH1F*)NormHisto(hXsecExtr);
  TH1F *hSigFracNorm     = (TH1F*)NormHisto(hSigFrac);
  TH1F *hBkgFracNorm     = (TH1F*)NormHisto(hBkgFrac);

  hYieldSigNorm->Write();
  hYieldBkgNorm->Write(); 
  hYieldSigExpNorm->Write();
  hAccSigExpNorm->Write();
  hChi2Norm->Write();
  hXsecFidNorm->Write();
  hXsecFidExpNorm->Write();
  hXsecExtrNorm->Write();
  hSigFracNorm->Write();
  hBkgFracNorm->Write();
}








