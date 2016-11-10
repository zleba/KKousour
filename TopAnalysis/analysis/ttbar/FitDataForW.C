void FitDataForW(TString CAT,TString CUT,int REBIN)
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  
  TFile *inf = TFile::Open("Histo_JetHT.root");
  TH1F *h = (TH1F*)inf->Get("boosted/h_mW_"+CUT+"_"+CAT);
  h->Rebin(REBIN);
  // -----------------------------------------
   
  TFile *fTemplates = TFile::Open("templates_W_"+CUT+"_"+CAT+"_workspace.root");
  RooWorkspace *wTemplates = (RooWorkspace*)fTemplates->Get("w");
  RooRealVar *x = (RooRealVar*)wTemplates->var("mW");
    
  RooDataHist *roohist_data = new RooDataHist("roohist_data","roohist_data",RooArgList(*x),h);
  RooRealVar *nFitBkg = new RooRealVar("nFitBkg","nFitBkg",2000,0,2000000);
  RooRealVar *nFitSig = new RooRealVar("nFitSig","nFitSig",2000,0,50000); 

  RooPlot *frame = x->frame();
  roohist_data->plotOn(frame);
 
  RooRealVar *meanMC  = (RooRealVar*)wTemplates->var("ttbar_mean");
  RooRealVar *sigmaMC = (RooRealVar*)wTemplates->var("ttbar_sigma");
  meanMC->setConstant(false);
  sigmaMC->setConstant(false);
  float mPreFit  = meanMC->getVal();
  float sPreFit  = sigmaMC->getVal();
  float emPreFit = meanMC->getError();
  float esPreFit = sigmaMC->getError();

  RooAbsPdf *pdf_signal = (RooAbsPdf*)wTemplates->pdf("ttbar_pdf");
  RooAbsPdf *pdf_bkg = (RooAbsPdf*)wTemplates->pdf("qcd_pdf"); 
  //---- QCD correction factor ---------------------------
  RooRealVar kBkg("bkgCor","bkgCor",0.5,-1,1);
  RooFormulaVar eff("eff","1+@0*@1",RooArgList(*x,kBkg));
  //---- corrected QCD -----------------------------------
  RooEffProd qcdEff("qcdEff","qcdEff",*pdf_bkg,eff);
  RooAddPdf *model = new RooAddPdf("model","model",RooArgList(*pdf_signal,qcdEff),RooArgList(*nFitSig,*nFitBkg)); 
 //kBkg.setConstant(true);
  RooFitResult *res = model->fitTo(*roohist_data,RooFit::Save(),RooFit::Extended(kTRUE));
  res->Print();

  model->plotOn(frame);
  float chi2 = frame->chiSquare(4);//chi2 with 4 free parameters

  cout<<"chi2   = "<<chi2<<endl;
  cout<<"kMean  = "<<meanMC->getVal()/mPreFit<<" +/- "<<sqrt(pow(meanMC->getError(),2)+pow(emPreFit,2))/mPreFit<<endl;
  cout<<"kSigma = "<<sigmaMC->getVal()/sPreFit<<" +/- "<<sqrt(pow(sigmaMC->getError(),2)+pow(esPreFit,2))/sPreFit<<endl;

  RooHist *hpull = frame->pullHist();
  model->plotOn(frame,RooFit::Components("qcd_pdf"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(kDashed));
  model->plotOn(frame,RooFit::Components("ttbar_pdf"),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(1));

  TCanvas *can = new TCanvas("Fitted_"+CUT+"_"+CAT+"_mW","Fitted_"+CUT+"_"+CAT+"_mW",900,600);
  can->cd(1)->SetBottomMargin(0.3);
  frame->SetMinimum(0.1);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetXaxis()->SetTitle("");
  frame->GetXaxis()->SetLabelSize(0.0);
  frame->Draw(); 

  RooHist  *hist_data = (RooHist*)frame->findObject("h_roohist_data",RooHist::Class());
  RooCurve *curve_all = (RooCurve*)frame->findObject("model_Norm[mW]",RooCurve::Class());
  RooCurve *curve_sig = (RooCurve*)frame->findObject("model_Norm[mW]_Comp[ttbar_pdf]",RooCurve::Class());
  RooCurve *curve_bkg = (RooCurve*)frame->findObject("model_Norm[mW]_Comp[qcd_pdf]",RooCurve::Class());

  TLegend *leg = new TLegend(0.75,0.7,0.92,0.9);
  leg->AddEntry(hist_data,"data","P");
  leg->AddEntry(curve_all,"signal + qcd","L");
  leg->AddEntry(curve_sig,"signal","L");
  leg->AddEntry(curve_bkg,"qcd","L");
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
  frame2->GetXaxis()->SetTitle("m_{W} (GeV)");
  frame2->Draw();   

  can->Print("plots/"+TString(can->GetName())+".pdf");
  can->Print("plots/"+TString(can->GetName())+".png");
}








