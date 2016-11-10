void FitData2D(TString CAT,TString CUT)
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  
  TFile *inf = TFile::Open("Histo_JetHT.root");
  TH2F *h = (TH2F*)inf->Get("boosted/h_mTopVsmW_"+CUT+"_"+CAT);
  h->Draw();
  // -----------------------------------------
  TFile *fTemplatesW = TFile::Open("templates_W_"+CUT+"_"+CAT+"_workspace.root");
  TFile *fTemplatesT = TFile::Open("templates_"+CUT+"_"+CAT+"_workspace.root");
  
  RooWorkspace *wTemplatesW = (RooWorkspace*)fTemplatesW->Get("w");
  RooWorkspace *wTemplatesT = (RooWorkspace*)fTemplatesT->Get("w");
  RooRealVar *xW = (RooRealVar*)wTemplatesW->var("mW");
  RooRealVar *xT = (RooRealVar*)wTemplatesT->var("mT");
  
  RooAbsPdf *pdf_signal_W = (RooAbsPdf*)wTemplatesW->pdf("ttbar_pdf");
  RooAbsPdf *pdf_bkg_W    = (RooAbsPdf*)wTemplatesW->pdf("qcd_pdf");
  
  RooAbsPdf *pdf_signal_T = (RooAbsPdf*)wTemplatesT->pdf("ttbar_pdf_Nominal");
  RooAbsPdf *pdf_bkg_T    = (RooAbsPdf*)wTemplatesT->pdf("qcd_pdf");
  
  RooProdPdf pdf_signal("signal","signal",RooArgSet(*pdf_signal_W,*pdf_signal_T));
  
  RooProdPdf pdf_bkg("bkg","bkg",RooArgSet(*pdf_bkg_W,*pdf_bkg_T));
  

  RooRealVar *nFitBkg = new RooRealVar("nFitBkg","nFitBkg",2000,0,2000000);
  RooRealVar *nFitSig = new RooRealVar("nFitSig","nFitSig",2000,0,50000);
  
  RooAddPdf *model = new RooAddPdf("model","model",RooArgList(pdf_signal,pdf_bkg),RooArgList(*nFitSig,*nFitBkg));

  RooDataHist *roohist_data = new RooDataHist("roohist_data","roohist_data",RooArgList(*xW,*xT),h);
  
  RooFitResult *res = model->fitTo(*roohist_data,RooFit::Save(),RooFit::Extended(kTRUE));
  res->Print();

  RooPlot *frameW = xW->frame();
  roohist_data->plotOn(frameW);
  model->plotOn(frameW);
  /*
  RooPlot *frameT = xT->frame();
  roohist_data->plotOn(frameT);
  model->plotOn(frameT);
  */
  TCanvas *canW = new TCanvas("cW","cW",900,600);
  frameW->Draw();
  /*
  TCanvas *canT = new TCanvas("cT","cT",900,600);
  frameT->Draw();
  */
}








