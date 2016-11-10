void FitBias(TString CAT,TString CUT,float SIG,float BKG,int NTOYS)
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);
  
  // -----------------------------------------
  TFile *fTemplates = TFile::Open("templates_"+CUT+"_"+CAT+"_workspace.root");
  RooWorkspace *wTemplates = (RooWorkspace*)fTemplates->Get("w");
  RooRealVar *x            = (RooRealVar*)wTemplates->var("mTop");
  RooAbsPdf *pdf_signal    = (RooAbsPdf*)wTemplates->pdf("ttbar_pdf_Nominal");
  RooAbsPdf *pdf_bkg       = (RooAbsPdf*)wTemplates->pdf("qcdCor_pdf"); 
  TRandom *rnd = new TRandom();
  rnd->SetSeed(0);
  x->setBins(250);   
  RooPlot *frame;

  TFile *outf;

  if (NTOYS > 1) { 
    outf = TFile::Open("FitBiasToys_"+CUT+"_"+CAT+".root","RECREATE");
  }

  float nSigInj,nBkgInj,nSigFit,nBkgFit,eSigFit,eBkgFit,nll;

  TTree *tr = new TTree("toys","toys");
  
  tr->Branch("nSigInj",&nSigInj,"nSigInj/F");
  tr->Branch("nSigFit",&nSigFit,"nSigFit/F");
  tr->Branch("nBkgInj",&nBkgInj,"nBkgInj/F");
  tr->Branch("nBkgFit",&nBkgFit,"nBkgFit/F");
  tr->Branch("eSigFit",&eSigFit,"eSigFit/F");
  tr->Branch("eBkgFit",&eBkgFit,"eBkgFit/F");
  tr->Branch("nll"    ,&nll    ,"nll/F");

  for(int itoy=0;itoy<NTOYS;itoy++) {
    // generate pseudodataset
    nSigInj = rnd->Poisson(SIG);
    nBkgInj = rnd->Poisson(BKG);
    RooRealVar *nSig = new RooRealVar("nSig","nSig",nSigInj);
    RooRealVar *nBkg = new RooRealVar("nBkg","nBkg",nBkgInj);
    RooAddPdf *model = new RooAddPdf("model","model",RooArgList(*pdf_signal,*pdf_bkg),RooArgList(*nSig,*nBkg)); 
    RooDataSet *data = model->generate(*x,nSigInj+nBkgInj);
    
    RooDataHist *roohist = new RooDataHist("roohist","roohist",RooArgList(*x),*data);
    // build fit model
    RooRealVar *nFitSig = new RooRealVar("nFitSig","nFitSig",SIG,0,10*SIG);
    RooRealVar *nFitBkg = new RooRealVar("nFitBkg","nFitBkg",BKG,0,10*BKG);
    RooAddPdf *modelFit = new RooAddPdf("modelFit","modelFit",RooArgList(*pdf_signal,*pdf_bkg),RooArgList(*nFitSig,*nFitBkg)); 
    // fit the pseudo dataset
    RooFitResult *res = modelFit->fitTo(*roohist,RooFit::Save(),RooFit::Extended(kTRUE));
    //res->Print();
    nSigFit = nFitSig->getVal();
    nBkgFit = nFitBkg->getVal();
    eSigFit = nFitSig->getError();
    eBkgFit = nFitBkg->getError();
    nll     = res->minNll();
    tr->Fill();
    if (itoy % 100 == 0) {
      cout<<"Toy #"<<itoy<<": injected = "<<nSigInj<<", fitted = "<<nSigFit<<", error = "<<eSigFit<<endl;
    }
    if (NTOYS == 1) {
      frame = x->frame();
      roohist->plotOn(frame); 
      model->plotOn(frame);
    }
  }
  if (NTOYS == 1) {
    TCanvas *can = new TCanvas("Toy","Toy",900,600);
    frame->Draw();
  }  
  else {
    outf->cd();
    tr->Write();
    outf->Close();
    fTemplates->Close();
  }  
}








