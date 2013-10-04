void CreateDataTemplates(double XMIN,double XMAX,double dX)
{
  gROOT->ForceStyle();
  RooMsgService::instance().setSilentMode(kTRUE);
  for(int i=0;i<2;i++) {
    RooMsgService::instance().setStreamStatus(i,kFALSE);
  }
 
  int NCAT[2] = {5,3};
  char name[1000];
  TString SELECTION[2] = {"NOM","VBF"};
  TString MASS_VAR[2] = {"mbbReg[1]","mbbReg[2]"};
   
  int NBINS   = (XMAX-XMIN)/dX; 
  //RooRealVar x("mbbReg","mbbReg",XMIN,XMAX);
  
  TFile *fBKG  = TFile::Open("bkg_shapes_workspace.root");
  RooWorkspace *wBkg = (RooWorkspace*)fBKG->Get("w");
  RooRealVar x(*(RooRealVar*)wBkg->var("mbbReg"));
  RooWorkspace *w = new RooWorkspace("w","workspace");
  TTree *tr;
  TH1F *h;
  TCanvas *canFit[5]; 
  RooDataHist *roohist[5];

  int counter(0);
  for(int iSEL=0;iSEL<2;iSEL++) {
    TFile *fDATA = TFile::Open("Fit_data_sel"+SELECTION[iSEL]+".root");
    for(int iCAT=0;iCAT<NCAT[iSEL];iCAT++) {
      sprintf(name,"FitData_sel%s_CAT%d",SELECTION[iSEL].Data(),iCAT);
      canFit[iCAT] = new TCanvas(name,name,900,600);
      canFit[iCAT]->cd(1)->SetBottomMargin(0.4);
      sprintf(name,"HbbCAT%d/events",iCAT);
      tr = (TTree*)fDATA->Get(name); 
      sprintf(name,"hMbb_%s_CAT%d",SELECTION[iSEL].Data(),iCAT);
      h = new TH1F(name,name,NBINS,XMIN,XMAX);
      tr->Draw(MASS_VAR[iSEL]+">>"+TString(name));
      sprintf(name,"yield_data_CAT%d",counter);
      RooRealVar *Yield = new RooRealVar(name,name,h->Integral());

      sprintf(name,"data_hist_CAT%d",counter);
      roohist[iCAT] = new RooDataHist(name,name,x,h);
    
      sprintf(name,"b0_CAT%d",counter);
      RooRealVar b0(name,name,0.5,0,1.);
      sprintf(name,"b1_CAT%d",counter);
      RooRealVar b1(name,name,0.5,0,1.);
      sprintf(name,"b2_CAT%d",counter);
      RooRealVar b2(name,name,0.5,0,1.);
      sprintf(name,"b3_CAT%d",counter);
      RooRealVar b3(name,name,0.5,0,1.);
      sprintf(name,"b4_CAT%d",counter);
      RooRealVar b4(name,name,0.5,0,1.);
      sprintf(name,"b5_CAT%d",counter);
      RooRealVar b5(name,name,0.5,0,1.);
    
      sprintf(name,"qcd_model_CAT%d",counter);
      RooBernstein qcd_pdf(name,name,x,RooArgSet(b0,b1,b2,b3,b4,b5));
    
      sprintf(name,"Z_model_CAT%d",counter);
      RooAbsPdf *z_pdf = (RooAbsPdf*)wBkg->pdf(name);
      sprintf(name,"Top_model_CAT%d",counter);
      RooAbsPdf *top_pdf = (RooAbsPdf*)wBkg->pdf(name);

      sprintf(name,"yield_ZJets_CAT%d",counter);
      RooRealVar *nZ = (RooRealVar*)wBkg->var(name);
      sprintf(name,"yield_Top_CAT%d",counter);
      RooRealVar *nT = (RooRealVar*)wBkg->var(name);
      sprintf(name,"yield_QCD_CAT%d",counter);
      RooRealVar nQCD(name,name,1000,0,1e+10);
      nZ->setConstant(kTRUE);
      nT->setConstant(kTRUE);
    
      sprintf(name,"bkg_model_CAT%d",counter);
      RooAddPdf model(name,name,RooArgList(*z_pdf,*top_pdf,qcd_pdf),RooArgList(*nZ,*nT,nQCD));

      RooFitResult *res = model.fitTo(*roohist[iCAT],RooFit::Save());
      res->Print();
    
      RooPlot* frame = x.frame();
      RooPlot* frame1 = x.frame();
      roohist[iCAT]->plotOn(frame);
      model.plotOn(frame,RooFit::LineWidth(2));
      RooHist *hresid = frame->residHist(); 
      //model.plotOn(frame,RooFit::VisualizeError(*res,1,kFALSE),RooFit::FillColor(kGray),RooFit::MoveToBack());
      model.plotOn(frame,RooFit::Components(qcd_pdf),RooFit::LineWidth(2),RooFit::LineColor(kBlack),RooFit::LineStyle(kDashed));
      model.plotOn(frame,RooFit::Components(*z_pdf),RooFit::LineWidth(2),RooFit::LineColor(kBlue));
      model.plotOn(frame,RooFit::Components(*top_pdf),RooFit::LineWidth(2),RooFit::LineColor(kGreen+1));
      frame->Draw(); 
      gPad->Update();
      TPad* pad = new TPad("pad", "pad", 0., 0., 1., 1.);
      pad->SetTopMargin(0.6);
      pad->SetFillColor(0);
      pad->SetFillStyle(0);
      pad->Draw();
      pad->cd(0);
      frame1->addPlotable(hresid,"p");
      frame1->Draw();

      w->import(*roohist[iCAT]);
      w->import(model);
      w->import(*Yield);
      counter++; 
    }// category loop
  }// selection loop
  w->Print();
  w->writeToFile("data_shapes_workspace.root");
}
