void CreateTemplatesForW(TString CAT, TString CUT)
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  TFile *infMC   = TFile::Open("Histo_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  TFile *infData = TFile::Open("Histo_JetHT.root");
    
  TH1F *hMC,*hData;
  RooDataHist *roohMC,*roohData;
  RooAddPdf *signal,*qcd;
  
  TString VAR,TAG;
  float XMIN,XMAX;

  VAR = "mW";
  TAG = CUT+"_"+CAT;
  XMIN = 20.;
  XMAX = 160.;

  RooWorkspace *w = new RooWorkspace("w","workspace");

  //---- define observable ------------------------
  RooRealVar *x = new RooRealVar("mW","mW",XMIN,XMAX);
  w->import(*x);
  //---- first do the data template ---------------
  
  hData = (TH1F*)infData->Get("boosted/h_mW_"+CUT+"_0btag");
  hData->Rebin(5);
  roohData = new RooDataHist("roohistData","roohistData",RooArgList(*x),hData);  

  //---- QCD -----------------------------------
  RooRealVar bBkg0("qcd_b0","qcd_b0",0.5,0,1);
  RooRealVar bBkg1("qcd_b1","qcd_b1",0.5,0,1);
  RooRealVar bBkg2("qcd_b2","qcd_b2",0.5,0,1);
  RooRealVar bBkg3("qcd_b3","qcd_b3",0.5,0,1);

  RooBernstein qcd1("qcd_brn","qcd_brn",*x,RooArgList(bBkg0,bBkg1,bBkg2,bBkg3));

  RooRealVar mQCD("qcd_mean2" ,"qcd_mean2",40,0,100);
  RooRealVar sQCD("qcd_sigma2","qcd_sigma2",20,0,50);
  RooGaussian qcd2("qcd_gaus" ,"qcd_gaus",*x,mQCD,sQCD);

  RooRealVar fqcd("qcd_f1","qcd_f1",0.5,0,1);

  qcd = new RooAddPdf("qcd_pdf","qcd_pdf",qcd1,qcd2,fqcd);

  //---- plots ---------------------------------------------------
  TCanvas *canB = new TCanvas("Template_W_QCD_"+CUT+"_"+CAT,"Template_W_QCD_"+CUT+"_"+CAT,900,600);

  RooFitResult *res = qcd->fitTo(*roohData,RooFit::Save());
  res->Print();
  RooPlot *frameB = x->frame();
  roohData->plotOn(frameB);
  qcd->plotOn(frameB);
  qcd->plotOn(frameB,RooFit::Components("qcd_brn"),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
  qcd->plotOn(frameB,RooFit::Components("qcd_gaus"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
  frameB->GetXaxis()->SetTitle("m_{W} (GeV)");
  frameB->Draw();
  gPad->Update();
  canB->Print("plots/"+TString(canB->GetName())+".pdf");

  RooArgSet *parsQCD = (RooArgSet*)qcd->getParameters(roohData);
  parsQCD->setAttribAll("Constant",true);

  w->import(*qcd);

  //---- then do the signal templates -------------
  hMC = (TH1F*)infMC->Get("boosted/hWt_"+VAR+"_"+TAG);
  roohMC = new RooDataHist("roohistTT","roohistTT",RooArgList(*x),hMC);
  normMC = ((TH1F*)infMC->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  hMC->Rebin(2);
      
  RooRealVar m("ttbar_mean","ttbar_mean",90,70,100);
  RooRealVar s("ttbar_sigma","ttbar_sigma",20,0,50);
    
  RooGaussian sig_core("ttbar_core","ttbar_core",*x,m,s);
  
  RooRealVar bSig0("ttbar_b0","ttbar_b0",0.5,0,1);
  RooRealVar bSig1("ttbar_b1","ttbar_b1",0.5,0,1);
  RooRealVar bSig2("ttbar_b2","ttbar_b2",0.5,0,1); 
  RooRealVar bSig3("ttbar_b3","ttbar_b3",0.5,0,1);
  RooRealVar bSig4("ttbar_b4","ttbar_b4",0.5,0,1);
  RooRealVar bSig5("ttbar_b5","ttbar_b5",0.5,0,1); 
  RooRealVar bSig6("ttbar_b6","ttbar_b6",0.5,0,1);
  RooRealVar bSig7("ttbar_b7","ttbar_b7",0.5,0,1);
  RooRealVar bSig8("ttbar_b8","ttbar_b8",0.5,0,1);

  RooBernstein sig_tail("ttbar_tail","ttbar_tail",*x,RooArgList(bSig0,bSig1,bSig2,bSig3,bSig4,bSig5,bSig6,bSig7,bSig8)); 

  RooRealVar fcore("ttbar_fcore","ttbar_fcore",0.5,0,1);

  signal = new RooAddPdf("ttbar_pdf","ttbar_pdf",RooArgList(sig_core,sig_tail),RooArgList(fcore));

  TCanvas *canS = new TCanvas("Template_W_TT_"+CAT+"_"+CUT,"Template_W_TT_"+CAT+"_"+CUT,900,600);

  res = signal->fitTo(*roohMC,RooFit::Save());

  cout<<"mean = "<<m.getVal()<<" +/ "<<m.getError()<<", sigma = "<<s.getVal()<<" +/- "<<s.getError()<<endl;

  //res->Print();
  RooPlot *frameS = x->frame();
  roohMC->plotOn(frameS);
  signal->plotOn(frameS);
  signal->plotOn(frameS,RooFit::Components("ttbar_core"),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
  signal->plotOn(frameS,RooFit::Components("ttbar_tail"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
  frameS->GetXaxis()->SetTitle("m_{W} (GeV)");
  frameS->Draw();
  gPad->Update();
  canS->Print("plots/"+TString(canS->GetName())+".pdf");

  RooArgSet *parsSig = (RooArgSet*)signal->getParameters(roohMC);
  parsSig->setAttribAll("Constant",true);

  w->import(*signal);
  
  //w->Print();
  w->writeToFile("templates_W_"+CUT+"_"+CAT+"_workspace.root");
}                            

