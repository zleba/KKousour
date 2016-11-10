void CreateTemplates(TString CAT, TString CUT)
{
  gROOT->ForceStyle();
  
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setStreamStatus(0,kFALSE);
  RooMsgService::instance().setStreamStatus(1,kFALSE);

  const int NMC = 1;

  TString SAMPLE[NMC] = {
    "TT_TuneCUETP8M1_13TeV-powheg-pythia8",
    /*
    "TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8",
    "TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8",
    "TT_TuneCUETP8M1_13TeV-powheg-scaledown-pythia8",
    "TT_TuneCUETP8M1_13TeV-powheg-scaleup-pythia8",
    "TT_TuneCUETP8M1_mtop1665_13TeV-powheg-pythia8",
    "TT_TuneCUETP8M1_mtop1695_13TeV-powheg-pythia8",
    "TT_TuneCUETP8M1_mtop1715_13TeV-powheg-pythia8",
    "TT_TuneCUETP8M1_mtop1735_13TeV-powheg-pythia8",
    "TT_TuneCUETP8M1_mtop1755_13TeV-powheg-pythia8",
    "TT_TuneCUETP8M1_mtop1785_13TeV-powheg-pythia8",
    "TT_TuneCUETP8M1mpiOFF_13TeV-powheg-pythia8",
    "TT_TuneEE5C_13TeV-powheg-herwigpp"
    */
  };
  TFile *infMC[NMC];
  float normMC[NMC];
  for(int k=0;k<NMC;k++) {
    infMC[k] = TFile::Open("Histo_"+SAMPLE[k]+".root");
  }
  TFile *infData = TFile::Open("Histo_JetHT.root");
  
  float XSEC(832.);
  float LUMI(15900.);
  
  const int NSRC = 1;
  TString ALIAS[NSRC] = {
    "Nominal",
    /*
    "NoWeight","BtagUp","BtagDown","TrigUp","TrigDown","JER","JERUp","JERDown","JESUp","JESDown",
    "aMC@NLO","Madgraph","ScaleDown","ScaleUp","mtop1665","mtop1695","mtop1715","mtop1735","mtop1755","mtop1785",
    "mpiOFF","herwigpp"
    */
  };
  RooRealVar *kMassScale = new RooRealVar("kMassScale","kMassScale",1.0,0.5,1.5);
  RooRealVar *kMassResol = new RooRealVar("kMassResol","kMassResol",1.0,0.5,1.5);
  kMassScale->setConstant(kTRUE);
  kMassResol->setConstant(kTRUE);

  RooRealVar *YieldTT[NSRC];
  RooRealVar *AccTT[NSRC];
  TH1F *hMC[NSRC];
  TH1F *hData;
  RooDataHist *roohMC[NSRC];
  RooDataHist *roohData,*roohDataCor,*roohDataCorUp,*roohDataCorDown;
  RooAddPdf *qcd,*qcdCor,*qcdCorUp,*qcdCorDown;
  RooAddPdf *signal[NSRC];
  
  TString VAR,TAG;
  float XMIN,XMAX;

  VAR = "mTop";
  TAG = CUT+"_"+CAT;
  XMIN = 70.;
  XMAX = 300.;

  RooWorkspace *w = new RooWorkspace("w","workspace");

  //---- define observable ------------------------
  RooRealVar *x = new RooRealVar("mTop","mTop",XMIN,XMAX);
  w->import(*x);
  //---- first do the data template ---------------
  
  hData = (TH1F*)infData->Get("boosted/h_mTop_"+CUT+"_0btag");
  hData->Rebin(5);
  roohData = new RooDataHist("roohistData","roohistData",RooArgList(*x),hData);  
  
  //---- QCD -----------------------------------
  RooRealVar bQCD0("qcd_b0","qcd_b0",0.5,0,1);
  RooRealVar bQCD1("qcd_b1","qcd_b1",0.5,0,1);
  RooRealVar bQCD2("qcd_b2","qcd_b2",0.5,0,1);
  RooRealVar bQCD3("qcd_b3","qcd_b3",0.5,0,1);
  RooBernstein qcd1("qcd_brn","qcd_brn",*x,RooArgList(bQCD0,bQCD1,bQCD2,bQCD3));

  RooRealVar mQCD("qcd_mean" ,"qcd_mean",140,130,300);
  RooRealVar sQCD("qcd_sigma","qcd_sigma",50,10,200);
  RooGaussian qcd2("qcd_gaus" ,"qcd_gaus",*x,mQCD,sQCD);

  RooRealVar fqcd("qcd_f","qcd_f",0.5,0,1);

  qcd = new RooAddPdf("qcd_pdf","qcd_pdf",qcd1,qcd2,fqcd);
  
  //---- plots ---------------------------------------------------
  TCanvas *canQCD = new TCanvas("Template_QCD_"+CUT+"_"+CAT,"Template_QCD_"+CUT+"_"+CAT,900,600);

  RooFitResult *res = qcd->fitTo(*roohData,RooFit::Save());
  res->Print();
  RooPlot *frameQCD = x->frame();
  roohData->plotOn(frameQCD);
  qcd->plotOn(frameQCD);
  qcd->plotOn(frameQCD,RooFit::Components("qcd_brn"),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
  qcd->plotOn(frameQCD,RooFit::Components("qcd_gaus"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
  frameQCD->GetXaxis()->SetTitle("m_{t} (GeV)");
  frameQCD->Draw();
  gPad->Update();
  canQCD->Print("plots/"+TString(canQCD->GetName())+".pdf");

  RooArgSet *parsQCD = (RooArgSet*)qcd->getParameters(roohData);
  parsQCD->setAttribAll("Constant",true);

  w->import(*qcd);

  //---- do the bkg templates -------------
  TFile *infW  = TFile::Open("Histo_WJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
  TFile *infDY = TFile::Open("Histo_DYJetsToQQ_HT180_13TeV-madgraphMLM-pythia8.root");
  TFile *infST_tW_top = TFile::Open("Histo_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
  TFile *infST_tW_antitop = TFile::Open("Histo_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root");
  TFile *infST_t_top = TFile::Open("Histo_ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");
  TFile *infST_t_antitop = TFile::Open("Histo_ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root");
  TH1F *hW = (TH1F*)infW->Get("boosted/hWt_"+VAR+"_"+TAG);
  TH1F *hDY = (TH1F*)infW->Get("boosted/hWt_"+VAR+"_"+TAG);
  TH1F *hST_tW_top = (TH1F*)infST_tW_top->Get("boosted/hWt_"+VAR+"_"+TAG);
  TH1F *hST_tW_antitop = (TH1F*)infST_tW_antitop->Get("boosted/hWt_"+VAR+"_"+TAG);
  TH1F *hST_t_top = (TH1F*)infST_t_top->Get("boosted/hWt_"+VAR+"_"+TAG);
  TH1F *hST_t_antitop = (TH1F*)infST_t_antitop->Get("boosted/hWt_"+VAR+"_"+TAG);
  float normW = ((TH1F*)infW->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  float normDY = ((TH1F*)infDY->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  float normST_tW_top = ((TH1F*)infST_tW_top->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  float normST_tW_antitop = ((TH1F*)infST_tW_antitop->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  float normST_t_top = ((TH1F*)infST_t_top->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  float normST_t_antitop = ((TH1F*)infST_t_antitop->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
  hW->Sumw2();
  hDY->Sumw2();
  hST_tW_top->Sumw2();
  hST_tW_antitop->Sumw2();
  hST_t_top->Sumw2();
  hST_t_antitop->Sumw2();
  hW->Scale(3539*LUMI/normW);
  hDY->Scale(1460.*LUMI/normDY);
  hST_tW_top->Scale(35.6*LUMI/normST_tW_top);
  hST_tW_antitop->Scale(35.6*LUMI/normST_tW_antitop);
  hST_t_top->Scale(136.02*LUMI/normST_t_top);
  hST_t_antitop->Scale(80.95*LUMI/normST_t_antitop);
  TH1F *hBkg = (TH1F*)hW->Clone("hBkg");
  
  hBkg->Add(hDY);
  hBkg->Add(hST_tW_top);
  hBkg->Add(hST_tW_antitop);
  hBkg->Add(hST_t_top);
  hBkg->Add(hST_t_antitop);
  
  hBkg->Rebin(10);
  RooDataHist *roohBkg = new RooDataHist("roohistBkg","roohistBkg",RooArgList(*x),hBkg); 
  RooRealVar mBkg1("bkg_mean1","bkg_mean1",172,150,180);
  RooRealVar sBkg1("bkg_sigma1","bkg_sigma1",20,0,50);

  RooGaussian bkg1("bkg_gaus1","bkg_gaus1",*x,mBkg1,sBkg1);

  RooRealVar mBkg2("bkg_mean2","bkg_mean2",90,70,100);
  RooRealVar sBkg2("bkg_sigma2","bkg_sigma2",5,0,20);

  RooGaussian bkg2("bkg_gaus2","bkg_gaus2",*x,mBkg2,sBkg2);
   
  RooRealVar bBkg0("bkg_b0","bkg_b0",0.5,0,1);
  RooRealVar bBkg1("bkg_b1","bkg_b1",0.5,0,1);
  RooRealVar bBkg2("bkg_b2","bkg_b2",0.5,0,1); 
  RooRealVar bBkg3("bkg_b3","bkg_b3",0.5,0,1);
  RooRealVar bBkg4("bkg_b4","bkg_b4",0.5,0,1);
  RooRealVar bBkg5("bkg_b5","bkg_b5",0.5,0,1); 
  RooRealVar bBkg6("bkg_b6","bkg_b6",0.5,0,1);
  RooRealVar bBkg7("bkg_b7","bkg_b7",0.5,0,1);
  RooRealVar bBkg8("bkg_b8","bkg_b8",0.5,0,1);

  RooBernstein bkg3("bkg_bkg","bkg_bkg",*x,RooArgList(bBkg0,bBkg1,bBkg2)); 

  RooRealVar fbkg1("bkg_f1","bkg_f1",0.9,0,1);
  RooRealVar fbkg2("bkg_f2","bkg_f2",0.01,0,1);

  RooAddPdf *bkg = new RooAddPdf("bkg_pdf","bkg_pdf",RooArgList(bkg1,bkg2,bkg3),RooArgList(fbkg1,fbkg2));
  res = bkg->fitTo(*roohBkg,RooFit::Save());  
  res->Print();

  TCanvas *canBkg = new TCanvas("Template_Bkg_"+CUT+"_"+CAT,"Template_Bkg_"+CUT+"_"+CAT,900,600);
  RooPlot *frameBkg = x->frame();
  roohBkg->plotOn(frameBkg);
  bkg->plotOn(frameBkg);
  bkg->plotOn(frameBkg,RooFit::Components("bkg_bkg"),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
  bkg->plotOn(frameBkg,RooFit::Components("bkg_gaus1"),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
  bkg->plotOn(frameBkg,RooFit::Components("bkg_gaus2"),RooFit::LineColor(kOrange+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
  frameBkg->GetXaxis()->SetTitle("m_{t} (GeV)");
  frameBkg->Draw();
  gPad->Update();
  canBkg->Print("plots/"+TString(canBkg->GetName())+".pdf");

  RooArgSet *parsBkg = (RooArgSet*)bkg->getParameters(roohData);
  parsBkg->setAttribAll("Constant",true);

  w->import(*bkg);
  
  //---- then do the signal templates -------------
  hMC[0]  = (TH1F*)infMC[0]->Get("boosted/hWt_"+VAR+"_"+TAG);
  
  /*
  hMC[1]  = (TH1F*)infMC[0]->Get(TYPE+"/h_"+VAR+"_"+TAG);
  hMC[2]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtBtagUp_"+VAR+"_"+TAG);
  hMC[4]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtTrigUp_"+VAR+"_"+TAG);
  if (TYPE == "resolved") {
    hMC[3]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtBtagDo_"+VAR+"_"+TAG);
    hMC[5]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtTrigDo_"+VAR+"_"+TAG);
  }
  else {
    hMC[3]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtBtagDown_"+VAR+"_"+TAG);
    hMC[5]  = (TH1F*)infMC[0]->Get(TYPE+"/hWtTrigDown_"+VAR+"_"+TAG);
  }
  hMC[6]  = (TH1F*)infMC[0]->Get(TYPE+"Smeared/hWt_"+VAR+"_"+TAG);
  hMC[7]  = (TH1F*)infMC[0]->Get(TYPE+"SmearedUp/hWt_"+VAR+"_"+TAG);
  hMC[8]  = (TH1F*)infMC[0]->Get(TYPE+"SmearedDown/hWt_"+VAR+"_"+TAG);
  hMC[9]  = (TH1F*)infMC[0]->Get(TYPE+"ShiftedUp/hWt_"+VAR+"_"+TAG);
  hMC[10] = (TH1F*)infMC[0]->Get(TYPE+"ShiftedDown/hWt_"+VAR+"_"+TAG);
  hMC[11] = (TH1F*)infMC[1]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[12] = (TH1F*)infMC[2]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[13] = (TH1F*)infMC[3]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[14] = (TH1F*)infMC[4]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[15] = (TH1F*)infMC[5]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[16] = (TH1F*)infMC[6]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[17] = (TH1F*)infMC[7]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[18] = (TH1F*)infMC[8]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[19] = (TH1F*)infMC[9]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[20] = (TH1F*)infMC[10]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[21] = (TH1F*)infMC[11]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  hMC[22] = (TH1F*)infMC[12]->Get(TYPE+"/hWt_"+VAR+"_"+TAG);
  */
  TCanvas *canS[NSRC];
  
  for(int k=0;k<NSRC;k++) {
    if (k < 11) {
      normMC[k] = ((TH1F*)infMC[0]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
      //cout<<ALIAS[k]<<" "<<normMC[k]<<" "<<((TH1F*)infMC[0]->Get("eventCounter/GenEventWeight"))->GetEntries()<<endl;
    }
    else {
      normMC[k] = ((TH1F*)infMC[k-10]->Get("eventCounter/GenEventWeight"))->GetSumOfWeights();
      //cout<<ALIAS[k]<<" "<<normMC[k]<<" "<<((TH1F*)infMC[k-10]->Get("eventCounter/GenEventWeight"))->GetEntries()<<endl;
    }
    hMC[k]->Rebin(2);
    double error(0.0);
    float signal_yield = XSEC*LUMI*hMC[k]->IntegralAndError(1,hMC[k]->GetNbinsX(),error)/normMC[k];
    float signal_error = XSEC*LUMI*error/normMC[k];
    
    YieldTT[k] = new RooRealVar("YieldTT_"+ALIAS[k],"YieldTT_"+ALIAS[k],signal_yield);
    YieldTT[k]->setError(signal_error);

    AccTT[k] = new RooRealVar("AccTT_"+ALIAS[k],"AccTT_"+ALIAS[k],signal_yield/(XSEC*LUMI));
    AccTT[k]->setError(signal_error/(XSEC*LUMI));

    cout<<"Yield "<<ALIAS[k]<<": "<<signal_yield<<" +/- "<<signal_error<<endl;
    
    roohMC[k] = new RooDataHist("roohistTT_"+ALIAS[k],"roohistTT_"+ALIAS[k],RooArgList(*x),hMC[k]);    
    RooRealVar m1("ttbar_mean1_"+ALIAS[k],"ttbar_mean1_"+ALIAS[k],172,150,180);
    RooRealVar s1("ttbar_sigma1_"+ALIAS[k],"ttbar_sigma1_"+ALIAS[k],20,0,50);

    RooFormulaVar m1Shift("ttbar_mean1Shifted_"+ALIAS[k],"@0*@1",RooArgList(m1,*(kMassScale)));
    RooFormulaVar s1Shift("ttbar_sigma1Shifted_"+ALIAS[k],"@0*@1",RooArgList(s1,*(kMassResol)));

    RooGaussian sig1("ttbar_gaus1_"+ALIAS[k],"ttbar_gaus2_"+ALIAS[k],*x,m1Shift,s1Shift);

    RooRealVar m2("ttbar_mean2_"+ALIAS[k],"ttbar_mean2_"+ALIAS[k],80,70,90);
    RooRealVar s2("ttbar_sigma2_"+ALIAS[k],"ttbar_sigma2_"+ALIAS[k],5,0,10);

    RooFormulaVar m2Shift("ttbar_mean2Shifted_"+ALIAS[k],"@0*@1",RooArgList(m2,*(kMassScale)));
    RooFormulaVar s2Shift("ttbar_sigma2Shifted_"+ALIAS[k],"@0*@1",RooArgList(s2,*(kMassResol)));
   
    RooGaussian sig2("ttbar_gaus2_"+ALIAS[k],"ttbar_gaus2_"+ALIAS[k],*x,m2Shift,s2Shift);

    RooRealVar bSig0("ttbar_b0_"+ALIAS[k],"ttbar_b0_"+ALIAS[k],0.5,0,1);
    RooRealVar bSig1("ttbar_b1_"+ALIAS[k],"ttbar_b1_"+ALIAS[k],0.5,0,1);
    RooRealVar bSig2("ttbar_b2_"+ALIAS[k],"ttbar_b2_"+ALIAS[k],0.5,0,1); 
    RooRealVar bSig3("ttbar_b3_"+ALIAS[k],"ttbar_b3_"+ALIAS[k],0.5,0,1);
    RooRealVar bSig4("ttbar_b4_"+ALIAS[k],"ttbar_b4_"+ALIAS[k],0.5,0,1);
    RooRealVar bSig5("ttbar_b5_"+ALIAS[k],"ttbar_b5_"+ALIAS[k],0.5,0,1); 
    RooRealVar bSig6("ttbar_b6_"+ALIAS[k],"ttbar_b6_"+ALIAS[k],0.5,0,1);
    RooRealVar bSig7("ttbar_b7_"+ALIAS[k],"ttbar_b7_"+ALIAS[k],0.5,0,1);
    RooRealVar bSig8("ttbar_b8_"+ALIAS[k],"ttbar_b8_"+ALIAS[k],0.5,0,1);

    RooBernstein sig3("ttbar_bkg_"+ALIAS[k],"ttbar_bkg_"+ALIAS[k],*x,RooArgList(bSig0,bSig1,bSig2,bSig3,bSig4,bSig5,bSig6,bSig7)); 

    RooRealVar mL("ttbar_mean3_"+ALIAS[k],"ttbar_mean3_"+ALIAS[k],90,0,200);
    RooRealVar sL("ttbar_sigma3_"+ALIAS[k],"ttbar_sigma3_"+ALIAS[k],10,0,100);
    //RooLandau sig3("ttbar_bkg_"+ALIAS[k],"ttbar_bkg_"+ALIAS[k],*x,mL,sL);

    RooRealVar fsig1("ttbar_f1_"+ALIAS[k],"ttbar_f1_"+ALIAS[k],0.5,0,1);
    RooRealVar fsig2("ttbar_f2_"+ALIAS[k],"ttbar_f2_"+ALIAS[k],0.1,0.01,1);

    RooAddPdf *signal = new RooAddPdf("ttbar_pdf_"+ALIAS[k],"ttbar_pdf_"+ALIAS[k],RooArgList(sig1,sig2,sig3),RooArgList(fsig1,fsig2));

    canS[k] = new TCanvas("Template_TT_"+CAT+"_"+CUT+"_"+ALIAS[k],"Template_TT_"+CAT+"_"+CUT+"_"+ALIAS[k],900,600);

    RooFitResult *res = signal->fitTo(*roohMC[k],RooFit::Save());

    cout<<"mean1 = "<<m1.getVal()<<", sigma1 = "<<s1.getVal()<<endl;
    cout<<"mean2 = "<<m2.getVal()<<", sigma2 = "<<s2.getVal()<<endl;

    res->Print();
    RooPlot *frameS = x->frame();
    roohMC[k]->plotOn(frameS);
    signal->plotOn(frameS);
    signal->plotOn(frameS,RooFit::Components("ttbar_gaus1_"+ALIAS[k]),RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::LineStyle(2));
    signal->plotOn(frameS,RooFit::Components("ttbar_gaus2_"+ALIAS[k]),RooFit::LineColor(kOrange),RooFit::LineWidth(2),RooFit::LineStyle(2));
    signal->plotOn(frameS,RooFit::Components("ttbar_bkg_"+ALIAS[k]),RooFit::LineColor(kGreen+1),RooFit::LineWidth(2),RooFit::LineStyle(2));
    frameS->GetXaxis()->SetTitle("m_{t} (GeV)");
    frameS->Draw();
    gPad->Update();
    canS[k]->Print("plots/"+TString(canS[k]->GetName())+".pdf");

    RooArgSet *parsSig = (RooArgSet*)signal->getParameters(roohMC[k]);
    parsSig->setAttribAll("Constant",true);

    w->import(*signal);
    w->import(*YieldTT[k]);
    w->import(*AccTT[k]);
  }
  
  //w->Print();
  w->writeToFile("templates_"+CUT+"_"+CAT+"_workspace.root");
}                            

