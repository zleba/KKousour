void CreateBkgTemplates(double XMIN,double XMAX)
{
  gROOT->ForceStyle();
  RooMsgService::instance().setSilentMode(kTRUE);
  for(int i=0;i<2;i++) {
    RooMsgService::instance().setStreamStatus(i,kFALSE);
  }
  int NCAT[2] = {5,3};
  float LUMI[2] = {19784,10683};
  TString SELECTION[2] = {"NOM","VBF"};
  TString MASS_VAR[2] = {"mbbReg[1]","mbbReg[2]"};
  TString PATH("");
  TFile *inf[9];
  TTree *tr;
  TH1F *hMbb[9],*hPass;
  TH1F *hZ,*hW,*hTT,*hST,*hTop;
  char name[1000];
  float LUMI;
  float XSEC[9]    = {56.4,11.1,3.79,30.7,11.1,1.76,245.8,650,1.2*1205};
  RooDataHist *roohist_Z[5],*roohist_T[5];
  RooRealVar x("mbbReg","mbbReg",XMIN,XMAX);
  RooWorkspace *w = new RooWorkspace("w","workspace");
  int counter(0);
  for(int iSEL=0;iSEL<2;iSEL++) {
    inf[0] = TFile::Open(PATH+"Fit_T_t-channel_sel"+SELECTION[iSEL]+".root");
    inf[1] = TFile::Open(PATH+"Fit_T_tW-channel_sel"+SELECTION[iSEL]+".root");
    inf[2] = TFile::Open(PATH+"Fit_T_s-channel_sel"+SELECTION[iSEL]+".root");
    inf[3] = TFile::Open(PATH+"Fit_Tbar_t-channel_sel"+SELECTION[iSEL]+".root");
    inf[4] = TFile::Open(PATH+"Fit_Tbar_tW-channel_sel"+SELECTION[iSEL]+".root");
    inf[5] = TFile::Open(PATH+"Fit_Tbar_s-channel_sel"+SELECTION[iSEL]+".root");
    inf[6] = TFile::Open(PATH+"Fit_TTJets_sel"+SELECTION[iSEL]+".root");
    inf[7] = TFile::Open(PATH+"Fit_ZJets_sel"+SELECTION[iSEL]+".root");
    inf[8] = TFile::Open(PATH+"Fit_WJets_sel"+SELECTION[iSEL]+".root");
     
    TCanvas *canZ = new TCanvas("canZ_"+SELECTION[iSEL],"canZ_"+SELECTION[iSEL],900,600); 
    TCanvas *canT = new TCanvas("canT_"+SELECTION[iSEL],"canT_"+SELECTION[iSEL],900,600);  
    if (iSEL == 0) {
      canZ->Divide(3,2);
      canT->Divide(3,2);
    }
    else {
      canZ->Divide(3,1);
     canT->Divide(3,1);
    }
    TCanvas *can = new TCanvas(); 
    for(int iCAT=0;iCAT<NCAT[iSEL];iCAT++) {
      for(int i=0;i<9;i++) {
        hPass = (TH1F*)inf[i]->Get("TriggerPass");
        sprintf(name,"HbbCAT%d/events",iCAT);
        tr = (TTree*)inf[i]->Get(name); 
        sprintf(name,"hMbb%d_sel%s_CAT%d",i,SELECTION[iSEL].Data(),iCAT);
        int NBINS(26);
        if (iCAT > 1 && iCAT<=2) NBINS = 13; 
        if (iCAT > 2) NBINS = 10;
        hMbb[i] = new TH1F(name,name,NBINS,XMIN,XMAX);
        can->cd();
        tr->Draw(MASS_VAR[iSEL]+">>"+TString(name));
        delete tr;
        hMbb[i]->Scale(LUMI[iSEL]*XSEC[i]/hPass->GetBinContent(1)); 
      }
      hZ  = (TH1F*)hMbb[7]->Clone("Z");
      hW  = (TH1F*)hMbb[8]->Clone("W");
      hTT = (TH1F*)hMbb[6]->Clone("TT");
      hST = (TH1F*)hMbb[0]->Clone("ST");
      hST->Add(hMbb[1]);
      hST->Add(hMbb[2]);
      hST->Add(hMbb[3]);
      hST->Add(hMbb[4]);
      hST->Add(hMbb[5]);
      hTop = (TH1F*)hTT->Clone("Top");
      hTop->Add(hST);
      //hZ->Add(hW);
 
      sprintf(name,"yield_ZJets_CAT%d",counter);
      RooRealVar *YieldZ = new RooRealVar(name,name,hZ->Integral());
      sprintf(name,"yield_Top_CAT%d",counter);
      RooRealVar *YieldT = new RooRealVar(name,name,hTop->Integral());

      sprintf(name,"roohist_Z_CAT%d",counter);
      roohist_Z[iCAT] = new RooDataHist(name,name,x,hZ);

      sprintf(name,"Z_peak_mean_CAT%d",counter);
      RooRealVar m(name,name,95,80,110);
      sprintf(name,"Z_peak_sigma_CAT%d",counter);
      RooRealVar s(name,name,12,9,20);

      sprintf(name,"Z_peak_a_CAT%d",counter);
      RooRealVar a(name,name,-1,-10,10);
      sprintf(name,"Z_peak_n_CAT%d",counter);
      RooRealVar n(name,name,1,0,10);
      
      sprintf(name,"Z_model_CAT%d",counter);
      RooCBShape modelZ(name,name,x,m,s,a,n);
    
      RooFitResult *resZ = modelZ.fitTo(*roohist_Z[iCAT],RooFit::Save(),RooFit::SumW2Error(kFALSE),"q"); 
  
      canZ->cd(iCAT+1);
      RooPlot* frame = x.frame();
      roohist_Z[iCAT]->plotOn(frame);
      //modelZ[iCAT]->plotOn(frame,RooFit::VisualizeError(*res,1,kFALSE),RooFit::FillColor(kGray));
      modelZ.plotOn(frame,RooFit::LineWidth(2));
      frame->Draw();
    
      sprintf(name,"roohist_T_CAT%d",counter);
      roohist_T[iCAT] = new RooDataHist(name,name,x,hTop);

      sprintf(name,"Top_peak_mean_CAT%d",counter);
      RooRealVar mT(name,name,120,30,170);
      sprintf(name,"Top_peak_sigma_CAT%d",counter);
      RooRealVar sT(name,name,15,0,100);

      sprintf(name,"Top_peak_a_CAT%d",counter);
      RooRealVar aT(name,name,-1,-20,20);
      sprintf(name,"Top_peak_n_CAT%d",counter);
      RooRealVar nT(name,name,1,0,20);

      sprintf(name,"Top_par0_CAT%d",counter);
      RooRealVar a0(name,name,0.1,0,10);
      sprintf(name,"Top_par1_CAT%d",counter);
      RooRealVar a1(name,name,0.1,0,10);
      sprintf(name,"Top_par2_CAT%d",counter);
      RooRealVar a2(name,name,0.1,0,10);
      sprintf(name,"Top_par3_CAT%d",counter);
      RooRealVar a3(name,name,0.1,0,10);
      sprintf(name,"Top_par4_CAT%d",counter);
      RooRealVar a4(name,name,0.1,0,10);

      sprintf(name,"c0_CAT%d",counter);
      RooRealVar c0(name,name,1,0.1,10);
      sprintf(name,"c1_CAT%d",counter);
      RooRealVar c1(name,name,0.1,0,10);
      sprintf(name,"c2_CAT%d",counter);
      RooRealVar c2(name,name,1,0,3);
      sprintf(name,"c3_CAT%d",counter);
      RooRealVar c3(name,name,30,20,100);

      sprintf(name,"Top_model_CAT%d",counter);
      RooAbsPdf *modelT; 

      if (iCAT < 0) {
        modelT = new RooCBShape(name,name,x,mT,sT,aT,nT);
      }
      else { 
        modelT = new RooGenericPdf(name,"pow(@0-30,@1)*exp(-@2*pow(@0,@3))",RooArgList(x,c0,c1,c2));
      }
    
      RooFitResult *resT = modelT->fitTo(*roohist_T[iCAT],RooFit::Save(),RooFit::SumW2Error(kFALSE),"q");

      canT->cd(iCAT+1);
      RooPlot* frame = x.frame();
      roohist_T[iCAT]->plotOn(frame);
      modelT->plotOn(frame,RooFit::LineWidth(2));
      frame->Draw();

      m.setConstant(kTRUE);
      s.setConstant(kTRUE);
      a.setConstant(kTRUE);
      n.setConstant(kTRUE);
      mT.setConstant(kTRUE);
      sT.setConstant(kTRUE);
      aT.setConstant(kTRUE);
      nT.setConstant(kTRUE);
      a0.setConstant(kTRUE);
      a1.setConstant(kTRUE);
      a2.setConstant(kTRUE);
      a3.setConstant(kTRUE);
      a4.setConstant(kTRUE);
      c0.setConstant(kTRUE);
      c1.setConstant(kTRUE);
      c2.setConstant(kTRUE);
      c3.setConstant(kTRUE);

      w->import(modelZ);
      w->import(*modelT);
      w->import(*YieldZ);
      w->import(*YieldT);
      YieldZ->Print();
      YieldT->Print();

      counter++;
    }// category loop
    delete can;
  }// selection loop
  w->Print();
  w->writeToFile("bkg_shapes_workspace.root");
}
