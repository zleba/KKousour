void CreateSigTemplates(double XMIN,double XMAX,double dX)
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
  const int NMASS(5);
  int   H_MASS[NMASS]   = {115,120,125,130,135};
  float XSEC_VBF[NMASS] = {1.215,1.069,0.911,0.746,0.585};
  float XSEC_GF[NMASS]  = {16.14,13.69,11.26,8.93,6.78};
  char name[1000];
  TFile *fVBF[NMASS],*fGF[NMASS];
  TTree *trVBF,*trGF;
  TH1F  *hGF[NMASS][5],*hVBF[NMASS][5],*hTOT[NMASS][5],*hPassGF,*hPassVBF;
  RooDataHist *RooHist[NMASS][5];
  RooAddPdf *model[NMASS][5];
  
  TCanvas *can[NMASS];
  TString PATH("");
  
  RooWorkspace *w = new RooWorkspace("w","workspace");
  int NBINS   = (XMAX-XMIN)/dX;
  
  RooRealVar x("mbbReg","mbbReg",XMIN,XMAX);
  RooRealVar kJES("CMS_scale_j","CMS_scale_j",1,0.9,1.1);
  RooRealVar kJER("CMS_res_j","CMS_res_j",1,0.8,1.2);
  kJES.setConstant(kTRUE);
  kJER.setConstant(kTRUE);

  for(int iMass=0;iMass<NMASS;iMass++) {
    int counter(0);
    for(int iSEL=0;iSEL<2;iSEL++) {
      cout<<"Mass = "<<H_MASS[iMass]<<" GeV"<<endl;
      sprintf(name,"Fit_VBFPowheg%d_sel%s.root",H_MASS[iMass],SELECTION[iSEL].Data());
      fVBF[iMass]  = TFile::Open(PATH+TString(name));
      hPassVBF = (TH1F*)fVBF[iMass]->Get("TriggerPass");
      sprintf(name,"Fit_GFPowheg%d_sel%s.root",H_MASS[iMass],SELECTION[iSEL].Data());
      fGF[iMass]  = TFile::Open(PATH+TString(name)); 
      hPassGF = (TH1F*)fGF[iMass]->Get("TriggerPass");
      sprintf(name,"HMassTemplate_%d_sel%s",H_MASS[iMass],SELECTION[iSEL].Data());
      can[iMass] = new TCanvas(name,name,1200,800);
      can[iMass]->Divide(3,2);
      for(int icat=0;icat<NCAT[iSEL];icat++) { 
        sprintf(name,"HbbCAT%d/events",icat);
        trVBF = (TTree*)fVBF[iMass]->Get(name);
        trGF  = (TTree*)fGF[iMass]->Get(name);
        can[iMass]->cd(icat+1);
        sprintf(name,"mass_VBF%d_sel%s_CAT%d",H_MASS[iMass],SELECTION[iSEL].Data(),icat);
        hVBF[iMass][icat] = new TH1F(name,name,NBINS,XMIN,XMAX);
        hVBF[iMass][icat]->Sumw2();
        trVBF->Draw(MASS_VAR[iSEL]+">>"+TString(name));
        sprintf(name,"mass_GF%d_sel%s_CAT%d",H_MASS[iMass],SELECTION[iSEL].Data(),icat);
        hGF[iMass][icat] = new TH1F(name,name,NBINS,XMIN,XMAX);
        hGF[iMass][icat]->Sumw2();
        trGF->Draw(MASS_VAR[iSEL]+">>"+TString(name));
        delete trVBF;
        delete trGF;
        hGF[iMass][icat]->Scale(LUMI[iSEL]*XSEC_GF[iMass]/hPassGF->GetBinContent(1));
        hVBF[iMass][icat]->Scale(LUMI[iSEL]*XSEC_VBF[iMass]/hPassVBF->GetBinContent(1));
        sprintf(name,"mass_Total%d_sel%s_CAT%d",H_MASS[iMass],SELECTION[iSEL].Data(),icat);
        hTOT[iMass][icat] = (TH1F*)hVBF[iMass][icat]->Clone(name);
        hTOT[iMass][icat]->Add(hGF[iMass][icat]);
      
        sprintf(name,"yield_signal_mass%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar *Yield = new RooRealVar(name,name,hTOT[iMass][icat]->Integral());
      
        sprintf(name,"roohist_mass%d_sel%s_CAT%d",H_MASS[iMass],SELECTION[iSEL].Data(),icat);
        RooHist[iMass][icat] = new RooDataHist(name,name,x,hTOT[iMass][icat]);

        sprintf(name,"mean_m%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar m(name,name,125,100,150);
        sprintf(name,"sigma_m%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar s(name,name,12,3,30);
      
        sprintf(name,"mean_shifted_m%d_CAT%d",H_MASS[iMass],counter);
        RooFormulaVar mShift(name,"@0*@1",RooArgList(m,kJES));
        sprintf(name,"sigma_shifted_m%d_CAT%d",H_MASS[iMass],counter);
        RooFormulaVar sShift(name,"@0*@1",RooArgList(s,kJER)); 
      
        sprintf(name,"alpha_m%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar a(name,name,1,-10,10);
        sprintf(name,"exp_m%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar n(name,name,1,0,100);

        sprintf(name,"b0_m%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar b0(name,name,0.5,0.,1.);
        sprintf(name,"b1_m%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar b1(name,name,0.5,0.,1.);
        sprintf(name,"b2_m%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar b2(name,name,0.5,0.,1.);
        sprintf(name,"b3_m%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar b3(name,name,0.5,0.,1.);

        sprintf(name,"signal_bkg_m%d_CAT%d",H_MASS[iMass],counter);
        RooBernstein bkg(name,name,x,RooArgSet(b0,b1,b2,b3));

        sprintf(name,"fsig_m%d_CAT%d",H_MASS[iMass],counter);
        RooRealVar fsig(name,name,0.7,0.,1.);
      
        sprintf(name,"signal_gauss_m%d_CAT%d",H_MASS[iMass],counter);
        RooCBShape sig(name,name,x,mShift,sShift,a,n);

        // model(x) = fsig*sig(x) + (1-fsig)*bkg(x)
        sprintf(name,"signal_model_m%d_CAT%d",H_MASS[iMass],counter);
        model[iMass][icat] = new RooAddPdf(name,name,RooArgList(sig,bkg),fsig);
      
        RooFitResult *res = model[iMass][icat]->fitTo(*RooHist[iMass][icat],RooFit::Save(),RooFit::SumW2Error(kFALSE),"q");
  
        //res->Print();
      
        RooPlot* frame = x.frame();
        RooHist[iMass][icat]->plotOn(frame);
        //model[iMass][icat]->plotOn(frame,RooFit::VisualizeError(*res,1,kFALSE),RooFit::FillColor(kGray));
        //RooHist[iMass][icat]->plotOn(frame);
        model[iMass][icat]->plotOn(frame);
        //model[iMass][icat]->plotOn(frame,RooFit::LineWidth(2));
        model[iMass][icat]->plotOn(frame,RooFit::Components(bkg),RooFit::LineColor(kBlue),RooFit::LineWidth(2),RooFit::LineStyle(kDashed)); 
        frame->GetXaxis()->SetNdivisions(505); 
        frame->GetXaxis()->SetTitle("M_{bb} (GeV)");
        frame->GetYaxis()->SetTitle("Events");
        frame->Draw();
        hGF[iMass][icat]->SetFillColor(kGreen-8); 
        hVBF[iMass][icat]->SetFillColor(kRed-10); 
        THStack *hs = new THStack("hs","hs");
        hs->Add(hGF[iMass][icat]);
        hs->Add(hVBF[iMass][icat]);
        hs->Draw("same hist");
        frame->Draw("same");
        gPad->RedrawAxis();
        
        TLegend *leg = new TLegend(0.65,0.35,0.9,0.45);  
        leg->AddEntry(hVBF[iMass][icat],"VBF","F");
        leg->AddEntry(hGF[iMass][icat],"GF","F");
        leg->SetFillColor(0);
        leg->SetBorderSize(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.05);
        leg->Draw();
      
        TPaveText *pave = new TPaveText(0.65,0.55,0.9,0.92,"NDC");
        sprintf(name,"M_{H} = %d GeV",H_MASS[iMass]);
        pave->AddText(name);
        sprintf(name,"%s selection",SELECTION[iSEL].Data());
        pave->AddText(name);
        sprintf(name,"CAT%d",icat);
        pave->AddText(name);
        sprintf(name,"m = %1.1f #pm %1.1f",m.getVal(),m.getError());
        pave->AddText(name);
        sprintf(name,"#sigma = %1.1f #pm %1.1f",s.getVal(),s.getError());
        pave->AddText(name);
        /*
        sprintf(name,"a = %1.2f #pm %1.2f",a.getVal(),a.getError());
        pave->AddText(name);
        sprintf(name,"n = %1.2f #pm %1.2f",n.getVal(),n.getError());
        pave->AddText(name);
        sprintf(name,"f = %1.2f #pm %1.2f",fsig.getVal(),fsig.getError());
        pave->AddText(name);
        */
        pave->SetFillColor(0);
        pave->SetBorderSize(0);
        pave->SetTextFont(42);
        pave->SetTextSize(0.05);
        pave->SetTextColor(kBlue);
        pave->Draw();
      
        b0.setConstant(kTRUE);
        b1.setConstant(kTRUE);
        b2.setConstant(kTRUE);
        b3.setConstant(kTRUE);
        m.setConstant(kTRUE);
        s.setConstant(kTRUE); 
        a.setConstant(kTRUE);
        n.setConstant(kTRUE);
        fsig.setConstant(kTRUE);
      
        w->import(*model[iMass][icat]);
        w->import(*RooHist[iMass][icat]);
        w->import(*res); 
        w->import(*Yield);    
       
        counter++;
      }// categories loop
    }// selection loop 
  }// mass loop
  w->Print();
  //x.Print();
  w->writeToFile("signal_shapes_workspace.root");
}
