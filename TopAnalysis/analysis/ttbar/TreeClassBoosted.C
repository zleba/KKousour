#define TreeClassBoosted_cxx
#include "TreeClassBoosted.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
//#include "CondFormats/BTauObjects/interface/BTagCalibrationReader.h"
//#include "../BtagSF_subjets.C"
#include <iostream>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDirectory.h>
using std::cin;
using std::cout;
using std::endl;

void TreeClassBoosted::Loop(TDirectory *DIR,bool isMC,bool doBtagSF,bool applyNLOWeight)
{
//---- mass window ------------------------
  float MASS_MIN = 140;
  float MASS_MAX = 200;
  //---- define histograms ------------------
  char name[1000];
  const int NVAR = 12;
  TString var[NVAR] = {"ht","nvtx","nJets","nBJets","met","mJJ","ptJJ","yJJ","dPhiJJ","mva","mTop","mW"}; 
  double dX[NVAR]   = {10,1,1,1,1,1,1,0.01,0.01,0.01,1.0,1.0};
  double XMIN[NVAR] = {0,0,0,0,0,0,0,-3,0,-1.5,0,0};
  double XMAX[NVAR] = {3000,50,15,10,150,13000,4500,3,3.142,1.5,1000,1000};
  //----- get the PU weight --------
  //TFile *puf  = TFile::Open("PileupWeight.root");
  //TH1F *hPuSF = (TH1F*)puf->Get("pileupSF");
  //hPuSF->SetDirectory(0);
  //puf->Close(); 
  //DIR->cd();
  //TH1F *hPuWt = new TH1F("hPUWt","hPUWt",1,-0.5,0.5);
  //----- get the NLO weight --------
  //TFile *nlof  = TFile::Open("NLOWeight.root");
  //TF1 *fNLOWeight = (TF1*)nlof->Get("nloWeight");
  //nlof->Close(); 
  //DIR->cd();
  //----- get the Trigger weight --------
  //TFile *trigf = TFile::Open("TriggerWeightBoosted.root");
  //TF1 *trigSF  = (TF1*)trigf->Get("triggerSF");
  //trigf->Close(); 
  //DIR->cd();
  //TH1F *hTrigWt     = new TH1F("hTrigWt","hTrigWt",100,0.5,1.5);
  //TH1F *hTrigWtUp   = new TH1F("hTrigWtUp","hTrigWtUp",100,0.5,1.5);
  //TH1F *hTrigWtDown = new TH1F("hTrigWtDown","hTrigWtDown",100,0.5,1.5);
  //----- get the Btag efficiencies --------
  //TFile *btagf = TFile::Open("BtagEfficiencyBoosted.root");
  //TH2F *heff_b = (TH2F*)btagf->Get("BtagEff_b");
  //TH2F *heff_c = (TH2F*)btagf->Get("BtagEff_c");
  //TH2F *heff_l = (TH2F*)btagf->Get("BtagEff_uds");
  //TH2F *heff_g = (TH2F*)btagf->Get("BtagEff_g");
  //heff_b->SetDirectory(0);
  //heff_c->SetDirectory(0);
  //heff_l->SetDirectory(0);
  //heff_g->SetDirectory(0);
  //btagf->Close(); 
  //DIR->cd();
  //--- btag scale factors ---------------
  //BTagCalibration calib("CSVv2","../CSVv2_subjets.csv");

  //BTagCalibrationReader readerBC(&calib,BTagEntry::OP_MEDIUM,"lt","central");
  //BTagCalibrationReader readerBCup(&calib,BTagEntry::OP_MEDIUM,"lt","up");
  //BTagCalibrationReader readerBCdo(&calib,BTagEntry::OP_MEDIUM,"lt","down");

  //BTagCalibrationReader readerUDSG(&calib,BTagEntry::OP_MEDIUM,"incl","central");
  //BTagCalibrationReader readerUDSGup(&calib,BTagEntry::OP_MEDIUM,"incl","up");
  //BTagCalibrationReader readerUDSGdo(&calib,BTagEntry::OP_MEDIUM,"incl","down");

  //TH1F *hBtagWt     = new TH1F("hBtagWt","hBtagWt",5000,0,5);
  //TH1F *hBtagWtUp   = new TH1F("hBtagWtUp","hBtagWtUp",5000,0,5);
  //TH1F *hBtagWtDown = new TH1F("hBtagWtDown","hBtagWtDown",5000,0,5);
  //--------------------------------------- 
  const int NCUTS = 5;
  float MVACUT[NCUTS] = {-100.0,0.0,0.1,0.2,0.3};
  TH1F *hVar[NVAR][3][NCUTS];
  TH1F *hVarWtNLO[NVAR][3][NCUTS];
  TH1F *hVarWt[NVAR][3][NCUTS];
  TH1F *hVarWtTrigUp[NVAR][3][NCUTS];
  TH1F *hVarWtTrigDown[NVAR][3][NCUTS];
  TH1F *hVarWtBtagUp[NVAR][3][NCUTS];
  TH1F *hVarWtBtagDown[NVAR][3][NCUTS]; 

  cout<<"Booking histograms.........."<<endl;
  for(int ivar=0;ivar<NVAR;ivar++) {
    int NBINS = (XMAX[ivar]-XMIN[ivar])/dX[ivar];
    for(int j=0;j<3;j++) {
      for(int icut=0;icut<NCUTS;icut++) {
        sprintf(name,"h_%s_Cut%d_%dbtag",var[ivar].Data(),icut,j); 
        hVar[ivar][j][icut] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
        sprintf(name,"hWtNLO_%s_Cut%d_%dbtag",var[ivar].Data(),icut,j); 
        hVarWtNLO[ivar][j][icut] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
        sprintf(name,"hWt_%s_Cut%d_%dbtag",var[ivar].Data(),icut,j); 
        hVarWt[ivar][j][icut] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
        sprintf(name,"hWtTrigUp_%s_Cut%d_%dbtag",var[ivar].Data(),icut,j); 
        hVarWtTrigUp[ivar][j][icut] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
        sprintf(name,"hWtTrigDown_%s_Cut%d_%dbtag",var[ivar].Data(),icut,j); 
        hVarWtTrigDown[ivar][j][icut] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
        sprintf(name,"hWtBtagUp_%s_Cut%d_%dbtag",var[ivar].Data(),icut,j); 
        hVarWtBtagUp[ivar][j][icut] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
        sprintf(name,"hWtBtagDown_%s_Cut%d_%dbtag",var[ivar].Data(),icut,j); 
        hVarWtBtagDown[ivar][j][icut] = new TH1F(name,name,NBINS,XMIN[ivar],XMAX[ivar]);
      }
    }
  }
  cout<<"Booking jet histograms.........."<<endl;
  const int NJETVAR = 10;
  TString varJet[NJETVAR] = {"jetMassSoftDrop","jetPt","jetEta","jetBtag","jetTau32","jetTau31",
                             "jetMassSub0","jetMassSub1","jetPtSub0","jetPtSub1"}; 
  double dXJET[NJETVAR]   = {0.5,1,0.1,0.01,0.01,0.01,1,1,1,1};
  double XMINJET[NJETVAR] = {0,0,-3,-10,0,0,0,0,0,0};
  double XMAXJET[NJETVAR] = {1000,2000,3,1,1,1,200,200,2000,2000};

  TH1F *hJetVar[NJETVAR][3][NCUTS];
  TH1F *hJetVarWt[NJETVAR][3][NCUTS];
  TH1F *hJetVarWtTrigUp[NJETVAR][3][NCUTS];
  TH1F *hJetVarWtTrigDown[NJETVAR][3][NCUTS];
  TH1F *hJetVarWtBtagUp[NJETVAR][3][NCUTS];
  TH1F *hJetVarWtBtagDown[NJETVAR][3][NCUTS];

  TH2F *hmTopVsmW[3][NCUTS];

  for(int ivar=0;ivar<NJETVAR;ivar++) {
    int NBINS = (XMAXJET[ivar]-XMINJET[ivar])/dXJET[ivar];   
    for(int j=0;j<3;j++) {
      for(int icut=0;icut<NCUTS;icut++) { 
        sprintf(name,"h_%s_Cut%d_%dbtag",varJet[ivar].Data(),icut,j);
        hJetVar[ivar][j][icut] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"hWt_%s_Cut%d_%dbtag",varJet[ivar].Data(),icut,j);
        hJetVarWt[ivar][j][icut] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"hWtTrigUp_%s_Cut%d_%dbtag",varJet[ivar].Data(),icut,j);
        hJetVarWtTrigUp[ivar][j][icut] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"hWtTrigDown_%s_Cut%d_%dbtag",varJet[ivar].Data(),icut,j);
        hJetVarWtTrigDown[ivar][j][icut] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"hWtBtagUp_%s_Cut%d_%dbtag",varJet[ivar].Data(),icut,j);
        hJetVarWtBtagUp[ivar][j][icut] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
        sprintf(name,"hWtBtagDown_%s_Cut%d_%dbtag",varJet[ivar].Data(),icut,j);
        hJetVarWtBtagDown[ivar][j][icut] = new TH1F(name,name,NBINS,XMINJET[ivar],XMAXJET[ivar]);
      }
    }
  }

  for(int j=0;j<3;j++) {
    for(int icut=0;icut<NCUTS;icut++) { 
      sprintf(name,"h_mTopVsmW_Cut%d_%dbtag",icut,j);
      hmTopVsmW[j][icut] = new TH2F(name,name,100,0,200,150,0,300);
    }
  }

  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"Reading "<<nentries<<" events"<<endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry % (nentries/10) == 0) cout<<100*jentry/nentries<<"%"<<endl;
    fChain->GetEntry(jentry);
    bool cut_trigger_sig = true;
    bool cut_trigger_ctl = true; 
    if (!isMC) {
      cut_trigger_sig = ((*triggerBit)[2]);
      cut_trigger_ctl = ((*triggerBit)[4]); 
    }
    bool cut_nJets   = (nJets > 1);
    bool cut_leptons = (nLeptons == 0); 
    bool cut_jetPt   = ((*jetPt)[1] > 450);
    bool cut_mass_L  = ((*jetMassSoftDrop)[0]>70 && (*jetMassSoftDrop)[0]<300);
    bool cut_mass_T  = ((*jetMassSoftDrop)[0]>MASS_MIN && (*jetMassSoftDrop)[0]<MASS_MAX);
    bool cut_mass    = cut_mass_L;

    if (!cut_trigger_sig && !cut_trigger_ctl)  continue;
    if (!cut_nJets)    continue;
    if (!cut_leptons)  continue;
    if (!cut_jetPt)    continue; 
    if (!cut_mass_L)   continue;

    bool cut_0btag = (((*jetNBSub)[0] == 0) && ((*jetNBSub)[1] == 0));
    bool cut_1btag = ((((*jetNBSub)[0] > 0) && ((*jetNBSub)[1] == 0)) || (((*jetNBSub)[0] == 0) && ((*jetNBSub)[1] > 0)));
    bool cut_2btag = (((*jetNBSub)[0] > 0) && ((*jetNBSub)[1] > 0));

    //---- nlo weight -----------------
    float wtNLO = 1.0;

    //---- gen weight -----------------
    float wtGen = 1.0;

    //---- pu weight ----------------
    float wtPU = 1.0; 

    //---- trigger weight ------------------
    float wtTrig = 1.0; 
    float wtTrigUp = 1.0;
    float wtTrigDown = 1.0;

    //---- btag weight ----------------
    float wtBtag = 1.0;
    float wtBtagErr = 0.0;

    if (isMC) {
      wtGen = genEvtWeight;
      /*
      wtTrig = trigSF->Eval((*jetPt)[0]);
      wtTrigUp = wtTrig+fabs(1-wtTrig)/2;
      wtTrigDown = wtTrig-fabs(1-wtTrig)/2;
      int pubin = hPuSF->FindBin(npu);
      if (pubin < 40) {
        wtPU = hPuSF->GetBinContent(pubin);
      }
      else {
        wtPU = hPuSF->GetBinContent(40);
      }
      if (doBtagSF) {
        std::vector<int>   subjet_flavor;
        std::vector<float> subjet_pt,subjet_eta,subjet_btag;
        for(int i=0;i<2;i++) {
          subjet_pt.push_back((*jetPtSub0)[i]);
          subjet_pt.push_back((*jetPtSub1)[i]);
          subjet_eta.push_back((*jetEtaSub0)[i]);
          subjet_eta.push_back((*jetEtaSub1)[i]);
          subjet_btag.push_back((*jetBtagSub0)[i]);
          subjet_btag.push_back((*jetBtagSub1)[i]);
          subjet_flavor.push_back((*jetFlavorSub0)[i]);
          subjet_flavor.push_back((*jetFlavorSub1)[i]);
        } 
        
        get_weight_btag_subjets(readerBC,readerBCup,readerBCdo,readerUDSG,readerUDSGup,readerUDSGdo,
                      subjet_pt,subjet_eta,subjet_flavor,subjet_btag,heff_b,heff_c,heff_l,heff_g,wtBtag,wtBtagErr);
        if (wtBtag != wtBtag) wtBtag = 1.0;
        if (wtBtagErr != wtBtagErr) wtBtagErr = 0.0;
      }
      */
    }
    /*
    hPuWt->Fill(0.0,wtPU);
    hTrigWt->Fill(wtTrig);
    hTrigWtUp->Fill(wtTrigUp);
    hTrigWtDown->Fill(wtTrigDown);
    hBtagWt->Fill(wtBtag);
    hBtagWtUp->Fill(wtBtag+wtBtagErr);
    hBtagWtDown->Fill(wtBtag-wtBtagErr);

    if (applyNLOWeight) {
      wtNLO = fNLOWeight->Eval(TMath::Max((*partonPt)[0],(*partonPt)[1]));
    }
    */

    float mTop = (*jetMassSoftDrop)[0];
    float mW   = (*jetMassSub0)[0];
    float x[NVAR] = {ht,float(nvtx),float(nJets),float(nBJets),met,mJJ,ptJJ,yJJ,dPhiJJ,mva,mTop,mW};
     
    //---- 2D histogram ----------------
    for(int icut=0;icut<NCUTS;icut++) {
      bool cut_mva = (mva > MVACUT[icut]);
      if (cut_mva) {
        if (cut_trigger_sig && cut_1btag) {
          hmTopVsmW[1][icut]->Fill(mW,mTop);
        } 
        if (cut_trigger_sig && cut_2btag) {
          hmTopVsmW[2][icut]->Fill(mW,mTop);
        }
        if (cut_trigger_ctl && cut_0btag) {
          hmTopVsmW[0][icut]->Fill(mW,mTop);
        }
      }
    }

    for(int ivar=0;ivar<NVAR;ivar++) {
      if (var[ivar].CompareTo("mTop") == 0) {
        cut_mass = cut_mass_L;
      }
      else {
        cut_mass = cut_mass_T;
      } 

      if (!cut_mass) continue;

      if (cut_trigger_sig && cut_1btag) {
        for(int icut=0;icut<NCUTS;icut++) {
          bool cut_mva = (mva > MVACUT[icut]);
          if (cut_mva) {
            hVar[ivar][1][icut]->Fill(x[ivar]);
            hVarWtNLO[ivar][1][icut]->Fill(x[ivar],wtNLO);
            hVarWt[ivar][1][icut]->Fill(x[ivar],wtGen*wtPU*wtTrig*wtBtag);
            hVarWtTrigUp[ivar][1][icut]->Fill(x[ivar],wtGen*wtPU*wtTrigUp*wtBtag);
            hVarWtTrigDown[ivar][1][icut]->Fill(x[ivar],wtGen*wtPU*wtTrigDown*wtBtag);
            hVarWtBtagUp[ivar][1][icut]->Fill(x[ivar],wtGen*wtPU*wtTrig*(wtBtag+wtBtagErr));
            hVarWtBtagDown[ivar][1][icut]->Fill(x[ivar],wtGen*wtPU*wtTrig*(wtBtag-wtBtagErr));
          }
        }
      }
      if (cut_trigger_sig && cut_2btag) {
        for(int icut=0;icut<NCUTS;icut++) {
          bool cut_mva = (mva > MVACUT[icut]);
          if (cut_mva) { 
            hVar[ivar][2][icut]->Fill(x[ivar]);
            hVarWtNLO[ivar][2][icut]->Fill(x[ivar],wtNLO);
            hVarWt[ivar][2][icut]->Fill(x[ivar],wtGen*wtPU*wtTrig*wtBtag);
            hVarWtTrigUp[ivar][2][icut]->Fill(x[ivar],wtGen*wtPU*wtTrigUp*wtBtag);
            hVarWtTrigDown[ivar][2][icut]->Fill(x[ivar],wtGen*wtPU*wtTrigDown*wtBtag);
            hVarWtBtagUp[ivar][2][icut]->Fill(x[ivar],wtGen*wtPU*wtTrig*(wtBtag+wtBtagErr));
            hVarWtBtagDown[ivar][2][icut]->Fill(x[ivar],wtGen*wtPU*wtTrig*(wtBtag-wtBtagErr));
          }
        } 
      }
      if (cut_trigger_ctl && cut_0btag) {
        for(int icut=0;icut<NCUTS;icut++) {
          bool cut_mva = (mva > MVACUT[icut]);
          if (cut_mva) {
            hVar[ivar][0][icut]->Fill(x[ivar]);
            hVarWtNLO[ivar][0][icut]->Fill(x[ivar],wtNLO);
            hVarWt[ivar][0][icut]->Fill(x[ivar],wtGen*wtPU*wtTrig*wtBtag);
            hVarWtTrigUp[ivar][0][icut]->Fill(x[ivar],wtGen*wtPU*wtTrigUp*wtBtag);
            hVarWtTrigDown[ivar][0][icut]->Fill(x[ivar],wtGen*wtPU*wtTrigDown*wtBtag);
            hVarWtBtagUp[ivar][0][icut]->Fill(x[ivar],wtGen*wtPU*wtTrig*(wtBtag+wtBtagErr));
            hVarWtBtagDown[ivar][0][icut]->Fill(x[ivar],wtGen*wtPU*wtTrig*(wtBtag-wtBtagErr));
          }
        }
      } 
    }
    
    for(int ivar=0;ivar<NJETVAR;ivar++) {
      if (varJet[ivar].CompareTo("jetMassSoftDrop") == 0) {
        cut_mass = cut_mass_L;
      }
      else {
        cut_mass = cut_mass_T;
      } 
      if (!cut_mass) continue;
      
      for(int j=0;j<2;j++) {
        float xJet[NJETVAR] = {(*jetMassSoftDrop)[j],(*jetPt)[j],(*jetEta)[j],(*jetBtag)[j],(*jetTau3)[j]/(*jetTau2)[j],(*jetTau3)[j]/(*jetTau1)[j],(*jetMassSub0)[j],(*jetMassSub1)[j],(*jetPtSub0)[j],(*jetPtSub1)[j]}; 

        if (cut_trigger_sig && cut_1btag) {  
          for(int icut=0;icut<NCUTS;icut++) {
            bool cut_mva = (mva > MVACUT[icut]);
            if (cut_mva) {
              hJetVar[ivar][1][icut]->Fill(xJet[ivar]);
              hJetVarWt[ivar][1][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrig*wtBtag);
              hJetVarWtTrigUp[ivar][1][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrigUp*wtBtag);
              hJetVarWtTrigDown[ivar][1][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrigDown*wtBtag);
              hJetVarWtBtagUp[ivar][1][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrig*(wtBtag+wtBtagErr));
              hJetVarWtBtagDown[ivar][1][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrig*(wtBtag-wtBtagErr));
            }
          }
        }
        if (cut_trigger_sig && cut_2btag) {
          for(int icut=0;icut<NCUTS;icut++) {
            bool cut_mva = (mva > MVACUT[icut]);
            if (cut_mva) {
              hJetVar[ivar][2][icut]->Fill(xJet[ivar]);
              hJetVarWt[ivar][2][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrig*wtBtag);
              hJetVarWtTrigUp[ivar][2][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrigUp*wtBtag);
              hJetVarWtTrigDown[ivar][2][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrigDown*wtBtag);
              hJetVarWtBtagUp[ivar][2][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrig*(wtBtag+wtBtagErr));
              hJetVarWtBtagDown[ivar][2][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrig*(wtBtag-wtBtagErr)); 
            }
          }
        }
        if (cut_trigger_ctl && cut_0btag) {
          for(int icut=0;icut<NCUTS;icut++) {
            bool cut_mva = (mva > MVACUT[icut]);
            if (cut_mva) {
              hJetVar[ivar][0][icut]->Fill(xJet[ivar]);
              hJetVarWt[ivar][0][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrig*wtBtag);
              hJetVarWtTrigUp[ivar][0][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrigUp*wtBtag);
              hJetVarWtTrigDown[ivar][0][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrigDown*wtBtag);
              hJetVarWtBtagUp[ivar][0][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrig*(wtBtag+wtBtagErr));
              hJetVarWtBtagDown[ivar][0][icut]->Fill(xJet[ivar],wtGen*wtPU*wtTrig*(wtBtag-wtBtagErr));
            }
          }
        }
      }
    }    
  }
  DIR->Write();
}
