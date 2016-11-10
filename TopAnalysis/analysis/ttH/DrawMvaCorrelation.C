void DrawMvaCorrelation()
{
  gROOT->ForceStyle();
  const int N = 8;
  const int NCAT = 4;
  float XSEC[N] = {0.2934,832,347700,32100,6831,1207,119.9,25.24}; 
  TFile *inf[N]; 
  inf[0] = TFile::Open("Histo_ttHJetTobb_M125_13TeV_amcatnloFXFX_madspin_pythia8.root");
  inf[1] = TFile::Open("Histo_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root");
  inf[2] = TFile::Open("Histo_QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  inf[3] = TFile::Open("Histo_QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  inf[4] = TFile::Open("Histo_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  inf[5] = TFile::Open("Histo_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  inf[6] = TFile::Open("Histo_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");
  inf[7] = TFile::Open("Histo_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root");

  TH2F *h[N][NCAT];
  TH2F *hQCD[NCAT];
  
  TCanvas *can[NCAT];

  char name[100];

  for(int j=0;j<NCAT;j++) {
    sprintf(name,"ttH/h_MvaVsMva_CAT%d",j); 
    for(int i=0;i<N;i++) {
      h[i][j] = (TH2F*)inf[i]->Get(name);
      h[i][j]->Scale(XSEC[i]);
    }
    sprintf(name,"h_MvaVsMva_QCD_CAT%d",j);
    hQCD[j] = (TH2F*)h[2][j]->Clone(name);
    hQCD[j]->Add(h[3][j]);
    hQCD[j]->Add(h[4][j]);
    hQCD[j]->Add(h[5][j]);
    hQCD[j]->Add(h[6][j]);
    hQCD[j]->Add(h[7][j]);

    sprintf(name,"can_MvaVsMva_CAT%d",j);
    can[j] = new TCanvas(name,name,900,600);
    can[j]->Divide(2,2);
    can[j]->cd(1);
    h[0][j]->GetXaxis()->SetTitle("mvaQCD");
    h[0][j]->GetYaxis()->SetTitle("mvaTTbar");
    h[0][j]->Draw("colz");
    can[j]->cd(2);
    h[1][j]->GetXaxis()->SetTitle("mvaQCD");
    h[1][j]->GetYaxis()->SetTitle("mvaTTbar");
    h[1][j]->Draw("colz");
    can[j]->cd(3);
    hQCD[j]->GetXaxis()->SetTitle("mvaQCD");
    hQCD[j]->GetYaxis()->SetTitle("mvaTTbar");
    hQCD[j]->Draw("colz");
  }

  

  
}
