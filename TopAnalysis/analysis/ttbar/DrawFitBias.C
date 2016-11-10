void DrawFitBias()
{
  gROOT->ForceStyle();

  const int N = 5;
  TString TAG[N]   = {"Cut0_2btag","Cut1_2btag","Cut2_2btag","Cut3_2btag","Cut4_2btag"};
  TString ALIAS[N] = {"No cut","F > 0.0","F > 0.1","F > 0.2","F > 0.3"};
  int COLOR[N]   = {kBlack,kBlue,kGreen+1,kMagenta,kRed+1};
  TFile *inf[N];
  TTree *tr[N];
  TH1F *hPull[N],*hError[N];
  TH1F *hPullMean   = new TH1F("hPullMean","hPullMean",N,0,N);
  TH1F *hPullRMS    = new TH1F("hPullRMS","hPullRMS",N,0,N);
  TH1F *hErrorAll   = new TH1F("hErrorAll","hErrorAll",N,0,N);
  TCanvas *canPull  = new TCanvas("ToysPull","ToysPull",1000,500);
  TCanvas *canError = new TCanvas("ToysError","ToysError",900,600);
  for(int i=0;i<N;i++) {
    inf[i] = TFile::Open("FitBiasToys_"+TAG[i]+".root");
    tr[i]  = (TTree*)inf[i]->Get("toys");
    hPull[i]  = new TH1F("hPull_"+TAG[i],"hPull_"+TAG[i],100,-5,5);
    hError[i] = new TH1F("hError_"+TAG[i],"hError_"+TAG[i],100,0.02,0.08);
    hPull[i]->SetLineColor(COLOR[i]);
    hError[i]->SetLineColor(COLOR[i]);
    tr[i]->Draw("(nSigFit-nSigInj)/eSigFit>>"+TString(hPull[i]->GetName()));
    tr[i]->Draw("eSigFit/nSigFit>>"+TString(hError[i]->GetName()));
    float m_pull  = hPull[i]->GetMean();
    float s_pull  = hPull[i]->GetRMS();
    float em_pull = hPull[i]->GetMeanError();
    float es_pull = hPull[i]->GetRMSError();
    float m_error = hError[i]->GetMean();
    float e_error = hError[i]->GetMeanError();
    hPullMean->SetBinContent(i+1,m_pull);
    hPullMean->SetBinError(i+1,em_pull);
    hPullRMS->SetBinContent(i+1,s_pull);
    hPullRMS->SetBinError(i+1,es_pull);
    hErrorAll->SetBinContent(i+1,m_error);
    hErrorAll->SetBinError(i+1,e_error);
    hPullMean->GetXaxis()->SetBinLabel(i+1,ALIAS[i]);
    hPullRMS->GetXaxis()->SetBinLabel(i+1,ALIAS[i]);
    hErrorAll->GetXaxis()->SetBinLabel(i+1,ALIAS[i]);
  }
 
  canPull->Divide(2,1);
  canPull->cd(1);
  gPad->SetGridy();
  hPullMean->GetYaxis()->SetTitle("Pull mean"); 
  hPullMean->GetYaxis()->SetRangeUser(-0.08,0.08);
  hPullMean->Draw("E");
  canPull->cd(2);
  gPad->SetGridy();
  hPullRMS->GetYaxis()->SetTitle("Pull RMS");
  hPullRMS->GetYaxis()->SetRangeUser(0.85,1.1);
  hPullRMS->Draw("E");
  
  canError->cd();
  gPad->SetGridy();
  hErrorAll->GetYaxis()->SetTitle("Relative uncertainty");
  hErrorAll->Draw("E"); 
  
}
