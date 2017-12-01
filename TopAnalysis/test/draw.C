#define SF   TString::Format
void draw()
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    TFile *f = TFile::Open("/nfs/dust/cms/user/zlebcr/JEC/histos/runG.root");
    TTree *tr = dynamic_cast<TTree*>(f->Get("ak4/events"));
    assert(tr);

    vector<int> lhtTrs ={40,60,80,140,200,260,320,400,450,500};

    TH1D *hBefore = new TH1D("hBefore", "hBefore", 50, 30, 750);
    TH1D *hAfter  = new TH1D("hAfter", "hAfter", 50, 30, 750);

    int i = 5;
    tr->Draw("jetPt[0] >> hBefore", SF("triggerBit[%d] && fabs(jetEta[0]) < 0.5",i));
    tr->Draw("jetPt[0] >> hAfter", SF("triggerBit[%d] && HLTjetPt[0] >= %d && fabs(jetEta[0]) < 0.5",i, lhtTrs[i+1]));

    TCanvas *can = new TCanvas("can", "canvas");

    gPad->SetLogx();
    gPad->SetLogy();
    hBefore->SetTitle(SF("HLT_PFJet%d",lhtTrs[i]));
    hBefore->Draw();
    can->Print("test.pdf(");
    can->Clear();

    gPad->SetLogy(0);
    hRatio = dynamic_cast<TH1D*>(hBefore->Clone("ratio"));
    hRatio->SetTitle(SF("(HLT_PFJet%d (emul) && HLT_PFJet%d)/HLT_PFJet%d",lhtTrs[i+1], lhtTrs[i], lhtTrs[i] ));
    hRatio->Divide (hAfter, hBefore, 1, 1, "B" );
    hRatio->Draw();
    hRatio->GetYaxis()->SetRangeUser(0.9, 1.03);
    can->Print("testYcut.pdf)");
}
