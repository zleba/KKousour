#define SF   TString::Format
#include <vector>

using namespace std;

struct EffPlotter {

    const vector<int> Thresholds = {40, 60 , 80 , 140, 200, 260, 320, 400, 450, 500};
    vector<TH2F*> hDen, hNum;

    vector<vector<TH1D*>> hEff;

    void ReadHistos(TString fileName);
    void CalcEff();
    void PlotEffCurves(int ptFirst, int nPt);


};


void EffPlotter::ReadHistos(TString fileName)
{
    TFile *file = TFile::Open(fileName);

    for(auto tr : Thresholds) {
        TString nDen = SF("HLTJetPt%d_All",  tr);
        TString nNum = SF("HLTJetPt%d_Emulated",  tr);

        TH2F *hD = dynamic_cast<TH2F*>(file->Get(nDen));
        assert(hD);
        TH2F *hN = dynamic_cast<TH2F*>(file->Get(nNum));
        hDen.push_back(hD);
        hNum.push_back(hN);
    }
}

void EffPlotter::CalcEff()
{
    assert(hDen.size() == hNum.size());
    hEff.resize(hDen.size());
    for(int iPt = 0; iPt < hDen.size(); ++iPt) {
        hDen[iPt]->RebinY(6);
        hNum[iPt]->RebinY(6);
        for(int y = 1; y <= hDen[iPt]->GetNbinsY(); ++y) {
            TH1D *hD = hDen[iPt]->ProjectionX(SF("hEff%d%d",iPt,y-1), y, y);
            TH1D *hN = hNum[iPt]->ProjectionX(SF("hNum%d%d",iPt,y-1), y, y);
            hD->Divide(hN, hD, 1., 1., "B");
            hD->SetTitle(SF("Eff%d",Thresholds[iPt+1]));
            hEff[iPt].push_back(hD);
        }
    }
}

void EffPlotter::PlotEffCurves(int ptFirst, int nPt)
{
    double w = 400, h = 500;
    int nY = hEff[ptFirst].size();

    double wTot = w * nY * 1.2;
    double hTot = h * nPt * 1.2;

    TCanvas *can = new TCanvas(SF("cEff%d",ptFirst), "Efficiencies", wTot, hTot);

    cout << "sizes " << nPt << " " << nY << endl;

    can->Divide(nY, nPt, 0, 0);


    for(int pt = ptFirst; pt < ptFirst+nPt; ++pt) {
       for(int y = 0; y < nY; ++y) {
           //hEff[pt][y]->Draw(); 

           double trigMin = 1*Thresholds[pt];
           double trigMax = 3*Thresholds[pt];
           TF1 * f = new TF1(SF("fun%d%d",pt,y),
                   "[0]+0.5*(1-[0])*(1+erf((x-[1])/[2]))",
                   trigMin, trigMax);
           f->SetParameters(0,1.1*trigMin,15);

           int id = (pt-ptFirst) * nY + y + 1;
           can->cd(id);
           gPad->SetBottomMargin(0.1);


           hEff[pt][y]->Fit(f, "QERS", "same", trigMin, trigMax);
           // get turn-on
           double turnon = f->GetX(0.99);

           cout << "TurnOut " << pt <<" "<< y <<" "<< turnon << endl;

           cout << "Id is " << id << endl;

           hEff[pt][y]->Draw();
           hEff[pt][y]->GetXaxis()->SetRangeUser(trigMin, trigMax);
           hEff[pt][y]->GetYaxis()->SetRangeUser(0, 1.1);

           f->Draw("same");

       }
    }

    can->SaveAs("Jindrich.pdf");
}

void DrawEffs()
{
    TH1::SetDefaultSumw2();
    gStyle->SetOptStat(0);
    

    EffPlotter effPl;

    effPl.ReadHistos("/nfs/dust/cms/user/lidrychj/JEC/macros/Data2016output.root");
    effPl.CalcEff();
    effPl.PlotEffCurves(0, 9);



}
