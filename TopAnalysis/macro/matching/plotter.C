#define SF TString::Format

#include "/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/plottingHelper.h"
R__LOAD_LIBRARY(/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/libPlottingHelper.so)

using namespace PlottingHelper;

struct PLOTTER {

    const int nPer = 8;
    int perID;

    vector<TH3D *> hBalEtaPtAll;
    vector<TH1D *> hJetPtAll;

    TString outName = "newPlots.pdf";

    public:
    void Init(TString inFile, TString outFile) {
        outName = outFile;
        //TString fileName = "histoDir/jetsAll.root";
        TFile *fIn = TFile::Open(inFile, "READ");

        hBalEtaPtAll.resize(nPer);
        hJetPtAll.resize(nPer);

        for(int i = 0; i < nPer; ++i) {
            char per = 'A' + i;
            hBalEtaPtAll[i] = dynamic_cast<TH3D*>(fIn->Get(SF("hBalEtaPt_%c",per)));
            assert(hBalEtaPtAll[i]);
            hJetPtAll[i] = dynamic_cast<TH1D*>(fIn->Get(SF("hJetPt_%c",per)));
            assert(hJetPtAll[i]);
        }
    }

    void AsymmetryEtaDep();
    void AsymmetryPtDep();
    void AsymmetryEtaPtDep();
    void MeanAsym();

};


pair<TPad*, TPad*> TitleSpace(int perID)
{
    TVirtualPad *can = gPad;
    const double sp = 0.97;
    TPad *up = new TPad(SF("%d",rand()),"",   0, sp, 1, 1.0);

    up->Draw();

    up->cd();
    TLatex *lat = new TLatex;
    TString s;
    if(perID != 0) s =  SF("Run 2016%c", perID + 'A');
    else s =  SF("Run 2016 all");
    lat->SetTextSize(0.6);
    lat->DrawLatexNDC(0.8, 0.3, s);


    can->cd();


    TPad *dn = new TPad(SF("%d",rand()),"", 0, 0.0, 1, sp);
    dn->Draw();
    dn->cd();
    //return make_pair(up, dn);
    return make_pair(up, dn);
}




void PLOTTER::AsymmetryEtaDep()
{
    hBalEtaPt = hBalEtaPtAll[perID];

    TCanvas *can = new TCanvas("can", "can");

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    //can->
    dnPad->Divide(5, 4, 0.0001, 0.02);
    //DivideTransparent(

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 20; ++i) {
        dnPad->cd(i);
        int iS = 1*i + 1;
        int iE = 1*(i+1) -1 + 1;
        TH1D *hEta = hBalEtaPt->ProjectionZ(SF("ProjEta%d",i),  iS, iE, 0, -1 );
        double l = hBalEtaPt->GetXaxis()->GetBinLowEdge(iS);
        double u = hBalEtaPt->GetXaxis()->GetBinUpEdge(iE);
        hEta->SetTitle(SF("%1.1f < #eta < %1.1f", l, u));
        hEta->GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
        hEta->Draw();
    }
    
    can->Print(outName);
    //can->Print(outName);
    can->Clear();
    delete can;
}


void PLOTTER::AsymmetryPtDep()
{
    hBalEtaPt = hBalEtaPtAll[perID];

    TCanvas *can = new TCanvas("can", "can");
    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    dnPad->Divide(5, 3, 0.0001, 0.02);

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 15; ++i) {
        dnPad->cd(i);
        int iS = 3*i + 6;
        int iE = 3*(i+1) -1 + 6;
        TH1D *hPt = hBalEtaPt->ProjectionZ(SF("ProjPt%d",i),  0, -1, iS, iE );
        double l = hBalEtaPt->GetYaxis()->GetBinLowEdge(iS);
        double u = hBalEtaPt->GetYaxis()->GetBinUpEdge(iE);
        hPt->SetTitle(SF("%3.0f < p_{T} < %3.0f", l, u));
        hPt->GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
        hPt->Draw();
    }
    //can->Print(outName +"(");
    can->Print(outName);
    can->Clear();
    delete can;
}





void PLOTTER::AsymmetryEtaPtDep()
{
    hBalEtaPt = hBalEtaPtAll[perID];

    const vector<int> etaBins = {1, 4, 7, 10, 13, 16};
    const vector<int> ptBins = {6, 12, 18, 24, 30, 36};

    //TProfile2D * hProf2 = hBalEtaPt->Project3DProfile();
    //Print grid

    TCanvas *can = new TCanvas("can", "can");
    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);

    dnPad->Divide(5, 5, 0, 0);
    for(int i = 1; i <= 5; ++i) 
    for(int j = 1; j <= 5; ++j) {
        dnPad->cd((i-1)*5 + j);
        int iS = etaBins[i];
        int iE = etaBins[i+1]-1;
        int jS = ptBins[j];
        int jE = ptBins[j+1]-1;

        TH1D *hBoth = hBalEtaPt->ProjectionZ(SF("ProjBoth%d%d",i,j),  iS, iE, jS, jE);
        //cout <<"Radek "<< i << " "<< j<<" " << hBoth->GetStdDev() << endl;

        double lPt = hBalEtaPt->GetYaxis()->GetBinLowEdge(jS);
        double uPt = hBalEtaPt->GetYaxis()->GetBinUpEdge(jE);
        double lEt = hBalEtaPt->GetXaxis()->GetBinLowEdge(iS);
        double uEt = hBalEtaPt->GetXaxis()->GetBinUpEdge(iE);

        //cout << "RADEK " << j <<" "<< lPt << " "<< endl;

        //hBoth->SetTitle(SF("#splitline{%3.0f < p_{T} < %3.0f}{%1.1f < #eta < %1.1f}", lPt, uPt, lEt, uEt));
        //cout << "RADEKT " << SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt) << endl;
        hBoth->SetTitle(SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt));

        hBoth->GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
        hBoth->Draw();
        //hBoth->SetTitle(SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt));
    }

    can->Print(outName);
    can->Clear();
    delete can;
}


void PLOTTER::MeanAsym()
{
    TH3D *hBalEtaPt =dynamic_cast<TH3D*>( hBalEtaPtAll[perID]->Clone(SF("%d",rand())));
    assert(hBalEtaPt);

    hBalEtaPt->RebinY(4);
    gStyle->SetPaintTextFormat("2.2f");
    TCanvas *can = new TCanvas("can", "can");

    TPad *upPad, *dnPad;
    tie(upPad, dnPad) = TitleSpace(perID);


    dnPad->SetMargin(0.1, 0.1, 0.1, 0.1);
    dnPad->SetLogx();

    TProfile2D * hProf = hBalEtaPt->Project3DProfile();
    hProf->GetXaxis()->SetRangeUser(20, 2000);
    hProf->BuildOptions(0.8, 1.2, "s");

    hProf->SetTitle("Mean value of p_{T}^{PUPPI} / p_{T}^{CHS}");
    hProf->GetYaxis()->SetTitle("#eta");
    hProf->GetXaxis()->SetTitle("p_{T}");
    hProf->GetXaxis()->SetMoreLogLabels();
    hProf->Draw("text");

    can->Print(outName);
    can->Clear();
    delete hProf;
    delete can;
}







void plotter()
{
    PLOTTER plot;
    TString myOut = "newPlotsAug.pdf";
    plot.Init("histoDir/Aug/jetsAll.root", myOut+"(");
    plot.perID = 0;

    for(int i = 0; i < 8; ++i) {
        plot.perID = i;
        plot.AsymmetryEtaDep();
        plot.outName = myOut;
        plot.AsymmetryPtDep();
        plot.AsymmetryEtaPtDep();

        if(i == 7)
            plot.outName = myOut + ")";

        plot.MeanAsym();
    }

    /*

    TH2D *hSigma = hProf->ProjectionXY("sigma");  

    for(int i = 1; i <= hProf->GetNbinsX(); ++i)
    for(int j = 1; j <= hProf->GetNbinsY(); ++j) {
        TH1D *hTemp = hBalEtaPt->ProjectionZ(SF("test%d%d",i,j), i, i, j, j);
        double s = 0;
        for(int k = 1; k <= hTemp->GetNbinsX(); ++k)
            s += pow(hTemp->GetBinCenter(k) - hTemp->GetMean(), 2) * hTemp->GetBinContent(k);
        s /= hTemp->Integral();

        cout << "Holka " << hTemp->GetStdDev()  << endl;
        hSigma->SetBinContent(i, j, hTemp->GetStdDev() );
        //hSigma->SetBinContent(i, j, sqrt(s) );
        delete hTemp;
    }

    hSigma->SetTitle("STD of p_{T}^{PUPPI} / p_{T}^{CHS}");
    hSigma->GetYaxis()->SetTitle("#eta");
    hSigma->GetXaxis()->SetTitle("p_{T}");

    hSigma->Draw("text");

    hJetPt->Draw();
    can->Print("res.pdf)");

    */
}
