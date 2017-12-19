#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TStyle.h"

#define SF TString::Format

#include "/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/plottingHelper.h"
R__LOAD_LIBRARY(/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/libPlottingHelper.so)

using namespace PlottingHelper;

double dist2(double eta1, double phi1, double eta2, double phi2)
{
    double dPhi = abs(phi2 - phi1);
    dPhi = min(dPhi, 2*M_PI - dPhi);

    return ( (eta2-eta1)*(eta2-eta1) + dPhi*dPhi);
}


//#ifdef __MAKECINT__
//#pragma link C++ class vector<Jet>+;
//#endif

//gInterpreter->GenerateDictionary("Jet","match.h");
#include "match.h"


//#pragma link C++ class vector<Jet> +;


void match()
{
    TFile *f = TFile::Open("/nfs/dust/cms/user/zlebcr/JEC/ntuples2/merged/jetsG.root");
    TTreeReader chsTree("ak4/events", f);
    TTreeReader puppiTree("ak4PUPPI/events", f);

    TTreeReaderArray<float> CHSjetPt = {chsTree, "jetPt"};
    TTreeReaderArray<float> CHSjetEta = {chsTree, "jetEta"};
    TTreeReaderArray<float> CHSjetPhi = {chsTree, "jetPhi"};
    TTreeReaderValue<vector<bool>> CHStriggerBit = {chsTree, "triggerBit"};

    TTreeReaderArray<float> PUPPIjetPt = {puppiTree, "jetPt"};
    TTreeReaderArray<float> PUPPIjetEta = {puppiTree, "jetEta"};
    TTreeReaderArray<float> PUPPIjetPhi = {puppiTree, "jetPhi"};
    TTreeReaderValue<vector<bool>> PUPPItriggerBit = {puppiTree, "triggerBit"};

    TH1::SetDefaultSumw2();

    TH3D *hBalEtaPt = new TH3D("hBalEtaPt", "hBalEtaPt", 5, 0., 4.7, 5, 50., 500., 60, 0.8, 1.2);


    const double jetSel = 50;

    int iEv = -1;
    while(chsTree.Next() && puppiTree.Next()) {
        ++iEv;

        //Trigger 120
        if((*CHStriggerBit).at(3)) continue;

        if(CHSjetPt.GetSize() == 0 || PUPPIjetPt.GetSize() == 0)
            continue;
        if(CHSjetPt[0] < 160) continue;


        for(int i = 0; i < PUPPIjetPt.GetSize(); ++i) {
            if(PUPPIjetPt[i] < jetSel) continue;

            int m = -1;
            for(int j = 0; j < CHSjetPt.GetSize(); ++j) {
                double d2 = dist2(PUPPIjetEta[i], PUPPIjetPhi[i], CHSjetEta[j], CHSjetPhi[j]);
                if(d2 < 0.2*0.2) {
                    m = j;
                    break;
                }
            }

            if(m != -1) {
                //cout << "match " << PUPPIjetPt[i] << " "<< PUPPIjetPt[i] / CHSjetPt[m] << endl;
                double r = PUPPIjetPt[i] / CHSjetPt[m];
                hBalEtaPt->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], r);
            }
            //else
                //hBalEtaPt->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], 1.199);

                //cout << "NO    " << PUPPIjetPt[i] << " "<< endl;

            //cout << i << " "<< PUPPIjetPt[i] << endl;

            //if(i
        }

        //hBalEtaPt


        /*
        cout << "CHS" << endl;
        for(int i = 0; i < CHSjetPt.GetSize(); ++i)
            cout << i << " "<< CHSjetPt[i] << endl;
        cout << "PUPPI" << endl;
        for(int i = 0; i < PUPPIjetPt.GetSize(); ++i)
            cout << i << " "<< PUPPIjetPt[i] << endl;
        */


        if(iEv % 10000 == 0 ) cout << iEv << endl;
        if(iEv > 5e4) break;
        //if(CHSjetPt.GetSize() > 0 && PUPPIjetPt.GetSize() > 0)
            //cout <<iEv<<" "<< CHSjetPt[0] << " "<< PUPPIjetPt[0]  << endl;

    }


    //TH1D *hPt = hBalEtaPt->ProjectionZ("Proj");

    TCanvas *can = new TCanvas("can", "can");

    can->Divide(5, 2, 0.0001, 0.02);

    gStyle->SetOptStat(0);
    //DivideTransparent( 

    for(int i = 1; i <= 5; ++i) {
        can->cd(i);
        TH1D *hEta = hBalEtaPt->ProjectionZ(SF("ProjEta%d",i),  i, i, 0, -1 );
        double l = hBalEtaPt->GetXaxis()->GetBinLowEdge(i);
        double u = hBalEtaPt->GetXaxis()->GetBinUpEdge(i);
        hEta->SetTitle(SF("%1.1f < #eta < %1.1f", l, u));
        hEta->GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
        hEta->Draw();
    }
    for(int i = 1; i <= 5; ++i) {
        can->cd(i+5);
        TH1D *hPt = hBalEtaPt->ProjectionZ(SF("ProjPt%d",i),  0, -1, i, i );
        double l = hBalEtaPt->GetYaxis()->GetBinLowEdge(i);
        double u = hBalEtaPt->GetYaxis()->GetBinUpEdge(i);
        hPt->SetTitle(SF("%3.0f < p_{T} < %3.0f", l, u));
        hPt->GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
        hPt->Draw();
    }
    can->Print("res.pdf(");
    can->Clear();

    //TProfile2D * hProf2 = hBalEtaPt->Project3DProfile();
    //Print grid
    can->Divide(5, 5, 0, 0);
    for(int i = 1; i <= 5; ++i) 
    for(int j = 1; j <= 5; ++j) {
        can->cd((i-1)*5 + j);
        TH1D *hBoth = hBalEtaPt->ProjectionZ(SF("ProjBoth%d%d",i,j),  i, i, j, j);
        //cout <<"Radek "<< i << " "<< j<<" " << hBoth->GetStdDev() << endl;

        double lPt = hBalEtaPt->GetYaxis()->GetBinLowEdge(j);
        double uPt = hBalEtaPt->GetYaxis()->GetBinUpEdge(j);
        double lEt = hBalEtaPt->GetXaxis()->GetBinLowEdge(i);
        double uEt = hBalEtaPt->GetXaxis()->GetBinUpEdge(i);

        cout << "RADEK " << j <<" "<< lPt << " "<< endl;

        //hBoth->SetTitle(SF("#splitline{%3.0f < p_{T} < %3.0f}{%1.1f < #eta < %1.1f}", lPt, uPt, lEt, uEt));
        cout << "RADEKT " << SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt) << endl;
        hBoth->SetTitle(SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt));

        hBoth->GetXaxis()->SetTitle("p_{T}^{PUPPI} / p_{T}^{CHS}");
        hBoth->Draw();
        //hBoth->SetTitle(SF("%3.0f < p_{T} < %3.0f     %1.1f < #eta < %1.1f", lPt, uPt, lEt, uEt));
    }

    can->Print("res.pdf");
    can->Clear();

    can->SetMargin(0.1, 0.1, 0.1, 0.1);

    TProfile2D * hProf = hBalEtaPt->Project3DProfile();
    hProf->BuildOptions(0.8, 1.2, "s");

    hProf->SetTitle("Mean value of p_{T}^{PUPPI} / p_{T}^{CHS}");
    hProf->GetYaxis()->SetTitle("#eta");
    hProf->GetXaxis()->SetTitle("p_{T}");
    hProf->Draw("text");

    can->Print("res.pdf");
    can->Clear();

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

    can->Print("res.pdf)");


    //hPt->Draw();

}


double GetWeight()
{







}