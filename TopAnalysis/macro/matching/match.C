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
#include <map>
#include <set>
#include <vector>

#define SF TString::Format

#include "/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/plottingHelper.h"
R__LOAD_LIBRARY(/afs/desy.de/user/c/connorpa/Libraries/PlottingHelper/libPlottingHelper.so)

using namespace PlottingHelper;




class Luminosity {
    vector<map<int, double>> lumiMap;


    public:

    static map<int,double> GetLumis(TString fName)
    {
        ifstream file(fName);
        assert(file.good());
        
        map<int,double> lumiMap;
        while(1) {
            TString trigName;
            double lum;
            file >> trigName >> lum;
            if(!file.good()) break;
            trigName.ReplaceAll("HLT_PFJet", "");
            trigName = trigName(0, trigName.Length()-3);

            int trigTr = trigName.Atoi();
            if(trigTr != 500) {
                lumiMap[trigTr] += lum;
            }
        }
        file.close();
        return lumiMap;

    }

    void LoadLumis() {
        vector<vector<TString>> files;
        files.push_back({"processedLumis2016B.json.txt"});
        files.push_back({"processedLumis2016C.json.txt"});
        files.push_back({"processedLumis2016D.json.txt"});
        files.push_back({"processedLumis2016E.json.txt"});
        files.push_back({"processedLumis2016Fearly.json.txt", "processedLumis2016Flate.json.txt" }); //Merging F into one period
        files.push_back({"processedLumis2016G.json.txt"});
        files.push_back({"processedLumis2016H.json.txt"});

        lumiMap.resize(files.size() + 1); //Index 0 is dedicated for x-section summed over all periods

        TString path = "/nfs/dust/cms/user/connorpa/SMPJ/effective_luminosity_Run2016BtoH/";

        //Read luminosities to map
        for(int i = 0; i < files.size(); ++i) {
            for(int j = 0; j < files[i].size(); ++j) {
                auto lum = GetLumis(path + files[i][j]);
                for(auto v : lum) {
                    lumiMap[i+1][v.first] += v.second;
                    //cout << "Printout " << v.first <<" "<< v.second << endl;
                }
            }
        }

        //Evaluate the total luminosity
        for(int i = 1; i < lumiMap.size(); ++i)
            for(auto v : lumiMap[i])
                lumiMap[0][v.first] += v.second;

        /*
        //Sum B,C,D and make them identical
        map<int,double> lumiMapSingle;
        for(auto v : lumiMap[0]) {
            double sum = lumiMap[1][v.first] + lumiMap[2][v.first] + lumiMap[3][v.first];
            lumiMap[1][v.first] = lumiMap[2][v.first] = lumiMap[3][v.first] = sum;
        }
        */



        //Assert identical trigger configuration:
        for(int i = 1; i < lumiMap.size(); ++i) {
            for(auto m : lumiMap[0])
                assert(lumiMap[i].count(m.first));
            //assert(lumiMap[0].size() == lumiMap[i].size());
        }


        for(auto v : lumiMap[0]) {
            cout << v.first <<" "<< v.second << endl;
        }


    }

    //runLum, totLum, trigID
    tuple<double,double,int> GetWeightID(int periodId, double pT0)
    {
        //Trigger tresholds from Patrick
        static const map<int,int> trigTrsh = {
            {40 , 74},
            {60 , 84},
            {80 , 114},
            {140 , 196},
            {200 , 245},
            {260 , 330},
            {320 , 395},
            {400 , 468},
            {450 , 507},
        };

        assert(periodId < lumiMap.size());
        int id = trigTrsh.size() - 1;
        for(auto it = trigTrsh.rbegin(); it != trigTrsh.rend(); ++it, --id) 
            if(pT0 >= it->second) {
                return make_tuple(1. / lumiMap[0].at(it->first),  1. / lumiMap[periodId].at(it->first), id);
            }

        return make_tuple(0.,0., -1);

    }

};






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

struct Histogram {
   static const int nPer = 2 + 'H'-'B'; 

   TString name, title; 
   int nDim = 0;
   vector<TH1 *> hist;


   /*
   Histogram(TString n, TString t, vector<double> bins) {
       for(int i = 0; i < nPer; ++i) {
            TString s = i;
           hist.push_back(new TH1D(n+s, t+s, bins.size()-1, bins.data()));
       }
       nDim = 1;
   }
   Histogram(TString n, TString t, int nbins, double l, double h) {
       for(int i = 0; i < nPer; ++i) {
            TString s = i;
            hist.push_back(new TH1D(n+s, t+s, nbins, l, h));
       }
       nDim = 1;
   }

   Histogram(TString n, TString t, int n1, double l1, double h1, int n2, double l2, double h2, int n3, double l3, double h3) {
       for(int i = 0; i < nPer; ++i) {
            TString s = i;
            hist.push_back(new TH3D(n+s, t+s, n1, l1, h1, n2, l2, h2, n3, l3, h3));
       }
       nDim = 3;
   }


   void Fill(int id, double x1, double x2=-99, double x3=-99, double x4=-99)
   {
       if(nDim == 1) {
           hist[0]->Fill(x1, );
           //for(int i = 0; i < n
       }
   }
   */
   /*
   int nPar = 0;
   bool isType = false;

   vector<vector<double>> vecArgs;
   //vector<double>   numArgs;
   int nDim = 0;

   template<typename T>
   void Init(TString n, TString t, T value)
   {
       Add(value);

       vector<int> types;
       for(int i = 0; i < vecArgs.size(); ++i) {
            if(vecArgs[i].size() == 1) {
                if(types.size() == 0 || types.back() == 3 || types.back() == -1)
                    types.push_back(1);
                else
                    types.back() += 1;
            }
            else {
                types.push_back(-1);
            }
       }

       for(int i = 0; i < types.size(); ++i)
           cout << "Types " << i <<" "<< types[i] << endl;
   }

   void Add(vector<double> vec) {
        vecArgs.push_back(vec);
   }
   void Add(double num) {
        vecArgs.push_back({num});
   }

   template<typename T, typename... Targs>
   void Init(TString n, TString t, T value, Targs... Fargs)
   {
       //cout << value << endl;
       Add(value);
       Init(n, t, Fargs...);

   }

   void Print() {
        for(int i = 0; i < vecArgs.size(); ++i)
            for(int j = 0; j < vecArgs[i].size(); ++j)
                cout << i<<" "<<j<<" : "<<vecArgs[i][j] << endl;
   }
   */

};

struct HistoManager {



};

void match(const char  *fileName, const char *outDir, int nevMax)
{
    //TFile *f = TFile::Open("/nfs/dust/cms/user/zlebcr/JEC/ntuples2/merged/jetsG.root");
    TFile *f = TFile::Open(fileName, "READ");
    TTreeReader chsTree("ak4/events", f);
    TTreeReader puppiTree("ak4PUPPI/events", f);


    TTreeReaderArray<float> CHSjetPt = {chsTree, "jetPt"};
    TTreeReaderArray<float> CHSjetUnc = {chsTree, "jetUnc"};
    TTreeReaderArray<float> CHSjetEta = {chsTree, "jetEta"};
    TTreeReaderArray<float> CHSjetPhi = {chsTree, "jetPhi"};
    TTreeReaderArray<float> CHSjetJEC = {chsTree, "jetJEC"};

    TTreeReaderValue<vector<bool>> CHStriggerBit = {chsTree, "triggerBit"};

    TTreeReaderArray<float> PUPPIjetPt = {puppiTree, "jetPt"};
    TTreeReaderArray<float> PUPPIjetUnc = {puppiTree, "jetUnc"};
    TTreeReaderArray<float> PUPPIjetEta = {puppiTree, "jetEta"};
    TTreeReaderArray<float> PUPPIjetPhi = {puppiTree, "jetPhi"};
    TTreeReaderArray<float> PUPPIjetJEC = {puppiTree, "jetJEC"};
    TTreeReaderValue<vector<bool>> PUPPItriggerBit = {puppiTree, "triggerBit"};

    TH1::SetDefaultSumw2();

    const vector<double> etaBins2  = { 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

    const vector<double> Ptbinning = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890,  2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717,    7000};

    const int nAsym = 60;
    vector<double> asymBins;
    for(int i = 0; i <= nAsym; ++i) asymBins.push_back( 0.8 + 0.4 * i /(nAsym+0.0));

    const int nPer = 8;
    vector<TH3D *> hBalEtaPt(nPer);
    vector<TH1D *> hJetPt(nPer);
    vector<TH3D *> hBalEtaPtUp(nPer);
    vector<TH1D *> hJetPtUp(nPer);
    vector<TH3D *> hBalEtaPtDn(nPer);
    vector<TH1D *> hJetPtDn(nPer);
    vector<TH3D *> hJECpuppi(nPer);
    vector<TH3D *> hJECchs(nPer);

    vector<TH2D *> hEtaPtCHS(nPer);
    vector<TH2D *> hEtaPtPUPPI(nPer);
    vector<TH2D *> hEtaPtPUPPIalone(nPer);


    vector<double> jecBins;
    for(double s = 0.9; s <= 1.6; s += 0.7/70)
        jecBins.push_back(s);

    for(int i = 0; i < nPer; ++i) {
        char per = 'A' + i;
        hBalEtaPt[i] = new TH3D(SF("hBalEtaPt_%c",per), SF("hBalEtaPt_%c",per), etaBins2.size()-1, etaBins2.data(), Ptbinning.size()-1, Ptbinning.data(),  asymBins.size()-1, asymBins.data()); //,  60, 0.8, 1.2);
        hJetPt[i] = new TH1D(SF("hJetPt_%c",per), SF("hJetPt_%c",per), Ptbinning.size()-1, Ptbinning.data());

        hBalEtaPtUp[i] = new TH3D(SF("hBalEtaPtUp_%c",per), SF("hBalEtaPtUp_%c",per), etaBins2.size()-1, etaBins2.data(), Ptbinning.size()-1, Ptbinning.data(),  asymBins.size()-1, asymBins.data()); //,  60, 0.8, 1.2);
        hJetPtUp[i] = new TH1D(SF("hJetPtUp_%c",per), SF("hJetPtUp_%c",per), Ptbinning.size()-1, Ptbinning.data());

        hBalEtaPtDn[i] = new TH3D(SF("hBalEtaPtDn_%c",per), SF("hBalEtaPtDn_%c",per), etaBins2.size()-1, etaBins2.data(), Ptbinning.size()-1, Ptbinning.data(),  asymBins.size()-1, asymBins.data()); //,  60, 0.8, 1.2);
        hJetPtDn[i] = new TH1D(SF("hJetPtDn_%c",per), SF("hJetPtDn_%c",per), Ptbinning.size()-1, Ptbinning.data());

        hJECpuppi[i] = new TH3D(SF("hJECpuppi_%c",per), SF("hJECpuppi_%c",per), etaBins2.size()-1, etaBins2.data(),
                                                                                Ptbinning.size()-1, Ptbinning.data(),
                                                                                jecBins.size()-1, jecBins.data());
        hJECchs[i] = new TH3D(SF("hJECchs_%c",per), SF("hJECchs_%c",per), etaBins2.size()-1, etaBins2.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data(),
                                                                          jecBins.size()-1, jecBins.data());

        hEtaPtCHS[i] =  new TH2D(SF("hEtaPtCHS_%c",per), SF("hEtaPtCHS_%c",per), etaBins2.size()-1, etaBins2.data(),Ptbinning.size()-1, Ptbinning.data() );
        hEtaPtPUPPI[i] = new TH2D(SF("hEtaPtPUPPI_%c",per), SF("hEtaPtPUPPI_%c",per), etaBins2.size()-1, etaBins2.data(),Ptbinning.size()-1, Ptbinning.data() );
        hEtaPtPUPPIalone[i] = new TH2D(SF("hEtaPtPUPPIalone_%c",per), SF("hEtaPtPUPPIalone_%c",per), etaBins2.size()-1, etaBins2.data(),Ptbinning.size()-1, Ptbinning.data() );

    }


    Luminosity lum;
    lum.LoadLumis();

    //Histogram hist("ah", "b", {1,2, 3});
    //hist.Init("ah", "b", 1, 2.2, 3.3, vector<double>{3.3, 4.3}, 5, 2.3, 5.1);
    //hist.Print();

    //return;
    const double jetSel = 5;

    int iEv = -1;
    while(chsTree.Next() && puppiTree.Next()) {
        ++iEv;

        if(iEv % 1000000 == 0 ) cout <<fileName<<" "<< iEv << endl;
        if(nevMax != -1) if(iEv >  nevMax) break;


        if(CHSjetPt.GetSize() < 1 || PUPPIjetPt.GetSize() < 1) continue;
        TString fileName = ((TChain*) chsTree.GetTree())->GetCurrentFile()->GetName();
        int fileId = fileName[fileName.Length()-6] -'A';

        double wgt, wgtTot;
        int id;
        tie(wgt,wgtTot,id) = lum.GetWeightID(fileId, CHSjetPt[0]);
        if(wgt == 0 || id < 0)  continue;

        //cout << "Info " << CHSjetPt[0]<<" "<< id <<" "<< (*CHStriggerBit).at(id) << endl;
        if((*CHStriggerBit).at(id) != 1) continue;


        hJetPt[fileId]->Fill(CHSjetPt[0], wgt);
        hJetPt[0]->Fill(CHSjetPt[0], wgtTot);

        hJetPtUp[fileId]->Fill(CHSjetPt[0]*(1+CHSjetUnc[0]), wgt);
        hJetPtUp[0]->Fill(CHSjetPt[0]*(1+CHSjetUnc[0]), wgtTot);

        hJetPtDn[fileId]->Fill(CHSjetPt[0]*(1-CHSjetUnc[0]), wgt);
        hJetPtDn[0]->Fill(CHSjetPt[0]*(1-CHSjetUnc[0]), wgtTot);

        //cout << CHSjetPt[0] << " "<< CHSjetPt[0]*(1+CHSjetUnc[0]) <<" "<< CHSjetPt[0]*(1-CHSjetUnc[0]) << endl;

        //exit(0);
        /*
        //Trigger 120
        if((*CHStriggerBit).at(3)) continue;

        if(CHSjetPt.GetSize() == 0 || PUPPIjetPt.GetSize() == 0)
            continue;
        if(CHSjetPt[0] < 160) continue;
        */

        for(int i = 0; i < PUPPIjetPt.GetSize(); ++i) {
            hJECpuppi[fileId]->Fill(abs(PUPPIjetEta[i]), PUPPIjetPt[i], PUPPIjetJEC[i], wgt);
            hJECpuppi[0]->Fill(abs(PUPPIjetEta[i]),  PUPPIjetPt[i], PUPPIjetJEC[i], wgtTot);


        }
        for(int i = 0; i < CHSjetPt.GetSize(); ++i) {
            hJECchs[fileId]->Fill(abs(CHSjetEta[i]), CHSjetPt[i], CHSjetJEC[i], wgt);
            hJECchs[0]->Fill(abs(CHSjetEta[i]), CHSjetPt[i], CHSjetJEC[i], wgtTot);

            if(CHSjetPt[i] >= jetSel) {
                hEtaPtCHS[fileId]->Fill(abs(CHSjetEta[i]), CHSjetPt[i], wgt);
                hEtaPtCHS[0]->Fill(abs(CHSjetEta[i]),  CHSjetPt[i], wgtTot);
            }
        }




        set<int> Indx;
        for(int i = 0; i < CHSjetPt.GetSize(); ++i)
            Indx.insert(i);


        for(int i = 0; i < PUPPIjetPt.GetSize(); ++i) {
            if(PUPPIjetPt[i] < jetSel) continue;

            hEtaPtPUPPI[fileId]->Fill(abs(PUPPIjetEta[i]), PUPPIjetPt[i], wgt);
            hEtaPtPUPPI[0]->Fill(abs(PUPPIjetEta[i]),  PUPPIjetPt[i], wgtTot);


            /*
            int m = -1;
            for(int j = 0; j < CHSjetPt.GetSize(); ++j) {
                double d2 = dist2(PUPPIjetEta[i], PUPPIjetPhi[i], CHSjetEta[j], CHSjetPhi[j]);
                if(d2 < 0.2*0.2) {
                    m = j;
                    break;
                }
            }
            */

            int m = -1;
            for(int ind :  Indx) {
                double d2 = dist2(PUPPIjetEta[i], PUPPIjetPhi[i], CHSjetEta[ind], CHSjetPhi[ind]);
                if(d2 < 0.2*0.2) {
                    m = ind;
                    Indx.erase(m);
                    break;
                }
            }


            if(m != -1) {
                //cout << "match " << PUPPIjetPt[i] << " "<< PUPPIjetPt[i] / CHSjetPt[m] << endl;
                double r = PUPPIjetPt[i] / CHSjetPt[m];

                double rUp = PUPPIjetPt[i] / (CHSjetPt[m] * (1 + CHSjetUnc[m]));
                double rDn = PUPPIjetPt[i] / (CHSjetPt[m] * (1 - CHSjetUnc[m]));

                hBalEtaPt[fileId]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], r, wgt);
                hBalEtaPt[0]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], r, wgtTot);

                hBalEtaPtUp[fileId]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], rUp, wgt);
                hBalEtaPtUp[0]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], rUp, wgtTot);
                hBalEtaPtDn[fileId]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], rDn, wgt);
                hBalEtaPtDn[0]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], rDn, wgtTot);

                //cout <<"Event " <<  PUPPIjetJEC[i] << " "<< CHSjetJEC[m] << endl;

            }
            else { //not matched
                hEtaPtPUPPIalone[fileId]->Fill(abs(PUPPIjetEta[i]), PUPPIjetPt[i], wgt);
                hEtaPtPUPPIalone[0]->Fill(abs(PUPPIjetEta[i]),  PUPPIjetPt[i], wgtTot);
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


        //if(CHSjetPt.GetSize() > 0 && PUPPIjetPt.GetSize() > 0)
            //cout <<iEv<<" "<< CHSjetPt[0] << " "<< PUPPIjetPt[0]  << endl;

    }
    //TFile *f = TFile::Open(fileName);
    TString outName = fileName;
    outName = "histoDir/"+TString(outDir) + "/"+outName(outName.Last('/')+1, 1000);
    cout << outName << endl;
    //return;

    TFile *fOut = TFile::Open(outName, "RECREATE");
    for(int i = 0; i < hBalEtaPt.size(); ++i) {
        hBalEtaPt[i]->Write();
        hJetPt[i]->Write();
        hBalEtaPtUp[i]->Write();
        hJetPtUp[i]->Write();
        hBalEtaPtDn[i]->Write();
        hJetPtDn[i]->Write();

        hJECchs[i]->Write();
        hJECpuppi[i]->Write();

        hEtaPtCHS[i]->Write();
        hEtaPtPUPPI[i]->Write();
        hEtaPtPUPPIalone[i]->Write();

    }
    fOut->Close();

    //TH1D *hPt = hBalEtaPt->ProjectionZ("Proj");


    //hPt->Draw();

}



