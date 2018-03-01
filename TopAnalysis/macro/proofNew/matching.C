#define matching_cxx
// The class definition in matching.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("matching.C")
// root> T->Process("matching.C","some options")
// root> T->Process("matching.C+")
//


#include "matching.h"
#include <TH2.h>
#include <TStyle.h>


double dist2(double eta1, double phi1, double eta2, double phi2)
{
    double dPhi = abs(phi2 - phi1);
    dPhi = min(dPhi, 2*M_PI - dPhi);

    return ( (eta2-eta1)*(eta2-eta1) + dPhi*dPhi);
}


const string jecTagCHS = "Summer16_07Aug2017";
const int versionCHS = 5;

const string jecTagPUPPI = "Summer16_07Aug2017";
const int versionPUPPI = 5;


//#ifdef __MAKECINT__
//#pragma link C++ class vector<Jet>+;
//#endif

//gInterpreter->GenerateDictionary("Jet","match.h");
//#include "match.h"


//#pragma link C++ class vector<Jet> +;

struct Histogram {
   static const int nPer = 2 + 'H'-'B'; 

   TString name, title; 
   int nDim = 0;
   vector<TH1 *> hist;

};

struct HistoManager {



};

#define SF TString::Format 





void matching::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void matching::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

    TH1::SetDefaultSumw2();

    const vector<double> etaBins2  = { 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191};

    const vector<double> Ptbinning = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890,  2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717,    7000};

    const int nAsym = 60;
    vector<double> asymBins;
    for(int i = 0; i <= nAsym; ++i) asymBins.push_back( 0.8 + 0.4 * i /(nAsym+0.0));





    vector<double> jecBins;
    for(double s = 0.9; s <= 1.6; s += 0.7/70)
        jecBins.push_back(s);

    for(int i = 0; i < nPer; ++i) {
        char per = 'A' + i;
        hBalEtaPt[i] = new TH3D(SF("hBalEtaPt_%c",per), SF("hBalEtaPt_%c",per), etaBins2.size()-1, etaBins2.data(), Ptbinning.size()-1, Ptbinning.data(),  asymBins.size()-1, asymBins.data()); //,  60, 0.8, 1.2);
        fOutput->Add(hBalEtaPt[i]);
        hJetPt[i] = new TH1D(SF("hJetPt_%c",per), SF("hJetPt_%c",per), Ptbinning.size()-1, Ptbinning.data());
        fOutput->Add(hJetPt[i]);

        hJECpuppi[i] = new TH3D(SF("hJECpuppi_%c",per), SF("hJECpuppi_%c",per), etaBins2.size()-1, etaBins2.data(),
                                                                                Ptbinning.size()-1, Ptbinning.data(),
                                                                                jecBins.size()-1, jecBins.data());
        fOutput->Add(hJECpuppi[i]);
        hJECchs[i] = new TH3D(SF("hJECchs_%c",per), SF("hJECchs_%c",per), etaBins2.size()-1, etaBins2.data(),
                                                                          Ptbinning.size()-1, Ptbinning.data(),
                                                                          jecBins.size()-1, jecBins.data());

        fOutput->Add(hJECchs[i]);
        hEtaPtCHS[i] =  new TH2D(SF("hEtaPtCHS_%c",per), SF("hEtaPtCHS_%c",per), etaBins2.size()-1, etaBins2.data(),Ptbinning.size()-1, Ptbinning.data() );
        fOutput->Add(hEtaPtCHS[i]);
        hEtaPtPUPPI[i] = new TH2D(SF("hEtaPtPUPPI_%c",per), SF("hEtaPtPUPPI_%c",per), etaBins2.size()-1, etaBins2.data(),Ptbinning.size()-1, Ptbinning.data() );
        fOutput->Add(hEtaPtPUPPI[i]);
        hEtaPtPUPPIalone[i] = new TH2D(SF("hEtaPtPUPPIalone_%c",per), SF("hEtaPtPUPPIalone_%c",per), etaBins2.size()-1, etaBins2.data(),Ptbinning.size()-1, Ptbinning.data() );
        fOutput->Add(hEtaPtPUPPIalone[i]);

    }

    lum.LoadLumis();

    //string jecTag = "Summer16_07Aug2017";
    //int version = 4;
    char period = 'B';
    string jetType = "AK4PFchs";
    vector<string> dumy;
    bool isMC = false;

    jetCorrsCHS.Init(isMC, jecTagCHS, period, versionCHS, "AK4PFchs", "", dumy);
    jetCorrsPUPPI.Init(isMC, jecTagPUPPI, period, versionPUPPI, "AK4PFPuppi", "", dumy);




}



/*
void CorrectJets(TTreeReaderValue<std::vector<QCDjet> > &Jets, double rho, JECs &jetCorrs ) {
    for (unsigned i = 0; i < Jets->size(); ++i) {
        double pt  = Jets->at(i).p4.Pt();
        double eta = Jets->at(i).p4.Eta();
        double area = Jets->at(i).area;

        cout <<"R "<<  i <<" "<< pt << " "<< eta <<" "<< area <<" "<< rho <<  endl;
        double CorFactor, Unc;
        vector<string> dumy;

        double corr = jetCorrs.JEC_CHScorrections(pt, eta, area, rho, dumy, CorFactor, Unc);
        cout << "Corr " << corr << endl;
    }
}
*/











Bool_t matching::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);


   //CorrectJets(chsJets, *rho, jetCorrsCHS);
   //CorrectJets(puppiJets, *rho, jetCorrsPUPPI);


    for (unsigned i = 0; i < chsJets->size(); ++i) {
        double pt  = chsJets->at(i).p4.Pt();
        double eta = chsJets->at(i).p4.Eta();
        double area = chsJets->at(i).area;

        //cout <<"R "<<  i <<" "<< pt << " "<< eta <<" "<< area <<" "<< *rho <<  endl;
        double CorFactor, Unc;
        vector<string> dumy;

        double corr = jetCorrsCHS.JEC_CHScorrections(pt, eta, area, *rho, dumy, CorFactor, Unc);
        chsJets->at(i).p4.SetPt(pt * corr);
        //cout << "Corr " << corr << endl;
    }
    sort(chsJets->begin(), chsJets->end(), [](const QCDjet &a, const QCDjet &b) {return a.p4.Pt() > b.p4.Pt();});

    //for (unsigned i = 0; i < chsJets->size(); ++i) {
        //cout << i <<" "<< chsJets->at(i).p4.Pt() << endl;
    //}

    for (unsigned i = 0; i < puppiJets->size(); ++i) {
        double pt  = puppiJets->at(i).p4.Pt();
        double eta = puppiJets->at(i).p4.Eta();
        double area = puppiJets->at(i).area;

        //cout <<"R "<<  i <<" "<< pt << " "<< eta <<" "<< area <<" "<< *rho <<  endl;
        double CorFactor, Unc;
        vector<string> dumy;

        double corr = jetCorrsPUPPI.JEC_CHScorrections(pt, eta, area, *rho, dumy, CorFactor, Unc);
        puppiJets->at(i).p4.SetPt(pt * corr);
    }
    sort(puppiJets->begin(), puppiJets->end(), [](const QCDjet &a, const QCDjet &b) {return a.p4.Pt() > b.p4.Pt();});










    const double jetSel = 5;

    /*
    for (unsigned i = 0; i < puppiJets->size(); ++i) {
        cout << i << " b "<< puppiJets->at(i).p4.Pt() <<" "<< puppiJets->at(i).jetJECtot <<  endl;
        //puppiJets->at(i).p4.SetPt(puppiJets->at(i).p4.Pt() * puppiJets->at(i).jetJECtot);
        //cout << i << " a "<< puppiJets->at(i).p4.Pt() << endl;
    }
    */


   //cout << *metEtPuppi <<" "<< chsJets->size() << endl;

    if(chsJets->size() < 1 || puppiJets->size() < 1) return 0;
    TString fileName =  fReader.GetTree()->GetCurrentFile()->GetName();
    //cout << "fileName is " << fileName  << endl;
        int fileId = fileName[fileName.Length()-6] -'A';

        double wgt, wgtTot;
        int id;
        tie(wgt,wgtTot,id) = lum.GetWeightID(fileId, chsJets->at(0).p4.Pt() );
        if(wgt == 0 || id < 0)  return 0;

        //cout << "Info " << CHSjetPt[0]<<" "<< id <<" "<< (*CHStriggerBit).at(id) << endl;
        if(triggerBit->at(id) != 1) return 0;


        hJetPt[fileId]->Fill(chsJets->at(0).p4.Pt(), wgt);
        hJetPt[0]->Fill(chsJets->at(0).p4.Pt(), wgtTot);

        //cout << CHSjetPt[0] << " "<< CHSjetPt[0]*(1+CHSjetUnc[0]) <<" "<< CHSjetPt[0]*(1-CHSjetUnc[0]) << endl;

        //exit(0);
        /*
        //Trigger 120
        if((*CHStriggerBit).at(3)) continue;

        if(CHSjetPt.GetSize() == 0 || PUPPIjetPt.GetSize() == 0)
            continue;
        if(CHSjetPt[0] < 160) continue;
        */

        for (unsigned i = 0; i < puppiJets->size(); ++i) {
            hJECpuppi[fileId]->Fill(abs(puppiJets->at(i).p4.Eta()), puppiJets->at(i).p4.Pt(), 1, wgt);
            hJECpuppi[0]->Fill(abs(puppiJets->at(i).p4.Eta()),  puppiJets->at(i).p4.Pt(), 1, wgtTot);
        }
        for (unsigned i = 0; i < chsJets->size(); ++i) {
            hJECchs[fileId]->Fill(abs(chsJets->at(i).p4.Eta()), chsJets->at(i).p4.Pt(), 1, wgt);
            hJECchs[0]->Fill(abs(chsJets->at(i).p4.Eta()),  chsJets->at(i).p4.Pt(), 1, wgtTot);

            if (chsJets->at(i).p4.Pt() >= jetSel) {
                hEtaPtCHS[fileId]->Fill(abs(chsJets->at(i).p4.Eta()), chsJets->at(i).p4.Pt(), wgt);
                hEtaPtCHS[0]->Fill(abs(chsJets->at(i).p4.Eta()), chsJets->at(i).p4.Pt(), wgtTot);
            }
        }

        set<int> Indx;
        for(unsigned i = 0; i < chsJets->size(); ++i)
            Indx.insert(i);


        for(unsigned i = 0; i < puppiJets->size(); ++i) {
            if(puppiJets->at(i).p4.Pt() < jetSel) continue;

            hEtaPtPUPPI[fileId]->Fill(abs(puppiJets->at(i).p4.Eta()), puppiJets->at(i).p4.Pt(), wgt);
            hEtaPtPUPPI[0]->Fill(abs(puppiJets->at(i).p4.Eta()), puppiJets->at(i).p4.Pt(), wgtTot);


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
                double d2 = dist2(puppiJets->at(i).p4.Eta(), puppiJets->at(i).p4.Phi(), chsJets->at(ind).p4.Eta(), chsJets->at(ind).p4.Phi());
                if(d2 < 0.2*0.2) {
                    m = ind;
                    Indx.erase(m);
                    break;
                }
            }


            if(m != -1) {
                //cout << "match " << PUPPIjetPt[i] << " "<< PUPPIjetPt[i] / CHSjetPt[m] << endl;
                double r = puppiJets->at(i).p4.Pt() / chsJets->at(m).p4.Pt();

                //double rUp = PUPPIjetPt[i] / (CHSjetPt[m] * (1 + CHSjetUnc[m]));
                //double rDn = PUPPIjetPt[i] / (CHSjetPt[m] * (1 - CHSjetUnc[m]));

                hBalEtaPt[fileId]->Fill(abs(puppiJets->at(i).p4.Eta()), puppiJets->at(i).p4.Pt(), r, wgt);
                hBalEtaPt[0]->Fill(abs(puppiJets->at(i).p4.Eta()), puppiJets->at(i).p4.Pt() , r, wgtTot);

                //hBalEtaPtUp[fileId]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], rUp, wgt);
                //hBalEtaPtUp[0]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], rUp, wgtTot);
                //hBalEtaPtDn[fileId]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], rDn, wgt);
                //hBalEtaPtDn[0]->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], rDn, wgtTot);

                //cout <<"Event " <<  PUPPIjetJEC[i] << " "<< CHSjetJEC[m] << endl;

            }
            else { //not matched
                hEtaPtPUPPIalone[fileId]->Fill(abs(puppiJets->at(i).p4.Eta()), puppiJets->at(i).p4.Pt(), wgt);
                hEtaPtPUPPIalone[0]->Fill(abs(puppiJets->at(i).p4.Eta()), puppiJets->at(i).p4.Pt(), wgtTot);
            }



            //else
                //hBalEtaPt->Fill(abs(PUPPIjetEta[i]),PUPPIjetPt[i], 1.199);

                //cout << "NO    " << PUPPIjetPt[i] << " "<< endl;

            //cout << i << " "<< PUPPIjetPt[i] << endl;

            //if(i
        }






   return kTRUE;
}

void matching::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void matching::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.


    cout << "Saving histos" << endl;
    //return;

    TString outputFilename = SF("histos/%sV%d__%sV%d.root", jecTagCHS.c_str(), versionCHS,    jecTagPUPPI.c_str(), versionPUPPI);
    TFile* outputFile = new TFile(outputFilename,"recreate");

    for(int i = 0; i < nPer; ++i) {
        char per = 'A' + i;
        //cout << "Before " << hBalEtaPt[i] << endl;
        hBalEtaPt[i] = dynamic_cast<TH3D*>(fOutput->FindObject(SF("hBalEtaPt_%c",per)));
        assert(hBalEtaPt[i]);
        //cout << "after " << hBalEtaPt[i] << endl;
        hJetPt[i] = dynamic_cast<TH1D*>(fOutput->FindObject(SF("hJetPt_%c",per)));
        assert(hJetPt[i]);

        hJECpuppi[i] = dynamic_cast<TH3D*>(fOutput->FindObject(SF("hJECpuppi_%c",per)));
        hJECchs[i] = dynamic_cast<TH3D*>(fOutput->FindObject(SF("hJECchs_%c",per)));

        hEtaPtCHS[i] =  dynamic_cast<TH2D*>(fOutput->FindObject(SF("hEtaPtCHS_%c",per)));
        hEtaPtPUPPI[i] = dynamic_cast<TH2D*>(fOutput->FindObject(SF("hEtaPtPUPPI_%c",per)));
        hEtaPtPUPPIalone[i] = dynamic_cast<TH2D*>(fOutput->FindObject(SF("hEtaPtPUPPIalone_%c",per)));


        hBalEtaPt[i]->Write();
        hJetPt[i]->Write();
                            
        hJECpuppi[i]->Write();
        hJECchs[i]->Write();
                            
        hEtaPtCHS[i]->Write();
        hEtaPtPUPPI[i]->Write();
        hEtaPtPUPPIalone[i]->Write();
    }

    outputFile->Close();




}
