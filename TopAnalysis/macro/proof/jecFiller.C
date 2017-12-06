#define jecFiller_cxx
// The class definition in jecFiller.h has been generated automatically
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
// root> T->Process("jecFiller.C")
// root> T->Process("jecFiller.C","some options")
// root> T->Process("jecFiller.C+")
//


#include "jecFiller.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <strstream>
#include <TString.h>
#include <TFile.h>


//static const vector<int> trigger_threshold{40,60,80,140,200,260,320,450,500};
//static const size_t ntriggers = trigger_thresholds.size();
//extern const vector<int> trigger_threshold;

using namespace std;
void jecFiller::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void jecFiller::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   /*   
   vector<double> ptEdges = {18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588, 1684,1784,1890,2000,2116,2238,2366,2500, 2640,2787,2941,3103,3273,3450,3637,3832,4037,4252,4477,4713,4961,5220,5492,5777,6076,6389,6717,7000 };

   hJetPt  = new TH1D("hJetPt", "hist", ptEdges.size()-1, ptEdges.data());
   fOutput->Add(hJetPt);
   */
   //vector<double> JetPt {40.,60.,80.,140.,200.,260.,320.,400.,450.,500};
   double JetPt[10]={40.,60.,80.,140.,200.,260.,320.,400.,450.,500};
   double JetPtThreshold[10]={60.,80.,140.,200.,260.,320.,400.,450.,500.,600.};

   // assert(JetPt.size() == JetPtThreshold.size());

   histoEmulated = new TH2D[10];
   histoPtAll = new TH2D[10];
   
#define SF TString::Format

   //vector<TH2 *> histoPtAll(JetPt.size(),nullptr);
   //for (size_t i = 0; i < histoPtAll.size(); ++i) 
   //    histoPtAll[i] = new TH2D(SF("Jindrich%d", 42), "title", NBins, ...);

   //for (TH2 * h: histoPtAll) {
   //    const char * str = TString::Format("Jindrich%d", 42);
   //    h = new TH2D(...);
   //}

   stringstream TitEmulated, TitPtAll, TitCurrent;
   TitEmulated.str("");
   TitPtAll.str("");
   TitCurrent.str("");
   TitEmulated<<"Emulated";
   TitPtAll<<"All";
   
   int NBins = 2800;
   double PtMin = 0.;
   double PtMax = 7000.;
   int EtaBin = 50;
   double EtaMin = 0.;
   double EtaMax = 5.;

   for(int i = 0; i<10; ++i){
     TitCurrent.str("");
     TitCurrent<<"HLTJetPt"<<JetPt[i+1]<<"_"<<TitEmulated.str();
     histoEmulated[i].SetName(TitCurrent.str().c_str());
     histoEmulated[i].SetTitle(TitCurrent.str().c_str());
     histoEmulated[i].SetBins(NBins,PtMin, PtMax, EtaBin, EtaMin, EtaMax);
     histoEmulated[i].SetDirectory(0);
     histoEmulated[i].Sumw2();
     // TH1::SetDefaultSumw2(kTRUE);

     TitCurrent.str("");
     TitCurrent<<"HLTJetPt"<<JetPt[i]<<"_"<<TitPtAll.str();
     histoPtAll[i].SetName(TitCurrent.str().c_str());
     histoPtAll[i].SetTitle(TitCurrent.str().c_str());
     histoPtAll[i].SetBins(NBins,PtMin,PtMax,EtaBin, EtaMin, EtaMax);
     histoPtAll[i].SetDirectory(0);
     histoPtAll[i].Sumw2();
  
     fOutput->Add(&histoEmulated[i]);
     fOutput->Add(&histoPtAll[i]);
   }
}


Bool_t jecFiller::Process(Long64_t entry)
{
  double JetPtThreshold[10]={60.,80.,140., 200., 260., 320., 400., 450., 500., 600.};
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
   
   //for(unsigned i = 0; i < jetPt.GetSize(); ++i)
   //    hJetPt->Fill(jetPt[i]);
   
   double deltaPhiMatch = 0.;
   double hweight = 1.;
   //bool hltcut[10];
   int etabin = 0.;


   //vector<bool> hltcuts(ntriggers,false);
   for(int triggerID = 0; triggerID < 10; ++triggerID){

       //for(int i = 0; i<10;++i)hltcut[i] = false;
       bool hltcut = false;

       if(triggerBit->at(triggerID) == 0) continue;

       size_t it = max_element(jetPt.begin(), jetPt.end()).fIndex;
       double leading_pT = jetPt[it];
       double leading_Eta = jetEta[it];
       double leading_Phi = jetPhi[it];

       for(int hltiobj = 0; hltiobj < *nTriggerObjects ; ++hltiobj){
           deltaPhiMatch = abs(HLTjetPhi[hltiobj]-leading_Phi);
           deltaPhiMatch = min(2*M_PI-deltaPhiMatch, deltaPhiMatch);


           //deltaPhiMatch = HLTjetPhi[hltiobj]-leading_Phi;
           //if(deltaPhiMatch < -TMath::Pi())deltaPhiMatch = deltaPhiMatch+2*TMath::Pi();
           //if(deltaPhiMatch > +TMath::Pi())deltaPhiMatch = deltaPhiMatch-2*TMath::Pi();
           //deltaPhiMatch = fabs(deltaPhiMatch);
           hltcut =  hltcut || HLTjetPt[hltiobj] > JetPtThreshold[triggerID];
       }//hltiobj

       histoPtAll[triggerID].Fill(leading_pT,fabs(leading_Eta),hweight);
       if(hltcut)histoEmulated[triggerID].Fill(leading_pT,fabs(leading_Eta),hweight);
   }//triggerID
   
   return kTRUE;
}

void jecFiller::SlaveTerminate ()
{
  //SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void jecFiller::Terminate ()
{
  //    hJetPt  = dynamic_cast<TH1D*>(fOutput->FindObject("hJetPt"));
  //  assert(hJetPt);
  //  TCanvas *can = new TCanvas("can", "can");
  //  hJetPt->Draw();
  //  can->SaveAs("radek.pdf");
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
    
  TString outputFilename;
  outputFilename = "Data2016test.root";
  TFile* outputFile = new TFile(outputFilename,"recreate");
  for(int i = 0; i < 10; ++i){
    histoPtAll[i].Write();
    histoEmulated[i].Write();
  }
  
  outputFile->Close();

}
