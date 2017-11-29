#define MyClass_cxx
#include "MyClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <strstream>
#include "TString.h"
#include "TFile.h"

using namespace std;
void MyClass::Loop()
{
//   In a ROOT session, you can do:
//      root> .L MyClass.C
//      root> MyClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   
   double JetPtThreshold[10] = {60.,80.,140.,200.,260.,320.,400.,450.,500.,600.};
   double JetPt[10] = {40., 60., 80., 140., 200., 260., 320., 400., 450., 500.};
   TH2F *histoEmulated = new TH2F[10];
   TH2F *histoAll = new TH2F[10];

   stringstream TitEmulated, TitAll, TitCurrent;
   TitEmulated.str("");
   TitAll.str("");
   TitCurrent.str("");
   TitEmulated<<"Emulated";
   TitAll<<"All";
   
   int Ptbins = 80;
   double Ptbinning[81] = {0, 1, 5, 6, 8, 10, 12, 15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84,97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};

   int Etabins = 47;
   double Etabinning[48];
   double leading_pT = 0.;
   double leading_Eta = 0.;
   double leading_phi = 0.;
   double deltaPhiMatch = 0.;
   double hweight = 1.;
   bool hltcut[10];
   for(int i = 0; i<(Etabins+1);i++)Etabinning[i]=i*0.1;
   for(int i = 0; i < 10; i++)hltcut[i]=false;
   for(int i = 0; i<10;i++){
     TitCurrent.str("");
     TitCurrent<<"HLTJetPt"<<JetPt[i]<<"_"<<TitEmulated.str();
     histoEmulated[i].SetName(TitCurrent.str().c_str());
     histoEmulated[i].SetTitle(TitCurrent.str().c_str());
     histoEmulated[i].SetBins(Ptbins, Ptbinning, Etabins, Etabinning);
     histoEmulated[i].SetDirectory(0);
     histoEmulated[i].Sumw2();

     TitCurrent.str("");
     TitCurrent<<"HLTJetPt"<<JetPt[i]<<"_"<<TitAll.str();
     histoAll[i].SetName(TitCurrent.str().c_str());
     histoAll[i].SetTitle(TitCurrent.str().c_str());
     histoAll[i].SetBins(Ptbins, Ptbinning, Etabins, Etabinning);
     histoAll[i].SetDirectory(0);
     histoAll[i].Sumw2();
   }

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;
   int decade = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //--- Progress report ---
     double progress = 100.0*jentry/(1.0*nentries);
     int a = TMath::FloorNint(progress);
     if( a > decade)
       cout<<a<<" %"<<endl;
     decade = a;
     //--- Progress report ---
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      for(int triggerID = 0; triggerID < 10; triggerID++){
	hltcut[triggerID] = false;
	if(triggerBit->at(triggerID)){
	  for (int i = 0; i < nJets; i++){
	    if(i==0){
	      leading_pT = jetPt->at(i);
	      leading_Eta = jetEta->at(i);
	      leading_Phi = jetPhi->at(i);
	    }
	    if( i>=1 && leading_pT < jetPt->at(i)){
	      leading_pT = jetPt->at(i);
	      leading_Eta = jetEta->at(i);
	      leading_Phi = jetPhi->at(i);
	    }
	  }// loop for nJets
	  
	  for(int hltiobj = 0; hltiobj < nTriggerObject; hltiobj++){
	    deltaPhiMatch =HLTjetPhi->at(hltiobj) - leading_Phi;
	    if(deltaPhiMatch < -TMath::Pi()) deltaPhiMatch = deltaPhiMatch+2*TMath::Pi();
	    if(deltaPhiMatch >  TMath::Pi()) deltaPhiMatch = deltaPhiMatch-2*TMath::Pi();
	    deltaPhiMatch = fabs(deltaPhiMatch);
	    if(HLTjetPt->at(hltiobj) > JetPtThreshold[triggerID] && sqrt(pow(HLTjetEta->at(hltiobj)-leading_Eta,2.)+pow(deltaPhiMatch,2))<0.2){
	      hltcut[triggerID] = true;
	    }
	    histoAll[triggerID].Fill(leading_pT,fabs(leading_Eta),hweight);
	    if(hltcut[triggerID])histoEmulated[triggerID].Fill(leading_pT,fabs(leading_Eta),hweight);
	  }// hltiobj
	} // if(bit->at(triggerID)
      } //loop for triggerId
   }

   TString outputFilename;
   outputFilename = "Data2016output.root";
   TFile* outputFile = new TFile(outputFilename,"recreate");
   for(int i = 0; i<10; i++){
     histoAll[i].Write();
     histoEmulated[i].Write();
   }
   outputFile->Close();
}
