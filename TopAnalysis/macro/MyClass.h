//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Nov 29 12:13:20 2017 by ROOT version 6.06/01
// from TTree events/events
// found on file: jetCorr1.root
//////////////////////////////////////////////////////////

#ifndef MyClass_h
#define MyClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"

class MyClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxrun = 1;
   const Int_t kMaxevt = 1;
   const Int_t kMaxlumi = 1;
   const Int_t kMaxnVtx = 1;
   const Int_t kMaxnJets = 1;
   const Int_t kMaxpvRho = 1;
   const Int_t kMaxpvz = 1;
   const Int_t kMaxpvchi2 = 1;
   const Int_t kMaxpvndof = 1;
   const Int_t kMaxrho = 1;
   const Int_t kMaxht = 1;
   const Int_t kMaxmetEt = 1;
   const Int_t kMaxmetSumEt = 1;
   const Int_t kMaxmetEtNoHF = 1;
   const Int_t kMaxmetSumEtNoHF = 1;
   const Int_t kMaxmetEtPuppi = 1;
   const Int_t kMaxmetSumEtPuppi = 1;
   const Int_t kMaxnTriggerObjects = 1;

   // Declaration of leaf types
   Int_t           runNo;
   Int_t           evtNo;
   Int_t           lumi;
   Int_t           nvtx;
   Int_t           nJets;
   Float_t         pvRho;
   Float_t         pvz;
   Float_t         pvchi2;
   Float_t         pvndof;
   Float_t         rho;
   Float_t         ht;
   Float_t         metEt;
   Float_t         metSumEt;
   Float_t         metEtNoHF;
   Float_t         metSumEtNoHF;
   Float_t         metEtPuppi;
   Float_t         metSumEtPuppi;
   vector<float>   *jetPt;
   vector<float>   *jetUnc;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetMass;
   vector<float>   *jetEnergy;
   vector<float>   *jetChf;
   vector<float>   *jetNhf;
   vector<float>   *jetPhf;
   vector<float>   *jetMuf;
   vector<float>   *jetElf;
   vector<int>     *jetChm;
   vector<int>     *jetNhm;
   vector<int>     *jetPhm;
   vector<int>     *jetElm;
   vector<int>     *jetMum;
   vector<bool>    *jetIsBtag;
   vector<float>   *HLTjetPt;
   vector<float>   *HLTjetEta;
   vector<float>   *HLTjetPhi;
   vector<float>   *HLTjetMass;
   vector<bool>    *triggerBit;
   vector<int>     *triggerPre;
   Int_t           nTriggerObject;

   // List of branches
   TBranch        *b_run_;   //!
   TBranch        *b_evt_;   //!
   TBranch        *b_lumi_;   //!
   TBranch        *b_nVtx_;   //!
   TBranch        *b_nJets_;   //!
   TBranch        *b_pvRho_;   //!
   TBranch        *b_pvz_;   //!
   TBranch        *b_pvchi2_;   //!
   TBranch        *b_pvndof_;   //!
   TBranch        *b_rho_;   //!
   TBranch        *b_ht_;   //!
   TBranch        *b_metEt_;   //!
   TBranch        *b_metSumEt_;   //!
   TBranch        *b_metEtNoHF_;   //!
   TBranch        *b_metSumEtNoHF_;   //!
   TBranch        *b_metEtPuppi_;   //!
   TBranch        *b_metSumEtPuppi_;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetUnc;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetEnergy;   //!
   TBranch        *b_jetChf;   //!
   TBranch        *b_jetNhf;   //!
   TBranch        *b_jetPhf;   //!
   TBranch        *b_jetMuf;   //!
   TBranch        *b_jetElf;   //!
   TBranch        *b_jetChm;   //!
   TBranch        *b_jetNhm;   //!
   TBranch        *b_jetPhm;   //!
   TBranch        *b_jetElm;   //!
   TBranch        *b_jetMum;   //!
   TBranch        *b_jetIsBtag;   //!
   TBranch        *b_HLTjetPt;   //!
   TBranch        *b_HLTjetEta;   //!
   TBranch        *b_HLTjetPhi;   //!
   TBranch        *b_HLTjetMass;   //!
   TBranch        *b_triggerBit;   //!
   TBranch        *b_triggerPre;   //!
   TBranch        *b_nTriggerObjects_;   //!

   MyClass(TTree *tree=0);
   virtual ~MyClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MyClass_cxx
MyClass::MyClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("runG.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("runG.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("runG.root:/ak4");
      dir->GetObject("events",tree);

   }
   Init(tree);
}

MyClass::~MyClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MyClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MyClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MyClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jetPt = 0;
   jetUnc = 0;
   jetEta = 0;
   jetPhi = 0;
   jetMass = 0;
   jetEnergy = 0;
   jetChf = 0;
   jetNhf = 0;
   jetPhf = 0;
   jetMuf = 0;
   jetElf = 0;
   jetChm = 0;
   jetNhm = 0;
   jetPhm = 0;
   jetElm = 0;
   jetMum = 0;
   jetIsBtag = 0;
   HLTjetPt = 0;
   HLTjetEta = 0;
   HLTjetPhi = 0;
   HLTjetMass = 0;
   triggerBit = 0;
   triggerPre = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNo", &runNo, &b_run_);
   fChain->SetBranchAddress("evtNo", &evtNo, &b_evt_);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi_);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nVtx_);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets_);
   fChain->SetBranchAddress("pvRho", &pvRho, &b_pvRho_);
   fChain->SetBranchAddress("pvz", &pvz, &b_pvz_);
   fChain->SetBranchAddress("pvchi2", &pvchi2, &b_pvchi2_);
   fChain->SetBranchAddress("pvndof", &pvndof, &b_pvndof_);
   fChain->SetBranchAddress("rho", &rho, &b_rho_);
   fChain->SetBranchAddress("ht", &ht, &b_ht_);
   fChain->SetBranchAddress("metEt", &metEt, &b_metEt_);
   fChain->SetBranchAddress("metSumEt", &metSumEt, &b_metSumEt_);
   fChain->SetBranchAddress("metEtNoHF", &metEtNoHF, &b_metEtNoHF_);
   fChain->SetBranchAddress("metSumEtNoHF", &metSumEtNoHF, &b_metSumEtNoHF_);
   fChain->SetBranchAddress("metEtPuppi", &metEtPuppi, &b_metEtPuppi_);
   fChain->SetBranchAddress("metSumEtPuppi", &metSumEtPuppi, &b_metSumEtPuppi_);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetUnc", &jetUnc, &b_jetUnc);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetEnergy", &jetEnergy, &b_jetEnergy);
   fChain->SetBranchAddress("jetChf", &jetChf, &b_jetChf);
   fChain->SetBranchAddress("jetNhf", &jetNhf, &b_jetNhf);
   fChain->SetBranchAddress("jetPhf", &jetPhf, &b_jetPhf);
   fChain->SetBranchAddress("jetMuf", &jetMuf, &b_jetMuf);
   fChain->SetBranchAddress("jetElf", &jetElf, &b_jetElf);
   fChain->SetBranchAddress("jetChm", &jetChm, &b_jetChm);
   fChain->SetBranchAddress("jetNhm", &jetNhm, &b_jetNhm);
   fChain->SetBranchAddress("jetPhm", &jetPhm, &b_jetPhm);
   fChain->SetBranchAddress("jetElm", &jetElm, &b_jetElm);
   fChain->SetBranchAddress("jetMum", &jetMum, &b_jetMum);
   fChain->SetBranchAddress("jetIsBtag", &jetIsBtag, &b_jetIsBtag);
   fChain->SetBranchAddress("HLTjetPt", &HLTjetPt, &b_HLTjetPt);
   fChain->SetBranchAddress("HLTjetEta", &HLTjetEta, &b_HLTjetEta);
   fChain->SetBranchAddress("HLTjetPhi", &HLTjetPhi, &b_HLTjetPhi);
   fChain->SetBranchAddress("HLTjetMass", &HLTjetMass, &b_HLTjetMass);
   fChain->SetBranchAddress("triggerBit", &triggerBit, &b_triggerBit);
   fChain->SetBranchAddress("triggerPre", &triggerPre, &b_triggerPre);
   fChain->SetBranchAddress("nTriggerObject", &nTriggerObject, &b_nTriggerObjects_);
   Notify();
}

Bool_t MyClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MyClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MyClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MyClass_cxx
