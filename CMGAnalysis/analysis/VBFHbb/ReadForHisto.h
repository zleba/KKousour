//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 16 13:44:10 2013 by ROOT version 5.32/00
// from TTree events/events
// found on file: ../../prod/VBFHbb/config/flatTree.root
//////////////////////////////////////////////////////////

#ifndef ReadForHisto_h
#define ReadForHisto_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TString.h"
// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxrun = 1;
const Int_t kMaxevt = 1;
const Int_t kMaxlumi = 1;
const Int_t kMaxnVtx = 1;
const Int_t kMaxselNOM = 1;
const Int_t kMaxselVBF = 1;
const Int_t kMaxselNOMsoft = 1;
const Int_t kMaxselVBFsoft = 1;
const Int_t kMaxnSoftJets = 1;
const Int_t kMaxnJets = 1;
const Int_t kMaxnBJets = 1;
const Int_t kMaxb1 = 1;
const Int_t kMaxb2 = 1;
const Int_t kMaxq1 = 1;
const Int_t kMaxq2 = 1;
const Int_t kMaxmvaNOM = 1;
const Int_t kMaxmvaVBF = 1;
const Int_t kMaxsoftHt = 1;
const Int_t kMaxrho = 1;
const Int_t kMaxht = 1;
const Int_t kMaxmet = 1;
const Int_t kMaxmetPhi = 1;
const Int_t kMaxmetSig = 1;
const Int_t kMaxsphericity = 1;
const Int_t kMaxaplanarity = 1;
const Int_t kMaxdEtaMax = 1;
const Int_t kMaxx1 = 1;
const Int_t kMaxx2 = 1;
const Int_t kMaxcosTheta = 1;
const Int_t kMaxmqq = 1;
const Int_t kMaxmbb = 1;
const Int_t kMaxmbbReg = 1;
const Int_t kMaxdEtaqq = 1;
const Int_t kMaxdEtabb = 1;
const Int_t kMaxptbb = 1;
const Int_t kMaxdPhiqq = 1;
const Int_t kMaxdPhibb = 1;
const Int_t kMaxyBoostqq = 1;
const Int_t kMaxyStarqq = 1;
const Int_t kMaxnpu = 1;

class ReadForHisto {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           runNo;
   Int_t           evtNo;
   Int_t           lumi;
   Int_t           nvtx;
   Bool_t          selNOM[3];
   Bool_t          selVBF[3];
   Bool_t          selNOMsoft[3];
   Bool_t          selVBFsoft[3];
   Int_t           nSoftJets;
   Int_t           nSoftJets2;
   Int_t           nSoftJets5;
   Int_t           nJets;
   Int_t           nBJets;
   Int_t           b1[3];
   Int_t           b2[3];
   Int_t           q1[3];
   Int_t           q2[3];
   Float_t         mvaNOM;
   Float_t         mvaVBF;
   Float_t         softHt;
   Float_t         rho;
   Float_t         ht;
   Float_t         met;
   Float_t         metPhi;
   Float_t         metSig;
   Float_t         sphericity;
   Float_t         aplanarity;
   Float_t         dEtaMax;
   Float_t         x1;
   Float_t         x2;
   Float_t         cosTheta[3];
   Float_t         mqq[3];
   Float_t         mbb[3];
   Float_t         mbbReg[3];
   Float_t         dEtaqq[3];
   Float_t         dEtabb[3];
   Float_t         ptbb[3];
   Float_t         dPhiqq[3];
   Float_t         dPhibb[3];
   Float_t         yBoostqq[3];
   Float_t         yStarqq[3];
   vector<int>     *btagIdx;
   vector<int>     *blikNOMIdx;
   vector<int>     *blikVBFIdx;
   vector<int>     *etaIdx;
   vector<bool>    *jetPuIdL;
   vector<bool>    *jetPuIdM;
   vector<bool>    *jetPuIdT;
   vector<bool>    *jetIdL;
   vector<bool>    *jetIdM;
   vector<bool>    *jetIdT;
   vector<bool>    *jetBtagL;
   vector<bool>    *jetBtagM;
   vector<bool>    *jetBtagT;
   vector<float>   *jetPt;
   vector<float>   *jetBtag;
   vector<float>   *jetBlikNOM;
   vector<float>   *jetBlikVBF;
   vector<float>   *jetCSV;
   vector<float>   *jetPuMva;
   vector<float>   *jetJec;
   vector<float>   *jetRegPt;
   vector<float>   *jetRegE;
   vector<float>   *jetUnc;
   vector<float>   *jetQGL;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetMetPhi;
   vector<float>   *jetMass;
   vector<float>   *jetChf;
   vector<float>   *jetNhf;
   vector<float>   *jetPhf;
   vector<float>   *jetMuf;
   vector<float>   *jetElf;
   vector<bool>    *triggerResult;
   Int_t           npu;
   Float_t         puWt;

   // List of branches
   TBranch        *b_run_;   //!
   TBranch        *b_evt_;   //!
   TBranch        *b_lumi_;   //!
   TBranch        *b_nVtx_;   //!
   TBranch        *b_selNOM_;   //!
   TBranch        *b_selVBF_;   //!
   TBranch        *b_selNOMsoft_;   //!
   TBranch        *b_selVBFsoft_;   //!
   TBranch        *b_nSoftJets_;   //!
   TBranch        *b_nSoftJets2_;   //!
   TBranch        *b_nSoftJets5_;   //!
   TBranch        *b_nJets_;   //!
   TBranch        *b_nBJets_;   //!
   TBranch        *b_b1_;   //!
   TBranch        *b_b2_;   //!
   TBranch        *b_q1_;   //!
   TBranch        *b_q2_;   //!
   TBranch        *b_mvaNOM_;   //!
   TBranch        *b_mvaVBF_;   //!
   TBranch        *b_softHt_;   //!
   TBranch        *b_rho_;   //!
   TBranch        *b_ht_;   //!
   TBranch        *b_met_;   //!
   TBranch        *b_metPhi_;   //!
   TBranch        *b_metSig_;   //!
   TBranch        *b_sphericity_;   //!
   TBranch        *b_aplanarity_;   //!
   TBranch        *b_dEtaMax_;   //!
   TBranch        *b_x1_;   //!
   TBranch        *b_x2_;   //!
   TBranch        *b_cosTheta_;   //!
   TBranch        *b_mqq_;   //!
   TBranch        *b_mbb_;   //!
   TBranch        *b_mbbReg_;   //!
   TBranch        *b_dEtaqq_;   //!
   TBranch        *b_dEtabb_;   //!
   TBranch        *b_ptbb_;   //!
   TBranch        *b_dPhiqq_;   //!
   TBranch        *b_dPhibb_;   //!
   TBranch        *b_yBoostqq_;   //!
   TBranch        *b_yStarqq_;   //!
   TBranch        *b_btagIdx;   //!
   TBranch        *b_blikNOMIdx;   //!
   TBranch        *b_blikVBFIdx;   //!
   TBranch        *b_etaIdx;   //!
   TBranch        *b_jetPuIdL;   //!
   TBranch        *b_jetPuIdM;   //!
   TBranch        *b_jetPuIdT;   //!
   TBranch        *b_jetIdL;   //!
   TBranch        *b_jetIdM;   //!
   TBranch        *b_jetIdT;   //!
   TBranch        *b_jetBtagL;   //!
   TBranch        *b_jetBtagM;   //!
   TBranch        *b_jetBtagT;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetBtag;   //!
   TBranch        *b_jetBlikNOM;   //!
   TBranch        *b_jetBlikVBF;   //!
   TBranch        *b_jetCSV;   //!
   TBranch        *b_jetPuMva;   //!
   TBranch        *b_jetJec;   //!
   TBranch        *b_jetRegPt;   //!
   TBranch        *b_jetRegE;   //!
   TBranch        *b_jetUnc;   //!
   TBranch        *b_jetQGL;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetMetPhi;   //!
   TBranch        *b_jetMass;   //!
   TBranch        *b_jetChf;   //!
   TBranch        *b_jetNhf;   //!
   TBranch        *b_jetPhf;   //!
   TBranch        *b_jetMuf;   //!
   TBranch        *b_jetElf;   //!
   TBranch        *b_triggerResult;   //!
   TBranch        *b_npu_;   //!
   TBranch        *b_wt;   //!

   ReadForHisto(TTree *tree=0);
   virtual ~ReadForHisto();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TFile *OUTF, TString SELECTION, bool isMC);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ReadForHisto_cxx
ReadForHisto::ReadForHisto(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("../../prod/VBFHbb/config/flatTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("../../prod/VBFHbb/config/flatTree.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("../../prod/VBFHbb/config/flatTree.root:/Hbb");
      dir->GetObject("events",tree);

   }
   Init(tree);
}

ReadForHisto::~ReadForHisto()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ReadForHisto::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ReadForHisto::LoadTree(Long64_t entry)
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

void ReadForHisto::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   btagIdx = 0;
   blikNOMIdx = 0;
   blikVBFIdx = 0;
   etaIdx = 0;
   jetPuIdL = 0;
   jetPuIdM = 0;
   jetPuIdT = 0;
   jetIdL = 0;
   jetIdM = 0;
   jetIdT = 0;
   jetBtagL = 0;
   jetBtagM = 0;
   jetBtagT = 0;
   jetPt = 0;
   jetBtag = 0;
   jetBlikNOM = 0;
   jetBlikVBF = 0;
   jetCSV = 0;
   jetPuMva = 0;
   jetJec = 0;
   jetRegPt = 0;
   jetRegE = 0;
   jetUnc = 0;
   jetQGL = 0;
   jetEta = 0;
   jetPhi = 0;
   jetMetPhi = 0;
   jetMass = 0;
   jetChf = 0;
   jetNhf = 0;
   jetPhf = 0;
   jetMuf = 0;
   jetElf = 0;
   triggerResult = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNo", &runNo, &b_run_);
   fChain->SetBranchAddress("evtNo", &evtNo, &b_evt_);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi_);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nVtx_);
   fChain->SetBranchAddress("selNOM", selNOM, &b_selNOM_);
   fChain->SetBranchAddress("selVBF", selVBF, &b_selVBF_);
   fChain->SetBranchAddress("selNOMsoft", selNOMsoft, &b_selNOMsoft_);
   fChain->SetBranchAddress("selVBFsoft", selVBFsoft, &b_selVBFsoft_);
   fChain->SetBranchAddress("nSoftJets", &nSoftJets, &b_nSoftJets_);
   fChain->SetBranchAddress("nSoftJets2", &nSoftJets2, &b_nSoftJets2_);
   fChain->SetBranchAddress("nSoftJets5", &nSoftJets5, &b_nSoftJets5_);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets_);
   fChain->SetBranchAddress("nBJets", &nBJets, &b_nBJets_);
   fChain->SetBranchAddress("b1", b1, &b_b1_);
   fChain->SetBranchAddress("b2", b2, &b_b2_);
   fChain->SetBranchAddress("q1", q1, &b_q1_);
   fChain->SetBranchAddress("q2", q2, &b_q2_);
   fChain->SetBranchAddress("mvaNOM", &mvaNOM, &b_mvaNOM_);
   fChain->SetBranchAddress("mvaVBF", &mvaVBF, &b_mvaVBF_);
   fChain->SetBranchAddress("softHt", &softHt, &b_softHt_);
   fChain->SetBranchAddress("rho", &rho, &b_rho_);
   fChain->SetBranchAddress("ht", &ht, &b_ht_);
   fChain->SetBranchAddress("met", &met, &b_met_);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi_);
   fChain->SetBranchAddress("metSig", &metSig, &b_metSig_);
   fChain->SetBranchAddress("sphericity", &sphericity, &b_sphericity_);
   fChain->SetBranchAddress("aplanarity", &aplanarity, &b_aplanarity_);
   fChain->SetBranchAddress("dEtaMax", &dEtaMax, &b_dEtaMax_);
   fChain->SetBranchAddress("x1", &x1, &b_x1_);
   fChain->SetBranchAddress("x2", &x2, &b_x2_);
   fChain->SetBranchAddress("cosTheta", cosTheta, &b_cosTheta_);
   fChain->SetBranchAddress("mqq", mqq, &b_mqq_);
   fChain->SetBranchAddress("mbb", mbb, &b_mbb_);
   fChain->SetBranchAddress("mbbReg", mbbReg, &b_mbbReg_);
   fChain->SetBranchAddress("dEtaqq", dEtaqq, &b_dEtaqq_);
   fChain->SetBranchAddress("dEtabb", dEtabb, &b_dEtabb_);
   fChain->SetBranchAddress("ptbb", ptbb, &b_ptbb_);
   fChain->SetBranchAddress("dPhiqq", dPhiqq, &b_dPhiqq_);
   fChain->SetBranchAddress("dPhibb", dPhibb, &b_dPhibb_);
   fChain->SetBranchAddress("yBoostqq", yBoostqq, &b_yBoostqq_);
   fChain->SetBranchAddress("yStarqq", yStarqq, &b_yStarqq_);
   fChain->SetBranchAddress("btagIdx", &btagIdx, &b_btagIdx);
   fChain->SetBranchAddress("blikNOMIdx", &blikNOMIdx, &b_blikNOMIdx);
   fChain->SetBranchAddress("blikVBFIdx", &blikVBFIdx, &b_blikVBFIdx);
   fChain->SetBranchAddress("etaIdx", &etaIdx, &b_etaIdx);
   fChain->SetBranchAddress("jetPuIdL", &jetPuIdL, &b_jetPuIdL);
   fChain->SetBranchAddress("jetPuIdM", &jetPuIdM, &b_jetPuIdM);
   fChain->SetBranchAddress("jetPuIdT", &jetPuIdT, &b_jetPuIdT);
   fChain->SetBranchAddress("jetIdL", &jetIdL, &b_jetIdL);
   fChain->SetBranchAddress("jetIdM", &jetIdM, &b_jetIdM);
   fChain->SetBranchAddress("jetIdT", &jetIdT, &b_jetIdT);
   fChain->SetBranchAddress("jetBtagL", &jetBtagL, &b_jetBtagL);
   fChain->SetBranchAddress("jetBtagM", &jetBtagM, &b_jetBtagM);
   fChain->SetBranchAddress("jetBtagT", &jetBtagT, &b_jetBtagT);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetBtag", &jetBtag, &b_jetBtag);
   fChain->SetBranchAddress("jetBlikNOM", &jetBlikNOM, &b_jetBlikNOM);
   fChain->SetBranchAddress("jetBlikVBF", &jetBlikVBF, &b_jetBlikVBF);
   fChain->SetBranchAddress("jetCSV", &jetCSV, &b_jetCSV);
   fChain->SetBranchAddress("jetPuMva", &jetPuMva, &b_jetPuMva);
   fChain->SetBranchAddress("jetJec", &jetJec, &b_jetJec);
   fChain->SetBranchAddress("jetRegPt", &jetRegPt, &b_jetRegPt);
   fChain->SetBranchAddress("jetRegE", &jetRegE, &b_jetRegE);
   fChain->SetBranchAddress("jetUnc", &jetUnc, &b_jetUnc);
   fChain->SetBranchAddress("jetQGL", &jetQGL, &b_jetQGL);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetMetPhi", &jetMetPhi, &b_jetMetPhi);
   fChain->SetBranchAddress("jetMass", &jetMass, &b_jetMass);
   fChain->SetBranchAddress("jetChf", &jetChf, &b_jetChf);
   fChain->SetBranchAddress("jetNhf", &jetNhf, &b_jetNhf);
   fChain->SetBranchAddress("jetPhf", &jetPhf, &b_jetPhf);
   fChain->SetBranchAddress("jetMuf", &jetMuf, &b_jetMuf);
   fChain->SetBranchAddress("jetElf", &jetElf, &b_jetElf);
   fChain->SetBranchAddress("triggerResult", &triggerResult, &b_triggerResult);
   fChain->SetBranchAddress("npu", &npu, &b_npu_);
   fChain->SetBranchAddress("puWt", &puWt, &b_wt);
   Notify();
}

Bool_t ReadForHisto::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ReadForHisto::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ReadForHisto::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ReadForHisto_cxx
