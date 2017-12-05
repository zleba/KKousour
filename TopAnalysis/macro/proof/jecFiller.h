//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Dec  2 00:13:48 2017 by ROOT version 6.06/08
// from TTree events/events
// found on file: /nfs/dust/cms/user/zlebcr/JEC/histos/jets1.root
//////////////////////////////////////////////////////////

#ifndef jecFiller_h
#define jecFiller_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector
#include <vector>
#include <TH1D.h>
#include <TH2.h>


class jecFiller : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> runNo = {fReader, "runNo"};
   TTreeReaderValue<Int_t> evtNo = {fReader, "evtNo"};
   TTreeReaderValue<Int_t> lumi = {fReader, "lumi"};
   TTreeReaderValue<Int_t> nvtx = {fReader, "nvtx"};
   TTreeReaderValue<Int_t> nJets = {fReader, "nJets"};
   TTreeReaderValue<Float_t> pvRho = {fReader, "pvRho"};
   TTreeReaderValue<Float_t> pvz = {fReader, "pvz"};
   TTreeReaderValue<Float_t> pvchi2 = {fReader, "pvchi2"};
   TTreeReaderValue<Float_t> pvndof = {fReader, "pvndof"};
   TTreeReaderValue<Float_t> rho = {fReader, "rho"};
   TTreeReaderValue<Float_t> ht = {fReader, "ht"};
   TTreeReaderValue<Float_t> metEt = {fReader, "metEt"};
   TTreeReaderValue<Float_t> metSumEt = {fReader, "metSumEt"};
   TTreeReaderValue<Float_t> metpt = {fReader, "metpt"};
   TTreeReaderValue<Float_t> metphi = {fReader, "metphi"};
   TTreeReaderValue<Float_t> metEtNoHF = {fReader, "metEtNoHF"};
   TTreeReaderValue<Float_t> metSumEtNoHF = {fReader, "metSumEtNoHF"};
   TTreeReaderValue<Float_t> metNoHFpt = {fReader, "metNoHFpt"};
   TTreeReaderValue<Float_t> metNoHFphi = {fReader, "metNoHFphi"};
   TTreeReaderValue<Float_t> metEtPuppi = {fReader, "metEtPuppi"};
   TTreeReaderValue<Float_t> metSumEtPuppi = {fReader, "metSumEtPuppi"};
   TTreeReaderValue<Float_t> metPuppipt = {fReader, "metPuppipt"};
   TTreeReaderValue<Float_t> metPuppiphi = {fReader, "metPuppiphi"};
   TTreeReaderArray<float> jetPt = {fReader, "jetPt"};
   TTreeReaderArray<float> jetUnc = {fReader, "jetUnc"};
   TTreeReaderArray<float> jetEta = {fReader, "jetEta"};
   TTreeReaderArray<float> jetPhi = {fReader, "jetPhi"};
   TTreeReaderArray<float> jetMass = {fReader, "jetMass"};
   TTreeReaderArray<float> jetEnergy = {fReader, "jetEnergy"};
   TTreeReaderArray<float> jetChf = {fReader, "jetChf"};
   TTreeReaderArray<float> jetNhf = {fReader, "jetNhf"};
   TTreeReaderArray<float> jetPhf = {fReader, "jetPhf"};
   TTreeReaderArray<float> jetMuf = {fReader, "jetMuf"};
   TTreeReaderArray<float> jetElf = {fReader, "jetElf"};
   TTreeReaderArray<int> jetChm = {fReader, "jetChm"};
   TTreeReaderArray<int> jetNhm = {fReader, "jetNhm"};
   TTreeReaderArray<int> jetPhm = {fReader, "jetPhm"};
   TTreeReaderArray<int> jetElm = {fReader, "jetElm"};
   TTreeReaderArray<int> jetMum = {fReader, "jetMum"};
   TTreeReaderValue<vector<bool>> jetIsBtag = {fReader, "jetIsBtag"};
   TTreeReaderArray<float> jetBtag = {fReader, "jetBtag"};
   TTreeReaderArray<float> HLTjetPt = {fReader, "HLTjetPt"};
   TTreeReaderArray<float> HLTjetEta = {fReader, "HLTjetEta"};
   TTreeReaderArray<float> HLTjetPhi = {fReader, "HLTjetPhi"};
   TTreeReaderArray<float> HLTjetMass = {fReader, "HLTjetMass"};
   TTreeReaderValue<vector<bool>> triggerBit = {fReader, "triggerBit"};
   TTreeReaderArray<int> triggerPre = {fReader, "triggerPre"};
   TTreeReaderValue<Int_t> nTriggerObjects = {fReader, "nTriggerObjects"};

   TH1D *hJetPt;
   TH2D *histoEmulated;
   TH2D *histoPtAll;

   jecFiller(TTree * /*tree*/ =0) { }
   virtual ~jecFiller() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(jecFiller,0);

};

#endif

#ifdef jecFiller_cxx
void jecFiller::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t jecFiller::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef jecFiller_cxx
