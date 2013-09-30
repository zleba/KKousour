#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "TFile.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "KKousour/CMGAnalysis/plugins/VbfHbbFlatTreeProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/BaseMET.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

using namespace std;
using namespace reco;

VbfHbbFlatTreeProducer::VbfHbbFlatTreeProducer(edm::ParameterSet const& cfg) 
{
  srcJets_            = cfg.getParameter<edm::InputTag>             ("jets");
  srcSoftJets_        = cfg.getParameter<edm::InputTag>             ("softJets");
  srcMET_             = cfg.getParameter<edm::InputTag>             ("met");
  srcRho_             = cfg.getParameter<edm::InputTag>             ("rho"); 
  srcBtag_            = cfg.getParameter<std::string>               ("btagger");
  puTag_              = cfg.getParameter<std::string>               ("putag");
  shiftJES_           = cfg.getParameter<double>                    ("shiftJES");
  ptMin_              = cfg.getParameter<double>                    ("ptMin");
  dEtaMin_            = cfg.getParameter<double>                    ("dEtaMin");
  btagThresholds_     = cfg.getParameter<vector<double> >           ("btagThresholds");
  saveJetProperties_  = cfg.getUntrackedParameter<bool>             ("saveJetProperties",false);
  saveSoftJets_       = cfg.getUntrackedParameter<bool>             ("saveSoftJets",false);
  savePartons_        = cfg.getUntrackedParameter<bool>             ("savePartons",false);
  forceNOM_           = cfg.getUntrackedParameter<bool>             ("forceNOM",false);
  forceVBF_           = cfg.getUntrackedParameter<bool>             ("forceVBF",false);
  forceQ50_           = cfg.getUntrackedParameter<bool>             ("forceQ50",false);
  
  srcPU_              = cfg.getUntrackedParameter<std::string>      ("pu","");
  srcGenJets_         = cfg.getUntrackedParameter<edm::InputTag>    ("genjets",edm::InputTag(""));
  srcGenParticles_    = cfg.getUntrackedParameter<edm::InputTag>    ("genparticles",edm::InputTag(""));

  triggerCache_       = triggerExpression::Data(cfg.getParameterSet("triggerConfiguration"));
  vtriggerAlias_      = cfg.getParameter<std::vector<std::string> > ("triggerAlias");
  vtriggerSelection_  = cfg.getParameter<std::vector<std::string> > ("triggerSelection");

  if (vtriggerAlias_.size() != vtriggerSelection_.size()) {
    cout<<"ERROR: the number of trigger aliases does not match the number of trigger names !!!"<<endl;
    return;
  }

  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    vtriggerSelector_.push_back(triggerExpression::parse(vtriggerSelection_[i]));
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::beginJob() 
{
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    triggerNamesHisto_->Fill(vtriggerSelection_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  triggerPassHisto_->SetBit(TH1::kCanRebin);
  //--- book the pileup histogram -----------
  puHisto_ = fs_->make<TH1F>("pileup","pileup",40,0,40);
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("selNOM"               ,&selNOM_            ,"selNOM_[3]/O");
  outTree_->Branch("selVBF"               ,&selVBF_            ,"selVBF_[3]/O");
  outTree_->Branch("selNOMsoft"           ,&selNOMsoft_        ,"selNOMsoft_[3]/O");
  outTree_->Branch("selVBFsoft"           ,&selVBFsoft_        ,"selVBFsoft_[3]/O");
  outTree_->Branch("nSoftJets"            ,&nSoftJets_         ,"nSoftJets_/I");
  outTree_->Branch("nSoftJets2"           ,&nSoftJets2_        ,"nSoftJets2_/I");
  outTree_->Branch("nSoftJets5"           ,&nSoftJets5_        ,"nSoftJets5_/I");
  outTree_->Branch("nJets"                ,&nJets_             ,"nJets_/I");
  outTree_->Branch("nBJets"               ,&nBJets_            ,"nBJets_/I");
  outTree_->Branch("b1"                   ,&b1_                ,"b1_[3]/I");
  outTree_->Branch("b2"                   ,&b2_                ,"b2_[3]/I");
  outTree_->Branch("q1"                   ,&q1_                ,"q1_[3]/I");
  outTree_->Branch("q2"                   ,&q2_                ,"q2_[3]/I");
  outTree_->Branch("mvaNOM"               ,&mvaNOM_            ,"mvaNOM_/F");
  outTree_->Branch("mvaVBF"               ,&mvaVBF_            ,"mvaVBF_/F");
  outTree_->Branch("softHt"               ,&softHt_            ,"softHt_/F");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("ht"                   ,&ht_                ,"ht_/F");
  outTree_->Branch("met"                  ,&met_               ,"met_/F");
  outTree_->Branch("metPhi"               ,&metPhi_            ,"metPhi_/F");
  outTree_->Branch("metSig"               ,&metSig_            ,"metSig_/F");
  outTree_->Branch("sphericity"           ,&sphericity_        ,"sphericity_/F");
  outTree_->Branch("aplanarity"           ,&aplanarity_        ,"aplanarity_/F");
  outTree_->Branch("dEtaMax"              ,&dEtaMax_           ,"dEtaMax_/F");
  outTree_->Branch("x1"                   ,&x1_                ,"x1_/F");
  outTree_->Branch("x2"                   ,&x2_                ,"x2_/F");
  outTree_->Branch("cosTheta"             ,&cosTheta_          ,"cosTheta_[3]/F");
  outTree_->Branch("mqq"                  ,&mqq_               ,"mqq_[3]/F");
  outTree_->Branch("mbb"                  ,&mbb_               ,"mbb_[3]/F");
  outTree_->Branch("mbbReg"               ,&mbbReg_            ,"mbbReg_[3]/F");
  outTree_->Branch("dEtaqq"               ,&dEtaqq_            ,"dEtaqq_[3]/F");
  outTree_->Branch("dEtabb"               ,&dEtabb_            ,"dEtabb_[3]/F");
  outTree_->Branch("ptbb"                 ,&ptbb_              ,"ptbb_[3]/F");
  outTree_->Branch("dPhiqq"               ,&dPhiqq_            ,"dPhiqq_[3]/F");
  outTree_->Branch("dPhibb"               ,&dPhibb_            ,"dPhibb_[3]/F");
  outTree_->Branch("yBoostqq"             ,&yBoostqq_          ,"yBoostqq_[3]/F");
  outTree_->Branch("yStarqq"              ,&yStarqq_           ,"yStarqq_[3]/F"); 
  //------------------------------------------------------------------
  btagIdx_        = new std::vector<int>;
  etaIdx_         = new std::vector<int>;
  blikNOMIdx_     = new std::vector<int>;
  blikVBFIdx_     = new std::vector<int>;
  btagL_          = new std::vector<bool>;
  btagM_          = new std::vector<bool>;
  btagT_          = new std::vector<bool>;
  puIdL_          = new std::vector<bool>;
  puIdM_          = new std::vector<bool>;
  puIdT_          = new std::vector<bool>;
  idL_            = new std::vector<bool>;
  idM_            = new std::vector<bool>;
  idT_            = new std::vector<bool>;
  pt_             = new std::vector<float>;
  btag_           = new std::vector<float>;
  blikNOM_        = new std::vector<float>;
  blikVBF_        = new std::vector<float>;
  csv_            = new std::vector<float>;
  puMva_          = new std::vector<float>;
  jec_            = new std::vector<float>;
  regPt_          = new std::vector<float>;
  regE_           = new std::vector<float>;
  unc_            = new std::vector<float>;
  qgl_            = new std::vector<float>;
  eta_            = new std::vector<float>;
  phi_            = new std::vector<float>;
  jetMetPhi_      = new std::vector<float>;
  mass_           = new std::vector<float>;
  energy_         = new std::vector<float>;
  chf_            = new std::vector<float>;
  nhf_            = new std::vector<float>;
  phf_            = new std::vector<float>;
  muf_            = new std::vector<float>;
  elf_            = new std::vector<float>;
  outTree_->Branch("btagIdx"              ,"vector<int>"       ,&btagIdx_);
  outTree_->Branch("blikNOMIdx"           ,"vector<int>"       ,&blikNOMIdx_);
  outTree_->Branch("blikVBFIdx"           ,"vector<int>"       ,&blikVBFIdx_);
  outTree_->Branch("etaIdx"               ,"vector<int>"       ,&etaIdx_);
  outTree_->Branch("jetPuIdL"             ,"vector<bool>"      ,&puIdL_);
  outTree_->Branch("jetPuIdM"             ,"vector<bool>"      ,&puIdM_);
  outTree_->Branch("jetPuIdT"             ,"vector<bool>"      ,&puIdT_);
  outTree_->Branch("jetIdL"               ,"vector<bool>"      ,&idL_);
  outTree_->Branch("jetIdM"               ,"vector<bool>"      ,&idM_);
  outTree_->Branch("jetIdT"               ,"vector<bool>"      ,&idT_); 
  outTree_->Branch("jetBtagL"             ,"vector<bool>"      ,&btagL_);
  outTree_->Branch("jetBtagM"             ,"vector<bool>"      ,&btagM_);
  outTree_->Branch("jetBtagT"             ,"vector<bool>"      ,&btagT_);
  outTree_->Branch("jetPt"                ,"vector<float>"     ,&pt_);
  outTree_->Branch("jetBtag"              ,"vector<float>"     ,&btag_); 
  outTree_->Branch("jetBlikNOM"           ,"vector<float>"     ,&blikNOM_); 
  outTree_->Branch("jetBlikVBF"           ,"vector<float>"     ,&blikVBF_); 
  outTree_->Branch("jetCSV"               ,"vector<float>"     ,&csv_); 
  outTree_->Branch("jetPuMva"             ,"vector<float>"     ,&puMva_);
  outTree_->Branch("jetJec"               ,"vector<float>"     ,&jec_);
  outTree_->Branch("jetRegPt"             ,"vector<float>"     ,&regPt_);
  outTree_->Branch("jetRegE"              ,"vector<float>"     ,&regE_);
  outTree_->Branch("jetUnc"               ,"vector<float>"     ,&unc_);
  outTree_->Branch("jetQGL"               ,"vector<float>"     ,&qgl_);
  outTree_->Branch("jetEta"               ,"vector<float>"     ,&eta_);
  outTree_->Branch("jetPhi"               ,"vector<float>"     ,&phi_);
  outTree_->Branch("jetMetPhi"            ,"vector<float>"     ,&jetMetPhi_);
  outTree_->Branch("jetMass"              ,"vector<float>"     ,&mass_);
  outTree_->Branch("jetEnergy"            ,"vector<float>"     ,&energy_);
  outTree_->Branch("jetChf"               ,"vector<float>"     ,&chf_);
  outTree_->Branch("jetNhf"               ,"vector<float>"     ,&nhf_);
  outTree_->Branch("jetPhf"               ,"vector<float>"     ,&phf_);
  outTree_->Branch("jetMuf"               ,"vector<float>"     ,&muf_);
  outTree_->Branch("jetElf"               ,"vector<float>"     ,&elf_);
  //------------------------------------------------------------------- 
  if (saveJetProperties_) {
    ptD_            = new std::vector<float>;
    ptD_QC_         = new std::vector<float>;
    axisMinor_      = new std::vector<float>;
    axisMajor_      = new std::vector<float>;
    axisMinor_QC_   = new std::vector<float>;
    axisMajor_QC_   = new std::vector<float>;
    pull_           = new std::vector<float>;
    pull_QC_        = new std::vector<float>;
    jetR_           = new std::vector<float>;
    jetRChg_QC_     = new std::vector<float>;
    vtx3dL_         = new std::vector<float>;
    vtx3deL_        = new std::vector<float>;
    vtxPt_          = new std::vector<float>;
    part_           = new std::vector<int>;
    nChg_QC_        = new std::vector<int>;
    nChg_ptCut_     = new std::vector<int>;
    nNeutral_ptCut_ = new std::vector<int>; 
    outTree_->Branch("jetPtD"               ,"vector<float>"     ,&ptD_);
    outTree_->Branch("jetPtD_QC"            ,"vector<float>"     ,&ptD_QC_);
    outTree_->Branch("jetAxisMinor"         ,"vector<float>"     ,&axisMinor_);
    outTree_->Branch("jetAxisMajor"         ,"vector<float>"     ,&axisMajor_);
    outTree_->Branch("jetAxisMinor_QC"      ,"vector<float>"     ,&axisMinor_QC_);
    outTree_->Branch("jetAxisMajor_QC"      ,"vector<float>"     ,&axisMajor_QC_);
    outTree_->Branch("jetPull"              ,"vector<float>"     ,&pull_);
    outTree_->Branch("jetPull_QC"           ,"vector<float>"     ,&pull_QC_);
    outTree_->Branch("jetR"                 ,"vector<float>"     ,&jetR_);
    outTree_->Branch("jetRChg_QC"           ,"vector<float>"     ,&jetRChg_QC_);
    outTree_->Branch("jetVtx3dL"            ,"vector<float>"     ,&vtx3dL_);
    outTree_->Branch("jetVtx3deL"           ,"vector<float>"     ,&vtx3deL_);
    outTree_->Branch("jetVtxPt"             ,"vector<float>"     ,&vtxPt_);
    outTree_->Branch("jetPart"              ,"vector<int>"       ,&part_); 
    outTree_->Branch("jetnChg_QC"           ,"vector<int>"       ,&nChg_QC_);
    outTree_->Branch("jetnChg_ptCut"        ,"vector<int>"       ,&nChg_ptCut_);
    outTree_->Branch("jetnNeutral_ptCut"    ,"vector<int>"       ,&nNeutral_ptCut_);
  }
  //------------------------------------------------------------------
  triggerResult_ = new std::vector<bool>;
  outTree_->Branch("triggerResult","vector<bool>",&triggerResult_);
  //------------------------------------------------------------------
  if (saveSoftJets_) { 
    softJetPt_  = new std::vector<float>;
    softJetEta_ = new std::vector<float>;
    softJetPhi_ = new std::vector<float>;
    softJetE_   = new std::vector<float>;
    outTree_->Branch("softJetPt" ,"vector<float>" ,&softJetPt_);
    outTree_->Branch("softJetEta","vector<float>" ,&softJetEta_);
    outTree_->Branch("softJetPhi","vector<float>" ,&softJetPhi_);
    outTree_->Branch("softJetE"  ,"vector<float>" ,&softJetE_);
  }
  //------------------- MC ---------------------------------
  outTree_->Branch("npu"           ,&npu_          ,"npu_/I");
  if (savePartons_) {
    partonId_       = new std::vector<int>;
    partonSt_       = new std::vector<int>;
    partonMatchIdx_ = new std::vector<int>;
    partonMatchDR_  = new std::vector<float>;
    partonPt_       = new std::vector<float>;
    partonEta_      = new std::vector<float>;
    partonPhi_      = new std::vector<float>;
    partonE_        = new std::vector<float>;
    outTree_->Branch("partonId"       ,"vector<int>"   ,&partonId_);
    outTree_->Branch("partonSt"       ,"vector<int>"   ,&partonSt_);
    outTree_->Branch("partonMatchIdx" ,"vector<int>"   ,&partonMatchIdx_);
    outTree_->Branch("partonMatchDR"  ,"vector<float>" ,&partonMatchDR_);
    outTree_->Branch("partonPt"       ,"vector<float>" ,&partonPt_);
    outTree_->Branch("partonEta"      ,"vector<float>" ,&partonEta_);
    outTree_->Branch("partonPhi"      ,"vector<float>" ,&partonPhi_);
    outTree_->Branch("partonE"        ,"vector<float>" ,&partonE_);
    genjetPt_       = new std::vector<float>;
    genjetEta_      = new std::vector<float>;
    genjetPhi_      = new std::vector<float>;
    genjetE_        = new std::vector<float>; 
    outTree_->Branch("genjetPt"       ,"vector<float>" ,&genjetPt_);
    outTree_->Branch("genjetEta"      ,"vector<float>" ,&genjetEta_);
    outTree_->Branch("genjetPhi"      ,"vector<float>" ,&genjetPhi_);
    outTree_->Branch("genjetE"        ,"vector<float>" ,&genjetE_);
  }
  qglCalc_    = new QGLCalculator();
  //jetReg_     = new JetRegressor("KKousour/CMGAnalysis/data/factoryJetRegNewGenJets_MLP.weights.xml");
  //jetReg_     = new JetRegressor("KKousour/CMGAnalysis/data/factory_JetReg_MLP.weights.xml");
  jetReg_     = new JetRegressor("KKousour/CMGAnalysis/data/factory_JetReg_BDT.weights.xml");
  bjetLikNOM_ = new BJetLikelihood("KKousour/CMGAnalysis/data/factory_BCanditate_NOM_BDT_GRAD.weights.xml");
  bjetLikVBF_ = new BJetLikelihood("KKousour/CMGAnalysis/data/factory_BCanditate_VBF_BDT_GRAD.weights.xml");
  discrNOM_   = new DiscriminatorMVA("KKousour/CMGAnalysis/data/factory_mvaNOM_BDT_GRAD.weights.xml");
  discrVBF_   = new DiscriminatorMVA("KKousour/CMGAnalysis/data/factory_mvaVBF_BDT_GRAD.weights.xml");
  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::endJob() 
{  
  if (puHisto_->GetEntries() > 0) {
    fillPuWeights();
  }
  delete qglCalc_;
  delete jetReg_;
  delete bjetLikNOM_;
  delete bjetLikVBF_;
  delete discrNOM_;
  delete discrVBF_;
  if (savePartons_) {
    delete partonSt_;
    delete partonId_;
    delete partonMatchIdx_;
    delete partonMatchDR_;
    delete partonPt_;
    delete partonEta_;
    delete partonPhi_;
    delete partonE_;
    delete genjetPt_;
    delete genjetEta_;
    delete genjetPhi_;
    delete genjetE_; 
  }
  if (saveSoftJets_) {
    delete softJetPt_;
    delete softJetEta_;
    delete softJetPhi_;
    delete softJetE_;
  }
  delete triggerResult_;
  delete btagIdx_;
  delete blikNOMIdx_;
  delete blikVBFIdx_; 
  delete etaIdx_;
  delete puIdL_;
  delete puIdM_;
  delete puIdT_;
  delete idL_;
  delete idM_;
  delete idT_;
  delete btagL_;
  delete btagM_;
  delete btagT_;
  delete pt_;
  delete btag_;
  delete blikNOM_;
  delete blikVBF_;
  delete csv_;
  delete puMva_;
  delete jec_;
  delete regPt_;
  delete regE_;
  delete unc_;
  delete qgl_;
  delete eta_;
  delete phi_;
  delete jetMetPhi_;
  delete mass_;
  delete energy_;
  delete chf_;
  delete nhf_;
  delete phf_;
  delete muf_;
  delete elf_;
  if (saveJetProperties_) {
    delete ptD_;
    delete ptD_QC_;
    delete axisMinor_;
    delete axisMajor_;
    delete axisMinor_QC_;
    delete axisMajor_QC_;
    delete pull_;
    delete pull_QC_;
    delete jetR_;
    delete jetRChg_QC_;
    delete vtx3dL_;
    delete vtx3deL_;
    delete vtxPt_;
    delete part_;
    delete nChg_QC_;
    delete nChg_ptCut_;
    delete nNeutral_ptCut_;
  }
  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  edm::Handle<edm::View<cmg::PFJet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<cmg::PFJet> cmg_jets = *jets;

  edm::Handle<edm::View<reco::Jet> > softJets;
  iEvent.getByLabel(srcSoftJets_,softJets);

  edm::Handle<edm::View<cmg::BaseMET> >  met;
  iEvent.getByLabel(srcMET_,met);

  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);

  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);

  edm::Handle<GenEventInfoProduct> hEventInfo;
  edm::Handle<GenParticleCollection> genParticles;
  
  edm::Handle<std::vector<cmg::PhysicsObjectWithPtr<edm::Ptr<reco::GenJet> > > > genjets;
  //edm::Handle<reco::GenJetCollection> genjets;
  //---------- pu -----------------------
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (!iEvent.isRealData()) {
    if (srcPU_ != "") {
      iEvent.getByLabel(srcPU_,PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator PUI;
      for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
        if (PUI->getBunchCrossing() == 0) {
          npu_ = PUI->getTrueNumInteractions();
        }
      }
    }
    puHisto_->Fill(npu_); 
  }
  //-------------- Trigger Info -----------------------------------
  triggerPassHisto_->Fill("totalEvents",1);
  if (triggerCache_.setEvent(iEvent,iSetup)) {
    for(unsigned itrig=0;itrig<vtriggerSelector_.size();itrig++) {
      bool result(false);
      if (vtriggerSelector_[itrig]) {
        if (triggerCache_.configurationUpdated()) {
          vtriggerSelector_[itrig]->init(triggerCache_);
        }
        result = (*(vtriggerSelector_[itrig]))(triggerCache_);
      }
      if (result) {
        triggerPassHisto_->Fill(vtriggerAlias_[itrig].c_str(),1);
      }
      triggerResult_->push_back(result);
    }
  }
     
  //----- at least one good vertex -----------
  bool cut_vtx = (recVtxs->size() > 0);
  
  if (cut_vtx) {
    //----- soft jets ----------------------
    softHt_     = 0;
    nSoftJets_  = 0; 
    nSoftJets2_ = 0; 
    nSoftJets5_ = 0; 
    for(int it=0;it<int(softJets->size());it++) {
      softHt_ += (*softJets)[it].pt();
      if (saveSoftJets_) {
        if (nSoftJets_ < 3) {
          softJetPt_ ->push_back((*softJets)[it].pt());
          softJetEta_->push_back((*softJets)[it].eta());
          softJetPhi_->push_back((*softJets)[it].phi());
          softJetE_  ->push_back((*softJets)[it].energy());
        }  
      }
      nSoftJets_++;
      if ((*softJets)[it].pt() > 2) {
        nSoftJets2_++;
        if ((*softJets)[it].pt() > 5) {
          nSoftJets5_++;
        }
      }
    }
    //----- PF jets ------------------------------
    nJets_ = 0;
    nBJets_ = 0;
    float ht(0.0);
    float sumP2(0.0),sumPxx(0.0),sumPxy(0.0),sumPxz(0.0),sumPyy(0.0),sumPyz(0.0),sumPzz(0.0);
    vector<float> vaJES;
    vector<TLorentzVector> vP4;
    for(edm::View<cmg::PFJet>::const_iterator ijet = cmg_jets.begin();ijet != cmg_jets.end(); ++ijet) {  
      float chf = ijet->component(1).fraction();
      float nhf = ijet->component(5).fraction();
      float phf = ijet->component(4).fraction();
      float muf = ijet->component(3).fraction();
      float elf = ijet->component(2).fraction();
      int chm = ijet->component(1).number();
      int npr = ijet->nConstituents();
      float eta = fabs(ijet->eta());
      float pt = ijet->pt();
      bool puIdL = ijet->passPuJetId(puTag_.c_str(),PileupJetIdentifier::kLoose);
      bool puIdM = ijet->passPuJetId(puTag_.c_str(),PileupJetIdentifier::kMedium);
      bool puIdT = ijet->passPuJetId(puTag_.c_str(),PileupJetIdentifier::kTight);
      bool idL = (npr>1 && phf<0.99 && nhf<0.99);
      bool idM = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
      bool idT = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.7 && muf<0.7 && chf>0 && chm>0) || eta>2.4));
      bool btagL = (ijet->bDiscriminator(srcBtag_.c_str()) > btagThresholds_[0]);
      bool btagM = (ijet->bDiscriminator(srcBtag_.c_str()) > btagThresholds_[1]);
      bool btagT = (ijet->bDiscriminator(srcBtag_.c_str()) > btagThresholds_[2]);
      if (idL && puIdL && (pt > ptMin_)) {
        vP4.push_back(TLorentzVector(ijet->px(),ijet->py(),ijet->pz(),ijet->energy()));
        sumP2  += ijet->p()  * ijet->p();
        sumPxx += ijet->px() * ijet->px();
        sumPxy += ijet->px() * ijet->py();
        sumPxz += ijet->px() * ijet->pz();
        sumPyy += ijet->py() * ijet->py();
        sumPyz += ijet->py() * ijet->pz();
        sumPzz += ijet->pz() * ijet->pz();   
        double aJES = 1.+shiftJES_*ijet->uncOnFourVectorScale();
        vaJES.push_back(aJES);
        chf_           ->push_back(chf);
        nhf_           ->push_back(nhf);
        phf_           ->push_back(phf);
        elf_           ->push_back(elf);
        muf_           ->push_back(muf);
        unc_           ->push_back(ijet->uncOnFourVectorScale());
        jec_           ->push_back(aJES/ijet->rawFactor());
        float dphiJetMet = fabs(deltaPhi(ijet->phi(),(*met)[0].phi()));
        regPt_         ->push_back((jetReg_->getTarget(*ijet,dphiJetMet,*rho,(*met)[0].et()))[0]/ijet->pt());
        regE_          ->push_back((jetReg_->getTarget(*ijet,dphiJetMet,*rho,(*met)[0].et()))[0]/ijet->pt());
        jetMetPhi_     ->push_back(dphiJetMet);
        pt_            ->push_back(aJES*pt);
        phi_           ->push_back(ijet->phi());
        eta_           ->push_back(ijet->eta());
        mass_          ->push_back(aJES*ijet->mass());
        energy_        ->push_back(aJES*ijet->energy());
        btag_          ->push_back(ijet->bDiscriminator(srcBtag_.c_str()));
        csv_           ->push_back(ijet->bDiscriminator("combinedSecondaryVertexBJetTags"));
        qgl_           ->push_back(qglCalc_->getQGL(*ijet,*rho));
        puMva_         ->push_back(ijet->puMva(puTag_.c_str()));
        puIdL_         ->push_back(puIdL);
        puIdM_         ->push_back(puIdM);
        puIdT_         ->push_back(puIdT);
        idL_           ->push_back(idL);
        idM_           ->push_back(idM);
        idT_           ->push_back(idT);
        btagL_         ->push_back(btagL);
        btagM_         ->push_back(btagM);
        btagT_         ->push_back(btagT);
        if (saveJetProperties_) {
          ptD_           ->push_back(ijet->ptd());
          ptD_QC_        ->push_back(ijet->ptdQC());
          vtx3dL_        ->push_back(ijet->vtx3dL());
          vtx3deL_       ->push_back(ijet->vtx3deL());
          vtxPt_         ->push_back(ijet->vtxPt());
          part_          ->push_back(npr);
          axisMajor_     ->push_back(ijet->axisMajor());
          axisMinor_     ->push_back(ijet->axisMinor());
          axisMajor_QC_  ->push_back(ijet->axisMajorQC());
          axisMinor_QC_  ->push_back(ijet->axisMinorQC());
          pull_          ->push_back(ijet->pull());
          pull_QC_       ->push_back(ijet->pullQC());
          jetR_          ->push_back(ijet->fmax());
          jetRChg_QC_    ->push_back(ijet->fmaxCharged());
          nChg_ptCut_    ->push_back(ijet->nChargedPtCut());
          nNeutral_ptCut_->push_back(ijet->nNeutralPtCut());
          nChg_QC_       ->push_back(ijet->nChargedQC());
        }
        ht += aJES*pt;
        nJets_++;
        if (btagM) {
          nBJets_++; 
        }
      }
    }// jet loop
    // ------- MC -----------------------------------------
    if (!iEvent.isRealData()) {
      //---------- genjets ------------------
      if (savePartons_) {
        iEvent.getByLabel(srcGenJets_, genjets);
        //for(reco::GenJetCollection::const_iterator igen = genjets->begin();igen != genjets->end(); ++igen) {
        for(std::vector<cmg::PhysicsObjectWithPtr<edm::Ptr<reco::GenJet> > >::const_iterator igen = genjets->begin();igen != genjets->end(); ++igen) {
          genjetPt_ ->push_back(igen->pt());
          genjetEta_->push_back(igen->eta());
          genjetPhi_->push_back(igen->phi());
          genjetE_  ->push_back(igen->energy()); 
        }// genjet loop
        //---------- partons ------------------
        iEvent.getByLabel(srcGenParticles_, genParticles);
        for(unsigned ip = 0; ip < genParticles->size(); ++ ip) {
          const GenParticle &p = (*genParticles)[ip];
          if (p.status() != 3) continue;
          int ndst3(0);
	  for(unsigned k = 0; k < p.numberOfDaughters(); k++) {
	    if (p.daughter(k)->status() == 3) {
              ndst3++; 
            }
	  }
          if (ndst3 > 0) continue;
          partonId_ ->push_back(p.pdgId());
	  partonSt_ ->push_back(p.status()); 
          partonPt_ ->push_back(p.pt());
          partonEta_->push_back(p.eta());
          partonPhi_->push_back(p.phi());
          partonE_  ->push_back(p.energy()); 
          //----- match partons with jets ------------
          float dRmin(1000);  
          int imatch(0);
          for(int j=0;j<nJets_;j++) {
            float dR = deltaR(p.eta(),p.phi(),(*eta_)[j],(*phi_)[j]);
            if (dR < dRmin) {
              imatch = j;
              dRmin = dR;
            }
          }   
          partonMatchIdx_->push_back(imatch);
          partonMatchDR_->push_back(dRmin);
        }// parton loop
      }
    }// if MC
    //---- compute sphericity -----------------------
    if (sumP2 > 0) {
      float Txx = sumPxx/sumP2;
      float Tyy = sumPyy/sumP2;
      float Tzz = sumPzz/sumP2;  
      float Txy = sumPxy/sumP2; 
      float Txz = sumPxz/sumP2; 
      float Tyz = sumPyz/sumP2;
      TMatrixDSym T(3);
      T(0,0) = Txx;
      T(0,1) = Txy;
      T(0,2) = Txz;
      T(1,0) = Txy;
      T(1,1) = Tyy;
      T(1,2) = Tyz;
      T(2,0) = Txz;
      T(2,1) = Tyz;
      T(2,2) = Tzz;
      TMatrixDSymEigen TEigen(T);
      TVectorD eigenValues(TEigen.GetEigenValues());
      sphericity_ = 1.5*(eigenValues(1)+eigenValues(2));
      aplanarity_ = 1.5*eigenValues(2); 
    }   
    ht_     = ht;
    rho_    = *rho;
    met_    = (*met)[0].et();
    metPhi_ = (*met)[0].phi();
    if ((*met)[0].sumEt() > 0) {
      metSig_ = (*met)[0].et()/(*met)[0].sumEt();
    }
    nVtx_   = recVtxs->size();
    run_    = iEvent.id().run();
    evt_    = iEvent.id().event();
    lumi_   = iEvent.id().luminosityBlock();
    if (nJets_ > 3) {
      //----- order jets according to the btag value -----
      order(*btag_,btagIdx_,4);
      //----- order jets according to the eta value -----
      order(*eta_,etaIdx_,4);
      //----- compute the blikelihood -------------------
      for(int j=0;j<4;j++) {
        blikNOM_->push_back(bjetLikNOM_->eval((*btagIdx_)[j],(*etaIdx_)[j],(*btag_)[j],(*eta_)[j]));
        blikVBF_->push_back(bjetLikVBF_->eval((*btagIdx_)[j],(*etaIdx_)[j],(*btag_)[j],(*eta_)[j])); 
      } 
      //----- order jets according to the blik value -----
      order(*blikNOM_,blikNOMIdx_,4);
      order(*blikVBF_,blikVBFIdx_,4);

      b1_[0] = (*btagIdx_)[0]; 
      b2_[0] = (*btagIdx_)[1];
      q1_[0] = (*btagIdx_)[2];
      q2_[0] = (*btagIdx_)[3];
      b1_[1] = (*blikNOMIdx_)[0]; 
      b2_[1] = (*blikNOMIdx_)[1];
      q1_[1] = (*blikNOMIdx_)[2];
      q2_[1] = (*blikNOMIdx_)[3]; 
      b1_[2] = (*blikVBFIdx_)[0]; 
      b2_[2] = (*blikVBFIdx_)[1];
      q1_[2] = (*blikVBFIdx_)[2];
      q2_[2] = (*blikVBFIdx_)[3];
      float sumE(0.0),sumPz(0.0);
      vector<TLorentzVector> vRegP4;
      for(int i=0;i<4;i++) {
        sumE  += vaJES[i]*vP4[i].Energy();
        sumPz += vaJES[i]*vP4[i].Pz();
        TLorentzVector p4(0,0,0,0);
        p4.SetPtEtaPhiE((*regPt_)[i]*(*pt_)[i],(*eta_)[i],(*phi_)[i],(*regE_)[i]*(*energy_)[i]);
        vRegP4.push_back(p4);
      }
      x1_ = (sumE+sumPz)/8000;
      x2_ = (sumE-sumPz)/8000;
      dEtaMax_    = fabs((*eta_)[(*etaIdx_)[0]]-(*eta_)[(*etaIdx_)[3]]);
      for(int k=0;k<3;k++) {
        dEtabb_[k]     = fabs((*eta_)[b1_[k]]-(*eta_)[b2_[k]]);
        dEtaqq_[k]     = fabs((*eta_)[q1_[k]]-(*eta_)[q2_[k]]);
        mbb_[k]        = (vaJES[b1_[k]]*vP4[b1_[k]]+vaJES[b2_[k]]*vP4[b2_[k]]).M();
        //mbbReg_[k]     = (vaJES[b1_[k]]*(*reg_)[b1_[k]]*vP4[b1_[k]]+vaJES[b2_[k]]*(*reg_)[b2_[k]]*vP4[b2_[k]]).M();
        mbbReg_[k]     = (vaJES[b1_[k]]*vRegP4[b1_[k]]+vaJES[b2_[k]]*vRegP4[b2_[k]]).M();   
        mqq_[k]        = (vaJES[q1_[k]]*vP4[q1_[k]]+vaJES[q2_[k]]*vP4[q2_[k]]).M();
        ptbb_[k]       = (vaJES[b1_[k]]*vP4[b1_[k]]+vaJES[b2_[k]]*vP4[b2_[k]]).Pt();
        dPhibb_[k]     = fabs(deltaPhi((*phi_)[b1_[k]],(*phi_)[b2_[k]]));
        dPhiqq_[k]     = fabs(deltaPhi((*phi_)[q1_[k]],(*phi_)[q2_[k]])); 
        yStarqq_[k]    = 0.5*(vP4[q1_[k]].Rapidity()-vP4[q2_[k]].Rapidity());
        yBoostqq_[k]   = 0.5*(vP4[q1_[k]].Rapidity()+vP4[q2_[k]].Rapidity());
        //----- plane angle at the qqbb rest frame ------
        TLorentzVector cmP4 = vP4[b1_[k]]+vP4[b2_[k]]+vP4[q1_[k]]+vP4[q2_[k]];
        TVector3 vb = cmP4.BoostVector();
        TLorentzVector cmP4b1 = vP4[b1_[k]];
        TLorentzVector cmP4b2 = vP4[b2_[k]];
        TLorentzVector cmP4q1 = vP4[q1_[k]];
        TLorentzVector cmP4q2 = vP4[q2_[k]];
        cmP4b1.Boost(-vb);
        cmP4b2.Boost(-vb);
        cmP4q1.Boost(-vb);
        cmP4q2.Boost(-vb);
        cosTheta_[k] = cos(((cmP4q1.Vect()).Cross(cmP4q2.Vect())).Angle((cmP4b1.Vect()).Cross(cmP4b2.Vect())));
      }
      mvaNOM_ = discrNOM_->eval(mqq_[1],dEtaqq_[1],dPhiqq_[1],(*btag_)[b1_[1]],(*btag_)[b2_[1]],(*qgl_)[b1_[1]],(*qgl_)[b2_[1]],(*qgl_)[q1_[1]],(*qgl_)[q2_[1]],softHt_,nSoftJets2_,cosTheta_[1],x1_,x2_,ht_);
      mvaVBF_ = discrVBF_->eval(mqq_[2],dEtaqq_[2],dPhiqq_[2],(*btag_)[b1_[2]],(*btag_)[b2_[2]],(*qgl_)[b1_[2]],(*qgl_)[b2_[2]],(*qgl_)[q1_[2]],(*qgl_)[q2_[2]],softHt_,nSoftJets2_,cosTheta_[2],x1_,x2_,ht_);

      bool cut_Btag = ((*btagM_)[b1_[0]] && (*btagL_)[b2_[0]]);
      if ((*pt_)[0]>70 && (*pt_)[1]>55 && (*pt_)[2]>40 && (*pt_)[3]>20) {
        for (int i=0;i<3;i++) { 
          if (dEtaqq_[i]>2.0 && mqq_[i]>200) {
            selNOMsoft_[i] = true;
          }
        } 
      } 
      if ((*pt_)[0]>80 && (*pt_)[1]>65 && (*pt_)[2]>50 && (*pt_)[3]>35) {
        for (int i=0;i<3;i++) { 
          if (dEtaqq_[i]>2.5 && mqq_[i]>250) {
            selNOM_[i] = true;
          }
        }
      }
      if (ht_>150 && cut_Btag) {
        for (int i=0;i<3;i++) { 
          if (dEtaqq_[i]>3.0 && mqq_[i]>500) {
            selVBFsoft_[i] = true;
          }
        }
      }
      if (ht_>350 && cut_Btag) {
        for (int i=0;i<3;i++) { 
          if (dEtaqq_[i]>3.5 && mqq_[i]>800) {
            selVBF_[i] = true;
          }
        }
      }
      //--------------------------------------------------------------
      bool SELECTION(false);
      if (forceNOM_) {
        SELECTION += (selNOMsoft_[0] || selNOMsoft_[1] || selNOMsoft_[2]);
      }
      if (forceVBF_) {
        SELECTION += (selVBFsoft_[0] || selVBFsoft_[1] || selVBFsoft_[2]);
      }
      if (!forceNOM_ && !forceVBF_) {
        SELECTION = (dEtaqq_[0]>dEtaMin_ || dEtaqq_[1]>dEtaMin_ || dEtaqq_[2]>dEtaMin_);
      }
      if (SELECTION) { 
        outTree_->Fill();     
      }
    }// if nJets > 3
  }// if vtx
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::order(vector<float> const& v, vector<int>* idx, int Nmax)
{
  vector<float> tmp;
  int N = min(int(v.size()),Nmax);
  for(int i=0;i<N;i++) {
    idx->push_back(i);
    tmp.push_back(v[i]); 
  }
  for(int i=0;i<N;i++) {
    for(int j=i+1;j<N;j++) {
      if (tmp[j] > tmp[i]) {
        int k = (*idx)[i];
        (*idx)[i] = (*idx)[j];
        (*idx)[j] = k;
        float x = tmp[i]; 
        tmp[i] = tmp[j];
        tmp[j] = x;
      }
    } 
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::initialize()
{
  for(int i=0;i<3;i++) {
    selNOM_[i]     = false;
    selVBF_[i]     = false;
    selNOMsoft_[i] = false;
    selVBFsoft_[i] = false;
    cosTheta_[i]   = -999;
    mbb_[i]        = -999;
    mbbReg_[i]     = -999;
    dEtaqq_[i]     = -999;
    dEtabb_[i]     = -999;
    yBoostqq_[i]   = -999;
    yStarqq_[i]    = -999;
    dPhiqq_[i]     = -999;
    dPhibb_[i]     = -999;
    ptbb_[i]       = -999;
    b1_[i]         = -999;
    b2_[i]         = -999;
    q1_[i]         = -999;
    q2_[i]         = -999;
  } 
  run_            = -999;
  evt_            = -999;
  lumi_           = -999;
  nVtx_           = -999;
  nSoftJets_      = -999;
  nSoftJets2_     = -999;
  nSoftJets5_     = -999;
  nJets_          = -999;
  nBJets_         = -999;
  rho_            = -999;
  met_            = -999;
  metPhi_         = -999;
  metSig_         = -999;
  ht_             = -999;
  softHt_         = -999;
  sphericity_     = -999;
  aplanarity_     = -999;
  dEtaMax_        = -999;
  x1_             = -999;
  x2_             = -999;
  mvaNOM_         = -999;
  mvaVBF_         = -999;
  pt_             ->clear();
  eta_            ->clear();
  phi_            ->clear();
  jetMetPhi_      ->clear();
  mass_           ->clear();
  energy_         ->clear();
  chf_            ->clear();
  nhf_            ->clear();
  phf_            ->clear();
  elf_            ->clear();
  muf_            ->clear();
  jec_            ->clear();
  regPt_          ->clear();
  regE_           ->clear();
  csv_            ->clear();
  btag_           ->clear();
  btagIdx_        ->clear();
  blikNOM_        ->clear();
  blikVBF_        ->clear();
  blikNOMIdx_     ->clear();
  blikVBFIdx_     ->clear();
  etaIdx_         ->clear();
  qgl_            ->clear();
  puMva_          ->clear();
  unc_            ->clear();
  puIdL_          ->clear();
  puIdM_          ->clear();
  puIdT_          ->clear();
  idL_            ->clear();
  idM_            ->clear();
  idT_            ->clear();
  btagL_          ->clear();
  btagM_          ->clear();
  btagT_          ->clear();
  triggerResult_  ->clear();
  if (saveJetProperties_) {
    ptD_            ->clear();
    ptD_QC_         ->clear();
    vtx3dL_         ->clear();
    vtx3deL_        ->clear();
    vtxPt_          ->clear();
    part_           ->clear();
    axisMinor_      ->clear();
    axisMajor_      ->clear();
    axisMinor_QC_   ->clear();
    axisMajor_QC_   ->clear();
    pull_           ->clear(); 
    pull_QC_        ->clear();
    jetR_           ->clear();
    jetRChg_QC_     ->clear();
    nChg_QC_        ->clear();
    nChg_ptCut_     ->clear();
    nNeutral_ptCut_ ->clear();
  }
  if (saveSoftJets_) {
    softJetPt_      ->clear(); 
    softJetEta_     ->clear();
    softJetPhi_     ->clear();
    softJetE_       ->clear();
  }
  //----- MC -------
  npu_ = -999;
  if (savePartons_) {
    partonSt_      ->clear();
    partonId_      ->clear();
    partonMatchIdx_->clear();
    partonMatchDR_ ->clear();
    partonPt_      ->clear();
    partonEta_     ->clear();
    partonPhi_     ->clear();
    partonE_       ->clear();
    genjetPt_      ->clear();
    genjetEta_     ->clear();
    genjetPhi_     ->clear();
    genjetE_       ->clear();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void VbfHbbFlatTreeProducer::fillPuWeights()
{
  edm::FileInPath f1("KKousour/CMGAnalysis/data/pu2012DCONLY.root");
  TFile *inf = new TFile(TString(f1.fullPath()));
  TH1D *puRef = (TH1D*)inf->Get("pileup");
  puRef->Sumw2();
  puRef->Rebin(25);
  puRef->Scale(1./puRef->Integral());
  puHisto_->Scale(1./puHisto_->Integral());
  TH1D *puWeightHisto = (TH1D*)puRef->Clone("puWeight");
  puWeightHisto->Divide(puHisto_);
  float wt(1.0);
  int npu;
  TBranch *b = outTree_->Branch("puWt",&wt,"wt/F");
  outTree_->SetBranchAddress("npu",&npu);
  cout<<"Filling tree with pu weight"<<endl;
  for(int i=0;i<outTree_->GetEntries();i++) {
    outTree_->GetEvent(i);
    wt = puWeightHisto->GetBinContent(puWeightHisto->FindBin(npu));
    b->Fill();
  }
  inf->Close();
  delete inf;
  delete puRef; 
}
//////////////////////////////////////////////////////////////////////////////////////////
VbfHbbFlatTreeProducer::~VbfHbbFlatTreeProducer() 
{
}

DEFINE_FWK_MODULE(VbfHbbFlatTreeProducer);
















