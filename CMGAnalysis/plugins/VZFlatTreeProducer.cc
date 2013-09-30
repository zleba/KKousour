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
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TLorentzVector.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "KKousour/CMGAnalysis/plugins/VZFlatTreeProducer.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "AnalysisDataFormats/CMGTools/interface/PFJet.h"
#include "AnalysisDataFormats/CMGTools/interface/Muon.h"
#include "AnalysisDataFormats/CMGTools/interface/Electron.h"
#include "AnalysisDataFormats/CMGTools/interface/BaseMET.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

using namespace std;
using namespace reco;

VZFlatTreeProducer::VZFlatTreeProducer(edm::ParameterSet const& cfg) 
{
  srcJets_            = cfg.getParameter<edm::InputTag> ("jets");
  srcBtag_            = cfg.getParameter<std::string>   ("btagger");
  srcMET_             = cfg.getParameter<edm::InputTag> ("met");
  srcElectrons_       = cfg.getParameter<edm::InputTag> ("electrons");
  srcMuons_           = cfg.getParameter<edm::InputTag> ("muons");
  srcRho_             = cfg.getParameter<edm::InputTag> ("rho"); 
  minJetPt_           = cfg.getParameter<double>        ("minJetPt");
  maxJetEta_          = cfg.getParameter<double>        ("maxJetEta");
  minJetPtSel_        = cfg.getParameter<double>        ("minJetPtSel");
  maxJetEtaSel_       = cfg.getParameter<double>        ("maxJetEtaSel");
  minBJetPtSel_       = cfg.getParameter<double>        ("minBJetPtSel");
  maxBJetEtaSel_      = cfg.getParameter<double>        ("maxBJetEtaSel");
  minMuonPt_          = cfg.getParameter<double>        ("minMuonPt");     
  minElectronPt_      = cfg.getParameter<double>        ("minElectronPt"); 
  srcPU_              = cfg.getUntrackedParameter<std::string>      ("pu","");
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
void VZFlatTreeProducer::beginJob() 
{
  jetP4_        = new std::vector<LorentzVector>();
  jetPt_        = new std::vector<float>(); 
  jetQGL_       = new std::vector<float>(); 
  jetEta_       = new std::vector<float>();
  jetPhi_       = new std::vector<float>();
  jetE_         = new std::vector<float>();
  jetBtag_      = new std::vector<float>();
  jetChf_       = new std::vector<float>();
  jetNhf_       = new std::vector<float>();
  jetPhf_       = new std::vector<float>();
  jetMuf_       = new std::vector<float>();
  jetElf_       = new std::vector<float>(); 
  jetPuMva_     = new std::vector<float>();
  jetId_        = new std::vector<int>();
  jetPuIdL_     = new std::vector<bool>;
  jetPuIdM_     = new std::vector<bool>;
  jetPuIdT_     = new std::vector<bool>;
  elP4_         = new std::vector<LorentzVector>();
  elPt_         = new std::vector<float>();
  elEta_        = new std::vector<float>();
  elPhi_        = new std::vector<float>();
  elE_          = new std::vector<float>(); 
  elIso_        = new std::vector<float>(); 
  elMva_        = new std::vector<float>();
  elCh_         = new std::vector<int>();
  elEB_         = new std::vector<bool>();
  elEE_         = new std::vector<bool>();
  elIdL_        = new std::vector<bool>();
  elIdM_        = new std::vector<bool>();
  elIdT_        = new std::vector<bool>();
  muP4_         = new std::vector<LorentzVector>(); 
  muPt_         = new std::vector<float>();
  muEta_        = new std::vector<float>();
  muPhi_        = new std::vector<float>();
  muE_          = new std::vector<float>();
  muIso_        = new std::vector<float>(); 
  muId_         = new std::vector<int>();
  muCh_         = new std::vector<int>();
  triggerResult_= new std::vector<bool>;
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetBit(TH1::kCanRebin);
  for(unsigned i=0;i<vtriggerSelection_.size();i++) {
    triggerNamesHisto_->Fill(vtriggerSelection_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  triggerPassHisto_->SetBit(TH1::kCanRebin);
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"       ,&run_         ,"run_/I");
  outTree_->Branch("evtNo"       ,&evt_         ,"evt_/I");
  outTree_->Branch("lumi"        ,&lumi_        ,"lumi_/I");
  outTree_->Branch("nvtx"        ,&nVtx_        ,"nVtx_/I");
  outTree_->Branch("njets"       ,&njets_       ,"njets_/I");
  outTree_->Branch("njetsSel"    ,&njetsSel_    ,"njetsSel_/I");
  outTree_->Branch("nbjetsSel"   ,&nbjetsSel_   ,"nbjetsSel_/I");
  outTree_->Branch("nmuons"      ,&nmuons_      ,"nmuons_/I");
  outTree_->Branch("nelectrons"  ,&nelectrons_  ,"nelectrons_/I");
  outTree_->Branch("rho"         ,&rho_         ,"rho_/F");
  outTree_->Branch("met"         ,&met_         ,"met_/F");
  outTree_->Branch("metPhi"      ,&metPhi_      ,"metPhi_/F");
  outTree_->Branch("sphericity"  ,&sphericity_  ,"sphericity_/F");
  outTree_->Branch("aplanarity"  ,&aplanarity_  ,"aplanarity_/F");
  outTree_->Branch("triggerResult","vector<bool>",&triggerResult_);
  //---------- jets -----------------------------------------
  outTree_->Branch("jetPt"    ,"vector<float>"     ,&jetPt_);
  outTree_->Branch("jetEta"   ,"vector<float>"     ,&jetEta_);
  outTree_->Branch("jetPhi"   ,"vector<float>"     ,&jetPhi_);
  outTree_->Branch("jetE"     ,"vector<float>"     ,&jetE_);
  outTree_->Branch("jetBtag"  ,"vector<float>"     ,&jetBtag_);
  outTree_->Branch("jetChf"   ,"vector<float>"     ,&jetChf_);
  outTree_->Branch("jetNhf"   ,"vector<float>"     ,&jetNhf_);
  outTree_->Branch("jetPhf"   ,"vector<float>"     ,&jetPhf_);
  outTree_->Branch("jetMuf"   ,"vector<float>"     ,&jetMuf_);
  outTree_->Branch("jetElf"   ,"vector<float>"     ,&jetElf_);
  outTree_->Branch("jetQGL"   ,"vector<float>"     ,&jetQGL_);
  outTree_->Branch("jetPuMva" ,"vector<float>"     ,&jetPuMva_);
  outTree_->Branch("jetId"    ,"vector<int>"       ,&jetId_);
  outTree_->Branch("jetPuIdL" ,"vector<bool>"      ,&jetPuIdL_);
  outTree_->Branch("jetPuIdM" ,"vector<bool>"      ,&jetPuIdM_);
  outTree_->Branch("jetPuIdT" ,"vector<bool>"      ,&jetPuIdT_); 
  //---------- electrons -----------------------------------------
  outTree_->Branch("elPt"     ,"vector<float>"     ,&elPt_);
  outTree_->Branch("elEta"    ,"vector<float>"     ,&elEta_);
  outTree_->Branch("elPhi"    ,"vector<float>"     ,&elPhi_);
  outTree_->Branch("elE"      ,"vector<float>"     ,&elE_);
  outTree_->Branch("elIso"    ,"vector<float>"     ,&elIso_);
  outTree_->Branch("elMva"    ,"vector<float>"     ,&elMva_);
  outTree_->Branch("elCh"     ,"vector<int>"       ,&elCh_);
  outTree_->Branch("elIdL"    ,"vector<bool>"      ,&elIdL_);
  outTree_->Branch("elIdM"    ,"vector<bool>"      ,&elIdM_);
  outTree_->Branch("elIdT"    ,"vector<bool>"      ,&elIdT_);
  outTree_->Branch("elEB"     ,"vector<bool>"      ,&elEB_);
  outTree_->Branch("elEE"     ,"vector<bool>"      ,&elEE_);
  //---------- muons -----------------------------------------
  outTree_->Branch("muPt"     ,"vector<float>"     ,&muPt_);
  outTree_->Branch("muEta"    ,"vector<float>"     ,&muEta_);
  outTree_->Branch("muPhi"    ,"vector<float>"     ,&muPhi_);
  outTree_->Branch("muE"      ,"vector<float>"     ,&muE_);
  outTree_->Branch("muIso"    ,"vector<float>"     ,&muIso_);
  outTree_->Branch("muId"     ,"vector<int>"       ,&muId_); 
  outTree_->Branch("muCh"     ,"vector<int>"       ,&muCh_);
  //---------- di-electrons -----------------------------------------
  outTree_->Branch("eePt"     ,&eePt_              ,"eePt_/F");
  outTree_->Branch("eeEta"    ,&eeEta_             ,"eeEta_/F");
  outTree_->Branch("eePhi"    ,&eePhi_             ,"eePhi_/F");
  outTree_->Branch("eeE"      ,&eeE_               ,"eeE_/F");
  outTree_->Branch("eeM"      ,&eeM_               ,"eeM_/F");
  //---------- di-muons -----------------------------------------
  outTree_->Branch("mmPt"     ,&mmPt_              ,"mmPt_/F");
  outTree_->Branch("mmEta"    ,&mmEta_             ,"mmEta_/F");
  outTree_->Branch("mmPhi"    ,&mmPhi_             ,"mmPhi_/F");
  outTree_->Branch("mmE"      ,&mmE_               ,"mmE_/F");
  outTree_->Branch("mmM"      ,&mmM_               ,"mmM_/F");
  //---------- electron-muon -----------------------------------------
  outTree_->Branch("emPt"     ,&emPt_              ,"emPt_/F");
  outTree_->Branch("emEta"    ,&emEta_             ,"emEta_/F");
  outTree_->Branch("emPhi"    ,&emPhi_             ,"emPhi_/F");
  outTree_->Branch("emE"      ,&emE_               ,"emE_/F");
  outTree_->Branch("emM"      ,&emM_               ,"emM_/F");
  //---------- di-jets -----------------------------------------
  outTree_->Branch("jjPt"     ,&jjPt_              ,"jjPt_/F");
  outTree_->Branch("jjEta"    ,&jjEta_             ,"jjEta_/F");
  outTree_->Branch("jjPhi"    ,&jjPhi_             ,"jjPhi_/F");
  outTree_->Branch("jjE"      ,&jjE_               ,"jjE_/F");
  outTree_->Branch("jjM"      ,&jjM_               ,"jjM_/F");
  outTree_->Branch("jjDPhi"   ,&jjDPhi_            ,"jjDPhi_/F");
  outTree_->Branch("jjDEta"   ,&jjDEta_            ,"jjDEta_/F");
  //---------- di-electrons+di-jet -----------------------------------------
  outTree_->Branch("eejjPt"   ,&eejjPt_            ,"eejjPt_/F");
  outTree_->Branch("eejjY"    ,&eejjY_             ,"eejjY_/F");
  outTree_->Branch("eejjM"    ,&eejjM_             ,"eejjM_/F");
  //---------- di-muons+di-jet -----------------------------------------
  outTree_->Branch("mmjjPt"   ,&mmjjPt_            ,"mmjjPt_/F");
  outTree_->Branch("mmjjY"    ,&mmjjY_             ,"mmjjY_/F");
  outTree_->Branch("mmjjM"    ,&mmjjM_             ,"mmjjM_/F");

  cout<<"Tree booked ......... "<<endl;
  QGL_ = new QGLCalculator();
}
//////////////////////////////////////////////////////////////////////////////////////////
void VZFlatTreeProducer::endJob() 
{
  delete QGL_;
  delete elP4_;
  delete elPt_;
  delete elEta_;
  delete elPhi_;
  delete elE_;
  delete elIso_;
  delete elMva_;
  delete elCh_;
  delete elIdL_;
  delete elIdM_;
  delete elIdT_;
  delete elEB_;
  delete elEE_;
  delete muP4_;
  delete muPt_;
  delete muEta_;
  delete muPhi_;
  delete muE_;
  delete muIso_;
  delete muId_;
  delete muCh_;
  delete jetP4_;
  delete jetPt_;
  delete jetQGL_;
  delete jetPuMva_; 
  delete jetEta_;
  delete jetPhi_;
  delete jetE_;
  delete jetBtag_; 
  delete jetChf_;
  delete jetNhf_;
  delete jetPhf_;
  delete jetMuf_;
  delete jetElf_;
  delete jetId_;
  delete jetPuIdL_;
  delete jetPuIdM_;
  delete jetPuIdT_;
  delete triggerResult_;
  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void VZFlatTreeProducer::initialize()
{
  run_          = -999;
  evt_          = -999;
  lumi_         = -999;
  nVtx_         = -999;
  rho_          = -999;
  npu_          = -999;
  met_          = -999;
  metPhi_       = -999;
  njets_        = -999;
  njetsSel_     = -999;
  nbjetsSel_    = -999;
  nelectrons_   = -999;
  nmuons_       = -999; 
  eePt_         = -999;
  eeEta_        = -999;
  eePhi_        = -999;
  eeE_          = -999;
  eeM_          = -999;
  mmPt_         = -999;
  mmEta_        = -999;
  mmPhi_        = -999;
  mmE_          = -999;
  mmM_          = -999;
  emPt_         = -999;
  emEta_        = -999;
  emPhi_        = -999;
  emE_          = -999;
  emM_          = -999;
  jjPt_         = -999;
  jjEta_        = -999;
  jjPhi_        = -999;
  jjE_          = -999;
  jjM_          = -999;
  jjDPhi_       = -999;
  jjDEta_       = -999;
  eejjPt_       = -999;
  eejjY_        = -999;
  eejjM_        = -999;
  mmjjPt_       = -999;
  mmjjY_        = -999;
  mmjjM_        = -999; 
  sphericity_   = -999;
  aplanarity_   = -999;
  elP4_->clear();
  elPt_->clear();
  elEta_->clear();
  elPhi_->clear();
  elE_->clear();
  elIso_->clear();
  elMva_->clear();
  elCh_->clear();
  elIdL_->clear();
  elIdM_->clear();
  elIdT_->clear();
  elEB_->clear();
  elEE_->clear();
  muP4_->clear();
  muPt_->clear();
  muEta_->clear();
  muPhi_->clear();
  muE_->clear();
  muIso_->clear();
  muId_->clear();
  muCh_->clear();
  jetP4_->clear();
  jetPt_->clear();
  jetQGL_->clear();
  jetPuMva_->clear(); 
  jetEta_->clear();
  jetPhi_->clear();
  jetE_->clear();
  jetBtag_->clear(); 
  jetChf_->clear();
  jetNhf_->clear();
  jetPhf_->clear();
  jetMuf_->clear();
  jetElf_->clear();
  jetId_->clear();
  jetPuIdL_->clear();
  jetPuIdM_->clear();
  jetPuIdT_->clear();
  triggerResult_->clear();
}
//////////////////////////////////////////////////////////////////////////////////////////

void VZFlatTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();
  
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
  
  // ------- MC -----------------------------------------
  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  if (!iEvent.isRealData()) {
    //---------- pu -----------------------
    if (srcPU_ != "") {
      iEvent.getByLabel(srcPU_,PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator PUI;
      for(PUI = PupInfo->begin(); PUI != PupInfo->end(); ++PUI) {
        if (PUI->getBunchCrossing() == 0) {
          npu_ = PUI->getTrueNumInteractions();
        }
      }
    }
  }
  //---- rho ---------------------------------------------------------------
  edm::Handle<double> rho;
  iEvent.getByLabel(srcRho_,rho);
  rho_ = *rho;
  //---- reco vertices block --------------------------------------------
  edm::Handle<reco::VertexCollection> recVtxs;
  iEvent.getByLabel("offlinePrimaryVertices",recVtxs);
  nVtx_ = 0;
  for(VertexCollection::const_iterator i_vtx = recVtxs->begin(); i_vtx != recVtxs->end(); ++i_vtx) {  
    if (!i_vtx->isFake() && (fabs(i_vtx->z()) < 24) && (i_vtx->ndof() >= 4)) {
      nVtx_++;
    }
  }
  //---- met block ---------------------------------------------
  edm::Handle<edm::View<cmg::BaseMET> >  met;
  iEvent.getByLabel(srcMET_,met);
  met_ = (*met)[0].et();
  metPhi_ = (*met)[0].phi();
  float sumP2(0.0),sumPxx(0.0),sumPxy(0.0),sumPxz(0.0),sumPyy(0.0),sumPyz(0.0),sumPzz(0.0);
  //---- muons block --------------------------------------------  
  nmuons_ = 0;
  edm::Handle<edm::View<cmg::Muon> > muons;
  iEvent.getByLabel(srcMuons_,muons);
  edm::View<cmg::Muon> cmg_muons = *muons;
  for(edm::View<cmg::Muon>::const_iterator imu = cmg_muons.begin();imu != cmg_muons.end(); ++imu) { 
    if (imu->pt() < minMuonPt_) continue;
    int id(0);
    if (
      imu->isGlobalMuon() && 
      imu->isPF() &&
      fabs(imu->normalizedChi2()) < 10 &&
      imu->numberOfValidMuonHits() > 0 &&
      imu->numberOfMatches() > 1 && 
      fabs(imu->dxy()) < 0.2 &&
      fabs(imu->dz()) < 0.5 &&
      imu->numberOfValidPixelHits() > 0 && 
      imu->numberOfValidTrackerHits() > 5
    ) {
      id = 1;
    }
    muPt_ ->push_back(imu->pt());
    muEta_->push_back(imu->eta());
    muPhi_->push_back(imu->phi());	
    muE_  ->push_back(imu->energy());	
    muP4_ ->push_back(imu->p4());
    muId_->push_back(id);
    muIso_->push_back(imu->relIso(0.5,0,0.4));
    muCh_->push_back(imu->charge());      
    nmuons_++;
    if (id > 0) {
      sumP2  += imu->p() * imu->p();
      sumPxx += imu->px() * imu->px();
      sumPxy += imu->px() * imu->py();
      sumPxz += imu->px() * imu->pz();
      sumPyy += imu->py() * imu->py();
      sumPyz += imu->py() * imu->pz();
      sumPzz += imu->pz() * imu->pz();
    }
  }// muon loop 
  //---- electrons block --------------------------------------------  
  nelectrons_ = 0;
  edm::Handle<edm::View<cmg::Electron> > electrons;
  iEvent.getByLabel(srcElectrons_,electrons);
  for(edm::View<cmg::Electron>::const_iterator iel = electrons->begin();iel != electrons->end(); ++iel) { 
    if (iel->pt() < minElectronPt_) continue;
    float ecalE                          = iel->sourcePtr()->get()->ecalEnergy();
    float trkP                           = ecalE/iel->sourcePtr()->get()->eSuperClusterOverP(); 
    float dxy                            = iel->dxy();
    float dz                             = iel->dz();
    float fep                            = fabs(1./ecalE-1./trkP);                              
    float sigmaIetaIeta                  = iel->sigmaIetaIeta();
    float hadronicOverEm                 = iel->hadronicOverEm();
    float deltaPhiSuperClusterTrackAtVtx = iel->deltaPhiSuperClusterTrackAtVtx();
    float deltaEtaSuperClusterTrackAtVtx = iel->deltaEtaSuperClusterTrackAtVtx();
    float etaSC = fabs(iel->sourcePtr()->get()->superCluster()->eta());
    bool idL(false),idM(false),idT(false),isEB(false),isEE(false);
    if (etaSC < 1.442) {
      isEB = true;
      if (hadronicOverEm < 0.12 && sigmaIetaIeta < 0.01 && dxy < 0.02 && fep < 0.05) {
        if (deltaPhiSuperClusterTrackAtVtx < 0.15 && deltaEtaSuperClusterTrackAtVtx < 0.007 && dz < 0.2)  {
          idL = true;
        }
        if (deltaPhiSuperClusterTrackAtVtx < 0.06 && deltaEtaSuperClusterTrackAtVtx < 0.004 && dz < 0.1)  {
          idM = true;
        }
        if (deltaPhiSuperClusterTrackAtVtx < 0.03 && deltaEtaSuperClusterTrackAtVtx < 0.004 && dz < 0.1)  {
          idT = true;
        }
      }
    }// if EB
    if (etaSC > 1.566) {
      isEE = true; 
      if (hadronicOverEm < 0.1 && sigmaIetaIeta < 0.03 && dxy < 0.02 && fep < 0.05) {
        if (deltaPhiSuperClusterTrackAtVtx < 0.10 && deltaEtaSuperClusterTrackAtVtx < 0.009 && dz < 0.2)  {
          idL = true;
        }
        if (deltaPhiSuperClusterTrackAtVtx < 0.03 && deltaEtaSuperClusterTrackAtVtx < 0.007 && dz < 0.1)  {
          idM = true;
        }
        if (deltaPhiSuperClusterTrackAtVtx < 0.02 && deltaEtaSuperClusterTrackAtVtx < 0.005 && dz < 0.1)  {
          idT = true;
        }
      }
    }// if EE	
    if (!idL) continue;
    elEB_->push_back(isEB); 
    elEE_->push_back(isEE);
    elPt_ ->push_back(iel->pt());
    elEta_->push_back(iel->eta());
    elPhi_->push_back(iel->phi());	
    elE_  ->push_back(iel->energy());
    elP4_ ->push_back(iel->p4());
    elMva_->push_back(iel->mva());
    elIdL_->push_back(idL);
    elIdM_->push_back(idM);
    elIdT_->push_back(idT);
    float Aeff;
    if (etaSC < 1.0) Aeff = 0.13;
    else if (etaSC >= 1.0 && etaSC < 1.479) Aeff = 0.14;
    else if (etaSC >= 1.479 && etaSC < 2.0) Aeff = 0.07;
    else if (etaSC >= 2.0 && etaSC < 2.2) Aeff = 0.09;
    else if (etaSC >= 2.2 && etaSC < 2.3) Aeff = 0.11;
    else if (etaSC >= 2.3 && etaSC < 2.5) Aeff = 0.11;
    else Aeff = 0.14;
    float iso = iel->chargedHadronIso(0.3)+max((iel->neutralHadronIso(0.3)+iel->photonIso(0.3)-(*rho)*Aeff),0.0);
    //elIso_->push_back(iel->relIso());
    elIso_->push_back(iso/iel->pt());
    elCh_->push_back(iel->charge());
    nelectrons_++;
    sumP2  += iel->p()  * iel->p();
    sumPxx += iel->px() * iel->px();
    sumPxy += iel->px() * iel->py();
    sumPxz += iel->px() * iel->pz();
    sumPyy += iel->py() * iel->py();
    sumPyz += iel->py() * iel->pz();
    sumPzz += iel->pz() * iel->pz();	      
  }// electron loop
  //---- jets block --------------------------------------------  
  njets_     = 0;
  njetsSel_  = 0;
  nbjetsSel_ = 0;
  edm::Handle<edm::View<cmg::PFJet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<cmg::PFJet> cmg_jets = *jets;
  for(edm::View<cmg::PFJet>::const_iterator ijet = cmg_jets.begin();ijet != cmg_jets.end(); ++ijet) { 
    bool puIdL = ijet->passPuJetId("full",PileupJetIdentifier::kLoose);
    bool puIdM = ijet->passPuJetId("full",PileupJetIdentifier::kMedium);
    bool puIdT = ijet->passPuJetId("full",PileupJetIdentifier::kTight);
    //---- remove jets beyond acceptance --------
    if (ijet->pt() < minJetPt_ || fabs(ijet->eta()) > maxJetEta_ || !puIdL) continue;
    //---- cross clean with isolated electrons and muons --
    bool matched(false);
    for(int iel=0;iel<nelectrons_;iel++) {
      double dR = deltaR(ijet->p4(),(*elP4_)[iel]);
      if (dR < 0.5 && (*elIso_)[iel] < 0.15) {
        matched = true;
        break;
      }
    }
    for(int imu=0;imu<nmuons_;imu++) {
      double dR = deltaR(ijet->p4(),(*muP4_)[imu]);
      if (dR < 0.5 && (*muIso_)[imu] < 0.15) {
        matched = true;
        break;
      }
    }
    if (matched) continue;

    jetP4_->push_back(ijet->p4());
    jetPt_->push_back(ijet->pt());
    jetEta_->push_back(ijet->eta());
    jetPhi_->push_back(ijet->phi());	
    jetE_->push_back(ijet->energy());	
    jetBtag_->push_back(ijet->bDiscriminator(srcBtag_.c_str()));
    jetPuMva_->push_back(ijet->puMva("full"));
    float chf = ijet->component(1).fraction();
    float chm = ijet->component(1).number();
    float nhf = ijet->component(5).fraction();
    float phf = ijet->component(4).fraction();
    float muf = ijet->component(3).fraction();
    float elf = ijet->component(2).fraction();
    float eta = fabs(ijet->eta());
    int npr = ijet->nConstituents();
    int id(0);
    if (npr>1 && phf<0.99 && nhf<0.99 && ((eta<=2.4 && elf<0.99 && chf>0 && chm>0) || eta>2.4)) {
      id = 1;
    }
    if (npr>1 && phf<0.9 && nhf<0.9 && ((eta<=2.4 && elf<0.9 && muf<0.9 && chf>0 && chm>0) || eta>2.4)) {
      id = 2;
    }
    jetId_->push_back(id);
    jetPuIdL_->push_back(puIdL);
    jetPuIdM_->push_back(puIdM);
    jetPuIdT_->push_back(puIdT);
    jetChf_->push_back(chf);
    jetNhf_->push_back(nhf);
    jetPhf_->push_back(phf);
    jetMuf_->push_back(muf);
    jetElf_->push_back(elf);
    jetQGL_->push_back(QGL_->getQGL(*ijet,rho_));          	      
    njets_++;
    if (ijet->pt() > minJetPtSel_ && fabs(ijet->eta()) < maxJetEtaSel_ && (id > 0) && puIdL) {
      njetsSel_++; 
    }
    if (ijet->pt() > minBJetPtSel_ && fabs(ijet->eta()) < maxBJetEtaSel_ && (id > 0) && puIdL && ijet->bDiscriminator(srcBtag_.c_str()) > 0.679) {
        nbjetsSel_++; 
    }
    if (id > 0) {
      sumP2  += ijet->p()  * ijet->p();
      sumPxx += ijet->px() * ijet->px();
      sumPxy += ijet->px() * ijet->py();
      sumPxz += ijet->px() * ijet->pz();
      sumPyy += ijet->py() * ijet->py();
      sumPyz += ijet->py() * ijet->pz();
      sumPzz += ijet->pz() * ijet->pz();
    }
  }// jet loop
  //---------- compute sphericity --------------------
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
  //--------------------------------------------------
  run_    = iEvent.id().run();
  evt_    = iEvent.id().event();
  lumi_   = iEvent.id().luminosityBlock();
  if (nVtx_ > 0 && (nelectrons_ > 0 || nmuons_ > 0) && njets_ > 1) {
    LorentzVector jjP4(0,0,0,0);
    if (njets_ > 1) {
      jjP4 = (*jetP4_)[0] + (*jetP4_)[1]; 
      jjPt_   = jjP4.Pt();
      jjEta_  = jjP4.Eta();
      jjPhi_  = jjP4.Phi();
      jjE_    = jjP4.E();
      jjM_    = jjP4.M(); 
      jjDPhi_ = deltaPhi((*jetPhi_)[0],(*jetPhi_)[1]);
      jjDEta_ = fabs((*jetEta_)[0]-(*jetEta_)[1]);
    }
    if (nelectrons_ > 0 && nmuons_ > 0) {
      LorentzVector emP4 = (*elP4_)[0] + (*muP4_)[0]; 
      emPt_  = emP4.Pt();
      emEta_ = emP4.Eta();
      emPhi_ = emP4.Phi();
      emE_   = emP4.E();
      emM_   = emP4.M();
    }
    if (nelectrons_ > 1) {
      LorentzVector eeP4 = (*elP4_)[0] + (*elP4_)[1]; 
      eePt_  = eeP4.Pt();
      eeEta_ = eeP4.Eta();
      eePhi_ = eeP4.Phi();
      eeE_   = eeP4.E();
      eeM_   = eeP4.M();
      LorentzVector eejjP4 = eeP4 + jjP4; 
      eejjPt_ = eejjP4.Pt();
      eejjY_  = eejjP4.Rapidity();
      eejjM_  = eejjP4.M();
    }
    if (nmuons_ > 1) {
      LorentzVector mmP4 = (*muP4_)[0] + (*muP4_)[1]; 
      mmPt_  = mmP4.Pt();
      mmEta_ = mmP4.Eta();
      mmPhi_ = mmP4.Phi();
      mmE_   = mmP4.E();
      mmM_   = mmP4.M();
      LorentzVector mmjjP4 = mmP4 + jjP4; 
      mmjjPt_ = mmjjP4.Pt();
      mmjjY_  = mmjjP4.Rapidity();
      mmjjM_  = mmjjP4.M();
    }
    outTree_->Fill();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
VZFlatTreeProducer::~VZFlatTreeProducer() 
{
}

DEFINE_FWK_MODULE(VZFlatTreeProducer);
