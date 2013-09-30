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
#include "TMath.h"

#include "KKousour/CMGAnalysis/plugins/ZJetsFlatTreeProducer.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

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

using namespace std;
using namespace reco;

ZJetsFlatTreeProducer::ZJetsFlatTreeProducer(edm::ParameterSet const& cfg) 
{
  srcJets_            = cfg.getParameter<edm::InputTag> ("jets");
  srcBtag_            = cfg.getParameter<std::string>   ("btagger");
  srcMET_             = cfg.getParameter<edm::InputTag> ("met");
  srcElectrons_       = cfg.getParameter<edm::InputTag> ("electrons");
  srcMuons_           = cfg.getParameter<edm::InputTag> ("muons");
  srcRho_             = cfg.getParameter<edm::InputTag> ("rho"); 
  minJetPt_           = cfg.getParameter<double>        ("minJetPt");
  maxJetEta_          = cfg.getParameter<double>        ("maxJetEta");
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
void ZJetsFlatTreeProducer::beginJob() 
{
  jetP4_        = new std::vector<LorentzVector>();
  jetPt_        = new std::vector<float>(); 
  jetPtD_       = new std::vector<float>(); 
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
  jetVtxPt_     = new std::vector<float>();
  jetVtx3dL_    = new std::vector<float>();
  jetVtx3deL_   = new std::vector<float>();
  jetId_        = new std::vector<int>();
  elP4_         = new std::vector<LorentzVector>();
  elPt_         = new std::vector<float>();
  elEta_        = new std::vector<float>();
  elPhi_        = new std::vector<float>();
  elE_          = new std::vector<float>(); 
  elIso_        = new std::vector<float>(); 
  elMva_        = new std::vector<float>();
  elId_         = new std::vector<int>();
  elCh_         = new std::vector<int>();
  muP4_         = new std::vector<LorentzVector>(); 
  muPt_         = new std::vector<float>();
  muEta_        = new std::vector<float>();
  muPhi_        = new std::vector<float>();
  muE_          = new std::vector<float>();
  muIso_        = new std::vector<float>(); 
  muId_         = new std::vector<int>();
  muCh_         = new std::vector<int>();
  triggerResult_ = new std::vector<bool>;
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
  outTree_->Branch("nmuons"      ,&nmuons_      ,"nmuons_/I");
  outTree_->Branch("nelectrons"  ,&nelectrons_  ,"nelectrons_/I");
  outTree_->Branch("rho"         ,&rho_         ,"rho_/F");
  outTree_->Branch("met"         ,&met_         ,"met_/F");
  outTree_->Branch("metPhi"      ,&metPhi_      ,"metPhi_/F");
  outTree_->Branch("triggerResult","vector<bool>",&triggerResult_);
  //---------- jets -----------------------------------------
  outTree_->Branch("jetPt"    ,"vector<float>"     ,&jetPt_);
  outTree_->Branch("jetPtD"   ,"vector<float>"     ,&jetPtD_);
  outTree_->Branch("jetEta"   ,"vector<float>"     ,&jetEta_);
  outTree_->Branch("jetPhi"   ,"vector<float>"     ,&jetPhi_);
  outTree_->Branch("jetE"     ,"vector<float>"     ,&jetE_);
  outTree_->Branch("jetBtag"  ,"vector<float>"     ,&jetBtag_);
  outTree_->Branch("jetChf"   ,"vector<float>"     ,&jetChf_);
  outTree_->Branch("jetNhf"   ,"vector<float>"     ,&jetNhf_);
  outTree_->Branch("jetPhf"   ,"vector<float>"     ,&jetPhf_);
  outTree_->Branch("jetMuf"   ,"vector<float>"     ,&jetMuf_);
  outTree_->Branch("jetElf"   ,"vector<float>"     ,&jetElf_);
  outTree_->Branch("jetPuMva" ,"vector<float>"     ,&jetPuMva_);
  outTree_->Branch("jetVtxPt" ,"vector<float>"     ,&jetVtxPt_);
  outTree_->Branch("jetVtx3dL" ,"vector<float>"     ,&jetVtx3dL_);
  outTree_->Branch("jetVtx3deL" ,"vector<float>"     ,&jetVtx3deL_);
  outTree_->Branch("jetId"    ,"vector<int>"       ,&jetId_);
  //---------- electrons -----------------------------------------
  outTree_->Branch("elPt"     ,"vector<float>"     ,&elPt_);
  outTree_->Branch("elEta"    ,"vector<float>"     ,&elEta_);
  outTree_->Branch("elPhi"    ,"vector<float>"     ,&elPhi_);
  outTree_->Branch("elE"      ,"vector<float>"     ,&elE_);
  outTree_->Branch("elIso"    ,"vector<float>"     ,&elIso_);
  outTree_->Branch("elMva"    ,"vector<float>"     ,&elMva_);
  outTree_->Branch("elId"     ,"vector<int>"       ,&elId_);
  outTree_->Branch("elCh"     ,"vector<int>"       ,&elCh_);
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
  outTree_->Branch("eeMPF"    ,&eeMPF_             ,"eeMPF_/F");
  //---------- di-muons -----------------------------------------
  outTree_->Branch("mmPt"     ,&mmPt_              ,"mmPt_/F");
  outTree_->Branch("mmEta"    ,&mmEta_             ,"mmEta_/F");
  outTree_->Branch("mmPhi"    ,&mmPhi_             ,"mmPhi_/F");
  outTree_->Branch("mmE"      ,&mmE_               ,"mmE_/F");
  outTree_->Branch("mmM"      ,&mmM_               ,"mmM_/F");
  outTree_->Branch("mmMPF"    ,&mmMPF_             ,"mmMPF_/F");
  //---------- di-electrons+jet -----------------------------------------
  outTree_->Branch("eejPt"    ,&eejPt_             ,"eejPt_/F");
  outTree_->Branch("eejY"     ,&eejY_              ,"eejY_/F");
  outTree_->Branch("eejM"     ,&eejM_              ,"eejM_/F");
  //---------- di-muons+jet -----------------------------------------
  outTree_->Branch("mmjPt"    ,&mmjPt_             ,"mmjPt_/F");
  outTree_->Branch("mmjY"     ,&mmjY_              ,"mmjY_/F");
  outTree_->Branch("mmjM"     ,&mmjM_              ,"mmjM_/F");

  cout<<"Tree booked ......... "<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void ZJetsFlatTreeProducer::endJob() 
{
  delete elP4_;
  delete elPt_;
  delete elEta_;
  delete elPhi_;
  delete elE_;
  delete elIso_;
  delete elMva_;
  delete elId_;
  delete elCh_;
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
  delete jetPtD_;
  delete jetPuMva_; 
  delete jetVtxPt_;
  delete jetVtx3dL_;
  delete jetVtx3deL_;
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
  delete triggerResult_;
  for(unsigned i=0;i<vtriggerSelector_.size();i++) {
    delete vtriggerSelector_[i];
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void ZJetsFlatTreeProducer::initialize()
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
  nelectrons_   = -999;
  nmuons_       = -999; 
  eePt_         = -999;
  eeEta_        = -999;
  eePhi_        = -999;
  eeE_          = -999;
  eeM_          = -999;
  eeMPF_        = -999;
  mmPt_         = -999;
  mmEta_        = -999;
  mmPhi_        = -999;
  mmE_          = -999;
  mmM_          = -999;
  mmMPF_        = -999;
  eejPt_        = -999;
  eejY_         = -999;
  eejM_         = -999;
  mmjPt_        = -999;
  mmjY_         = -999;
  mmjM_         = -999; 
  elP4_->clear();
  elPt_->clear();
  elEta_->clear();
  elPhi_->clear();
  elE_->clear();
  elIso_->clear();
  elMva_->clear();
  elId_->clear();
  elCh_->clear();
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
  jetPtD_->clear();
  jetPuMva_->clear(); 
  jetVtxPt_->clear();
  jetVtx3dL_->clear();
  jetVtx3deL_->clear();
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
  triggerResult_->clear();
}
//////////////////////////////////////////////////////////////////////////////////////////

void ZJetsFlatTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
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
  
  //---- muons block --------------------------------------------  
  nmuons_ = 0;
  edm::Handle<edm::View<cmg::Muon> > muons;
  iEvent.getByLabel(srcMuons_,muons);
  edm::View<cmg::Muon> cmg_muons = *muons;
  for(edm::View<cmg::Muon>::const_iterator imu = cmg_muons.begin();imu != cmg_muons.end(); ++imu) { 
    if (imu->pt() < minMuonPt_) continue;
    muPt_ ->push_back(imu->pt());
    muEta_->push_back(imu->eta());
    muPhi_->push_back(imu->phi());	
    muE_  ->push_back(imu->energy());	
    muP4_ ->push_back(imu->p4());
    int id(0);
    if (
      imu->isGlobalMuon() && 
      fabs(imu->normalizedChi2()) < 10 &&
      imu->numberOfValidMuonHits() > 0 &&
      imu->numberOfMatches() >1 && 
      fabs(imu->dxy()) < 0.2 &&
      imu->numberOfValidPixelHits() > 0 && 
      imu->numberOfValidTrackerHits() > 8
    ) {
      id = 1;
    }
    muId_->push_back(id);
    muIso_->push_back(imu->relIso());
    muCh_->push_back(imu->charge());      
    nmuons_++;
  }// muon loop 
  //---- electrons block --------------------------------------------  
  nelectrons_ = 0;
  edm::Handle<edm::View<cmg::Electron> > electrons;
  iEvent.getByLabel(srcElectrons_,electrons);
  for(edm::View<cmg::Electron>::const_iterator iel = electrons->begin();iel != electrons->end(); ++iel) { 
    if (iel->pt() < minElectronPt_) continue;
    elPt_ ->push_back(iel->pt());
    elEta_->push_back(iel->eta());
    elPhi_->push_back(iel->phi());	
    elE_  ->push_back(iel->energy());
    elP4_ ->push_back(iel->p4());
    int id(0); 
    float sigmaIetaIeta                  = iel->sigmaIetaIeta();
    float hadronicOverEm                 = iel->hadronicOverEm();
    float deltaPhiSuperClusterTrackAtVtx = iel->deltaPhiSuperClusterTrackAtVtx();
    float deltaEtaSuperClusterTrackAtVtx = iel->deltaEtaSuperClusterTrackAtVtx();
    float etaSC = iel->sourcePtr()->get()->superCluster()->eta();
    if (etaSC < 1.4442) {
      if (sigmaIetaIeta < 0.01 && deltaPhiSuperClusterTrackAtVtx < 0.06 && deltaEtaSuperClusterTrackAtVtx < 0.004 && hadronicOverEm < 0.12) 
      id = 1;
    }// if EB
    if (etaSC > 1.5660) {
      if (sigmaIetaIeta < 0.03 && deltaPhiSuperClusterTrackAtVtx < 0.03 && deltaEtaSuperClusterTrackAtVtx < 0.007 && hadronicOverEm < 0.10) 
      id = 1;
    }// if EE	
    elMva_->push_back(iel->mva());
    elId_->push_back(id);
    elIso_->push_back(iel->relIso());
    elCh_->push_back(iel->charge());
    nelectrons_++;	      
  }// electron loop
  //---- jets block --------------------------------------------  
  njets_ = 0;
  edm::Handle<edm::View<cmg::PFJet> > jets;
  iEvent.getByLabel(srcJets_,jets);
  edm::View<cmg::PFJet> cmg_jets = *jets;
  for(edm::View<cmg::PFJet>::const_iterator ijet = cmg_jets.begin();ijet != cmg_jets.end(); ++ijet) { 
    //---- remove jets beyond acceptance --------
    if (ijet->pt() < minJetPt_ || fabs(ijet->eta()) > maxJetEta_) continue;
    //---- cross clean with electrons and muons --
    bool matched(false);
    for(int iel=0;iel<nelectrons_;iel++) {
      double dR = deltaR(ijet->p4(),elP4_->at(iel));
      if (dR < 0.25) {
        matched = true;
        break;
      }
    }
    for(int imu=0;imu<nmuons_;imu++) {
      double dR = deltaR(ijet->p4(),muP4_->at(imu));
      if (dR < 0.25) {
        matched = true;
        break;
      }
    }
    if (matched) continue;

    jetP4_->push_back(ijet->p4());
    jetPt_->push_back(ijet->pt());
    jetPtD_->push_back(ijet->ptd());
    jetEta_->push_back(ijet->eta());
    jetPhi_->push_back(ijet->phi());	
    jetE_->push_back(ijet->energy());	
    jetBtag_->push_back(ijet->bDiscriminator(srcBtag_.c_str()));
    jetPuMva_->push_back(ijet->puMva("full"));
    jetVtx3dL_->push_back(ijet->vtx3dL());
    jetVtx3deL_->push_back(ijet->vtx3deL());
    jetVtxPt_->push_back(ijet->vtxPt());
    float chf = ijet->component(1).fraction();
    float chm = ijet->component(1).number();
    float nhf = ijet->component(5).fraction();
    float phf = ijet->component(4).fraction();
    float muf = ijet->component(3).fraction();
    float elf = ijet->component(2).fraction();
    int npr = ijet->nConstituents();
    int id(0);
    if (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(ijet->eta())<=2.4 && elf<0.99 && chf>0 && chm>0) || fabs(ijet->eta())>2.4)) {
      id = 1;
    }
    if (npr>1 && phf<0.9 && nhf<0.9 && ((fabs(ijet->eta())<=2.4 && elf<0.9 && muf<0.9 && chf>0 && chm>0) || fabs(ijet->eta())>2.4)) {
      id = 2;
    }
    jetId_->push_back(id);
    jetChf_->push_back(chf);
    jetNhf_->push_back(nhf);
    jetPhf_->push_back(phf);
    jetMuf_->push_back(muf);
    jetElf_->push_back(elf);	      
    njets_++;
  }// jet loop
  run_    = iEvent.id().run();
  evt_    = iEvent.id().event();
  lumi_   = iEvent.id().luminosityBlock();
  if (nVtx_ > 0 && (nelectrons_ > 1 || nmuons_ > 1) && njets_ > 0) {
    if (nelectrons_ > 1) {
      LorentzVector eeP4 = elP4_->at(0) + elP4_->at(1); 
      eePt_  = eeP4.Pt();
      eeEta_ = eeP4.Eta();
      eePhi_ = eeP4.Phi();
      eeE_   = eeP4.E();
      eeM_   = eeP4.M();
      eeMPF_ = 1.0+met_*cos(deltaPhi(metPhi_,eePhi_))/eePt_;
      LorentzVector eejP4 = eeP4 + jetP4_->at(0); 
      eejPt_ = eejP4.Pt();
      eejY_  = eejP4.Rapidity();
      eejM_  = eejP4.M();
    }
    if (nmuons_ > 1) {
      LorentzVector mmP4 = muP4_->at(0) + muP4_->at(1); 
      mmPt_  = mmP4.Pt();
      mmEta_ = mmP4.Eta();
      mmPhi_ = mmP4.Phi();
      mmE_   = mmP4.E();
      mmM_   = mmP4.M();
      mmMPF_ = 1.0+met_*cos(deltaPhi(metPhi_,mmPhi_))/mmPt_;
      LorentzVector mmjP4 = mmP4 + jetP4_->at(0); 
      mmjPt_ = mmjP4.Pt();
      mmjY_  = mmjP4.Rapidity();
      mmjM_  = mmjP4.M();
    }
    outTree_->Fill();
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
ZJetsFlatTreeProducer::~ZJetsFlatTreeProducer() 
{
}

DEFINE_FWK_MODULE(ZJetsFlatTreeProducer);
