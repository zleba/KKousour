#include <string>
#include <cmath>
#include <functional>
#include <vector>
#include <cassert>
#include "Math/SpecFuncMathMore.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

#include "KKousour/TopAnalysis/plugins/BoostedTTbarFlatTreeProducer.h"

using namespace std;
using namespace reco;
using namespace fastjet;

BoostedTTbarFlatTreeProducer::BoostedTTbarFlatTreeProducer(edm::ParameterSet const& cfg)
{ 
  jetsToken             = consumes<pat::JetCollection>(cfg.getParameter<edm::InputTag>("jets"));
  genjetsToken          = consumes<GenJetCollection>(cfg.getUntrackedParameter<edm::InputTag>("genjets",edm::InputTag("")));
  //  muonsToken            = consumes<pat::MuonCollection>(cfg.getParameter<edm::InputTag>("muons"));
  //electronsToken        = consumes<pat::ElectronCollection>(cfg.getParameter<edm::InputTag>("electrons"));
  //metToken              = consumes<pat::METCollection>(cfg.getParameter<edm::InputTag>("met"));
  candsToken            = consumes<pat::PackedCandidateCollection>(cfg.getParameter<edm::InputTag>("candidates"));
  rhoToken              = consumes<double>(cfg.getParameter<edm::InputTag>("rho"));
  recVtxsToken          = consumes<reco::VertexCollection>(cfg.getParameter<edm::InputTag>("vertices"));
  triggerResultsToken   = consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"));
  triggerPrescalesToken = consumes<pat::PackedTriggerPrescales>(cfg.getParameter<edm::InputTag>("triggerPrescales"));
  pupInfoToken          = consumes<edm::View<PileupSummaryInfo> >(edm::InputTag("slimmedAddPileupInfo"));
  genEvtInfoToken       = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
  genParticlesToken     = consumes<edm::View<reco::GenParticle> >(edm::InputTag("prunedGenParticles"));
  lheEvtInfoToken       = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
  runInfoToken          = consumes<LHERunInfoProduct>(edm::InputTag("externalLHEProducer"));
  //srcBtag_              = cfg.getParameter<std::string>("btagger");
  //  xmlFile_              = cfg.getParameter<std::string>("xmlFile");
  triggerNames_         = cfg.getParameter<std::vector<std::string> >("triggerNames");
  etaMax_               = cfg.getParameter<double>("etaMax");
  ptMin_                = cfg.getParameter<double>("ptMin");
  ptMinLeading_         = cfg.getParameter<double>("ptMinLeading");
  //massMin_              = cfg.getParameter<double>("massMin");
  //btagMin_              = cfg.getParameter<double>("btagMin");
  //minMuPt_              = cfg.getParameter<double>("minMuPt");
  //minElPt_              = cfg.getParameter<double>("minElPt");
  isMC_                 = cfg.getUntrackedParameter<bool>("isMC",false);
  isPrint_              = cfg.getUntrackedParameter<bool>("isPrint",false);
  saveWeights_          = cfg.getUntrackedParameter<bool>("saveWeights",true);
  debug_                = cfg.getUntrackedParameter<bool>("debug",false);
  GenptMin_             = cfg.getUntrackedParameter<double>("GenptMin");
  GenetaMax_            = cfg.getUntrackedParameter<double>("GenetaMax");

  jetFlavourInfosToken_ = consumes<reco::JetFlavourInfoMatchingCollection>( cfg.getParameter<edm::InputTag>("jetFlavourInfos"));
  
  //Gen Jet information
  fAKJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, 0.8);
  int activeAreaRepeats = 1;
  double ghostArea      = 0.01;
  double ghostEtaMax    = 7.0;
  fActiveArea           = new fastjet::ActiveAreaSpec (ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition       = new fastjet::AreaDefinition (fastjet::active_area_explicit_ghosts, *fActiveArea );
  sd = new fastjet::contrib::SoftDrop(0.0,0.1,0.8);//beta_, zCut_, R0 );
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducer::beginJob() 
{
  //--- book the trigger histograms ---------
  triggerNamesHisto_ = fs_->make<TH1F>("TriggerNames","TriggerNames",1,0,1);
  triggerNamesHisto_->SetCanExtend(TH1::kAllAxes);
  for(unsigned i=0;i<triggerNames_.size();i++) {
    triggerNamesHisto_->Fill(triggerNames_[i].c_str(),1);
  }
  triggerPassHisto_ = fs_->make<TH1F>("TriggerPass","TriggerPass",1,0,1);
  triggerPassHisto_->SetCanExtend(TH1::kAllAxes);

  cutFlowHisto_ = fs_->make<TH1F>("CutFlow","CutFlow",1,0,1);
  cutFlowHisto_->SetCanExtend(TH1::kAllAxes);
 
  //--- book the tree -----------------------
  outTree_ = fs_->make<TTree>("events","events");
  outTree_->Branch("runNo"                ,&run_               ,"run_/I");
  outTree_->Branch("evtNo"                ,&evt_               ,"evt_/I");
  outTree_->Branch("lumi"                 ,&lumi_              ,"lumi_/I");
  outTree_->Branch("nvtx"                 ,&nVtx_              ,"nVtx_/I");
  outTree_->Branch("nJets"                ,&nJets_             ,"nJets_/I");
  outTree_->Branch("pvRho"                ,&pvRho_             ,"pvRho_/F");
  outTree_->Branch("pvz"                  ,&pvz_               ,"pvz_/F");
  outTree_->Branch("pvchi2"               ,&pvchi2_            ,"pvchi2_/F");
  outTree_->Branch("pvndof"               ,&pvndof_            ,"pvndof_/F");
  outTree_->Branch("rho"                  ,&rho_               ,"rho_/F");
  outTree_->Branch("ht"                   ,&ht_                ,"ht_/F");
  //------------------------------------------------------------------
  flavor_         = new std::vector<int>;
  flavorHadron_   = new std::vector<int>;
  pt_             = new std::vector<float>;
  unc_            = new std::vector<float>;
  cor_            = new std::vector<float>;
  eta_            = new std::vector<float>;
  phi_            = new std::vector<float>;
  mass_           = new std::vector<float>;
  energy_         = new std::vector<float>;
  chf_            = new std::vector<float>;
  nhf_            = new std::vector<float>;
  phf_            = new std::vector<float>;
  muf_            = new std::vector<float>;
  elf_            = new std::vector<float>;
  chm_            = new std::vector<int>;
  nhm_            = new std::vector<int>;
  phm_            = new std::vector<int>;
  elm_            = new std::vector<int>;
  mum_            = new std::vector<int>;
  outTree_->Branch("jetPt"                ,"vector<float>"     ,&pt_);
  outTree_->Branch("jetCorr"              ,"vector<float>"     ,&cor_);
  outTree_->Branch("jetUnc"               ,"vector<float>"     ,&unc_);
  outTree_->Branch("jetEta"               ,"vector<float>"     ,&eta_);
  outTree_->Branch("jetPhi"               ,"vector<float>"     ,&phi_);
  outTree_->Branch("jetMass"              ,"vector<float>"     ,&mass_);
  outTree_->Branch("jetEnergy"            ,"vector<float>"     ,&energy_);
  outTree_->Branch("jetChf"               ,"vector<float>"     ,&chf_);
  outTree_->Branch("jetNhf"               ,"vector<float>"     ,&nhf_);
  outTree_->Branch("jetPhf"               ,"vector<float>"     ,&phf_);
  outTree_->Branch("jetMuf"               ,"vector<float>"     ,&muf_);
  outTree_->Branch("jetElf"               ,"vector<float>"     ,&elf_);
  outTree_->Branch("jetChm"               ,"vector<int>"       ,&chm_);
  outTree_->Branch("jetNhm"               ,"vector<int>"       ,&nhm_);
  outTree_->Branch("jetPhm"               ,"vector<int>"       ,&phm_);
  outTree_->Branch("jetElm"               ,"vector<int>"       ,&elm_);
  outTree_->Branch("jetMum"               ,"vector<int>"       ,&mum_);
  //------------------------------------------------------------------
  triggerBit_ = new std::vector<bool>;
  triggerPre_ = new std::vector<int>;
  outTree_->Branch("triggerBit"           ,"vector<bool>"      ,&triggerBit_);
  outTree_->Branch("triggerPre"           ,"vector<int>"       ,&triggerPre_);
  cout<<"Begin job finished"<<endl;
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducer::endJob() 
{  
  delete flavor_;
  delete flavorHadron_;
  delete pt_;
  delete cor_;
  delete unc_;
  delete btag_;
  delete eta_;
  delete phi_;
  delete mass_;
  delete energy_;
  delete chf_;
  delete nhf_;
  delete phf_;
  delete muf_;
  delete elf_;
  delete chm_;
  delete nhm_;
  delete phm_;
  delete elm_;
  delete mum_;
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
  if (isMC_ && debug_) {
    iRun.getByToken(runInfoToken,runInfo);
    for(vector<LHERunInfoProduct::Header>::const_iterator it = runInfo->headers_begin();it != runInfo->headers_end(); it++) {
      cout<<it->tag()<<endl;
      vector<string> lines = it->lines();
      for(unsigned int iLine = 0; iLine < lines.size(); iLine++) {
        cout<< lines.at(iLine);
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducer::endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{
}
//////////////////////////////////////////////////////////////////////////////////////////
bool BoostedTTbarFlatTreeProducer::isGoodJet(const pat::Jet &jet)
{
  bool res  = true; // by default is good, unless fails a cut bellow
/*  float chf = jet.chargedHadronEnergyFraction();
  float nhf = jet.neutralHadronEnergyFraction();
  float phf = jet.photonEnergyFraction();
  float muf = jet.muonEnergyFraction();
  float elf = jet.electronEnergyFraction();
  int chm   = jet.chargedHadronMultiplicity();
  int npr   = jet.neutralMultiplicity()+jet.chargedMultiplicity();
  float eta = fabs(jet.eta());
  float pt  = jet.pt();
  bool idL  = (npr>1 && phf<0.99 && nhf<0.99);
  //bool idM = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.99 && muf<0.99 && chf>0 && chm>0) || eta>2.4));
  bool idT = (idL && ((eta<=2.4 && nhf<0.9 && phf<0.9 && elf<0.9 && muf<0.9 && chf>0 && chm>0) || eta>2.4));
  //if (!idT) res = false;
  if (pt < ptMin_) res = false;
  if (eta > etaMax_) res = false;
  cout << "Jet Parameters "<< idT << " "<< pt <<" "<< eta << endl;
  //if (jet.userFloat("ak8PFJetsCHSSoftDropMass") < massMin_) res = false;
  //if ((jet.subjets("SoftDrop")).size() < 2) res = false;*/
  return res;
}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducer::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  initialize();

  if(isPrint_) cout<<"**** EVENT ****"<<endl;

  iEvent.getByToken(jetsToken,jets);
  iEvent.getByToken(candsToken,cands);
  iEvent.getByToken(rhoToken,rho);
  iEvent.getByToken(recVtxsToken,recVtxs);  
  iEvent.getByToken(triggerResultsToken,triggerResults);  
  iEvent.getByToken(triggerPrescalesToken,triggerPrescales); 

  triggerBit_->clear();
  triggerPre_->clear();

  //-------------- Trigger Info -----------------------------------
  triggerPassHisto_->Fill("totalEvents",1);
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);  
  bool passTrigger(false);
  for(unsigned int k=0;k<triggerNames_.size();k++) {
    bool bit(false);
    int pre(1);
    for(unsigned int itrig=0;itrig<triggerResults->size();itrig++) {
      //cout<<"Trigger name of index "<<itrig<<" "<<string(names.triggerName(itrig))<<endl;
      string trigger_name = string(names.triggerName(itrig));
      //--- erase the last character, i.e. the version number----
      trigger_name.pop_back();
      if (trigger_name == triggerNames_[k]) {
        bit = true;//triggerResults->accept(itrig);
        pre = triggerPrescales->getPrescaleForIndex(itrig);
        if (bit) {
          triggerPassHisto_->Fill(triggerNames_[k].c_str(),1);
        } 
      }
    }
    //--- if at least one monitored trigger has fired passTrigger becomes true
    passTrigger += bit;
    triggerBit_->push_back(bit); 
    triggerPre_->push_back(pre);   
  }   
  //----- at least one good vertex -----------
  bool cut_vtx = (recVtxs->size() > 0);
  if (cut_vtx) {
    pvRho_    = (*recVtxs)[0].position().Rho();
    pvz_    = (*recVtxs)[0].z();
    pvndof_ = (*recVtxs)[0].ndof();
    pvchi2_ = (*recVtxs)[0].chi2();
  }// if vtx
  //----- PF jets ------------------------------
  nJets_  = 0;
  nGenJets_  = 0;
  nBJets_ = 0;
  ht_     = 0.0;

  edm::ESHandle<JetCorrectorParametersCollection> PFJetCorParCollCHS;
  iSetup.get<JetCorrectionsRecord>().get("AK8PFchs",PFJetCorParCollCHS);
  JetCorrectorParameters const& PFJetCorParCHS = (*PFJetCorParCollCHS)["Uncertainty"];
  mPFUncCHS = new JetCorrectionUncertainty(PFJetCorParCHS);//"Summer16_23Sep2016V4_MC_Uncertainty_AK8PFchs.txt");

  double unc=0.0;
  
  vector<LorentzVector> vP4; 
  if(jets->size() > 0)
      cout << "Number of Jets "<< jets->size() <<" "<< jets->begin()->pt() <<  endl;
  for(pat::JetCollection::const_iterator ijet =jets->begin();ijet != jets->end(); ++ijet) {  
    if (isGoodJet(*ijet)) {
        cout << "Good jet " << ijet->pt() << endl;
        flavor_        ->push_back(ijet->partonFlavour());
	flavorHadron_  ->push_back(ijet->hadronFlavour());
        chf_           ->push_back(ijet->chargedHadronEnergyFraction());
        nhf_           ->push_back(ijet->neutralHadronEnergyFraction());
        phf_           ->push_back(ijet->photonEnergyFraction());
        elf_           ->push_back(ijet->electronEnergyFraction());
        muf_           ->push_back(ijet->muonEnergyFraction());
	chm_           ->push_back(ijet->chargedHadronMultiplicity());
	nhm_           ->push_back(ijet->neutralHadronMultiplicity());
	phm_           ->push_back(ijet->photonMultiplicity());
	elm_           ->push_back(ijet->electronMultiplicity());
	mum_           ->push_back(ijet->muonMultiplicity());
        pt_            ->push_back(ijet->pt());

	mPFUncCHS->setJetEta(ijet->eta());
	mPFUncCHS->setJetPt(ijet->pt()); // here you must use the CORRECTED jet pt
	unc = mPFUncCHS->getUncertainty(true);
	cor_           ->push_back(1+unc);
	unc_           ->push_back(unc);
        phi_           ->push_back(ijet->phi());
        eta_           ->push_back(ijet->eta());
        mass_          ->push_back(ijet->mass());
        energy_        ->push_back(ijet->energy());
        vP4.push_back(ijet->p4());
        ht_ += ijet->pt();
        nJets_++;
    }// if good jet
  }// jet loop       
  rho_    = *rho;
  nVtx_   = recVtxs->size();
  run_    = iEvent.id().run();
  evt_    = iEvent.id().event();
  lumi_   = iEvent.id().luminosityBlock();

  bool cut_RECO = (nJets_ >= 1);  
 
  cutFlowHisto_->Fill("All",1);
  if (iEvent.isRealData()) {
    if (cut_RECO) {
      cutFlowHisto_->Fill("At least one jet",1);
      if (passTrigger || 1) {
        cutFlowHisto_->Fill("Trigger",1);
        outTree_->Fill();
      }
    }
  } 
  else {
    if (cut_RECO) {
      cutFlowHisto_->Fill("At least one jet (reco || gen)",1);
      outTree_->Fill();
    }
  }

}
//////////////////////////////////////////////////////////////////////////////////////////
void BoostedTTbarFlatTreeProducer::initialize()
{
  run_            = -1;
  evt_            = -1;
  lumi_           = -1;
  nVtx_           = -1;
  nJets_          = -1;
  rho_            = -1;
  pvRho_          = -999;
  pvz_            = -999;
  pvndof_         = -999;
  pvchi2_         = -999;
  flavor_         ->clear();
  flavorHadron_   ->clear();
  pt_             ->clear();
  cor_             ->clear();
  unc_             ->clear();
  eta_            ->clear();
  phi_            ->clear();
  mass_           ->clear();
  energy_         ->clear();
  chf_            ->clear();
  nhf_            ->clear();
  phf_            ->clear();
  elf_            ->clear();
  muf_            ->clear();
  chm_            ->clear();
  nhm_            ->clear();
  phm_            ->clear();
  elm_            ->clear();
  mum_            ->clear();
}
//////////////////////////////////////////////////////////////////////////////////////////
BoostedTTbarFlatTreeProducer::~BoostedTTbarFlatTreeProducer() 
{
}

DEFINE_FWK_MODULE(BoostedTTbarFlatTreeProducer);
















