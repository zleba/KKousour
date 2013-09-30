#ifndef VZFlatTreeProducer_h
#define VZFlatTreeProducer_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionData.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionEvaluator.h"
#include "HLTrigger/HLTcore/interface/TriggerExpressionParser.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "TTree.h"
#include "TH1F.h"
#include "KKousour/CMGAnalysis/plugins/QGLCalculator.h"

class VZFlatTreeProducer : public edm::EDAnalyzer 
{
  public:
    typedef reco::Particle::LorentzVector LorentzVector;
    explicit VZFlatTreeProducer(edm::ParameterSet const& cfg);
    virtual void beginJob();
    virtual void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup);
    virtual void endJob();
    virtual ~VZFlatTreeProducer();

  private:  
    void initialize();
    //---- configurable parameters --------   
    edm::InputTag srcJets_,srcMET_,srcElectrons_,srcMuons_,srcRho_;
    edm::Service<TFileService> fs_;
    std::string srcBtag_,srcPU_;
    TTree *outTree_; 
    //---- TRIGGER -------------------------
    triggerExpression::Data triggerCache_;
    std::vector<triggerExpression::Evaluator*> vtriggerSelector_;
    std::vector<std::string> vtriggerAlias_,vtriggerSelection_;
    TH1F *triggerPassHisto_,*triggerNamesHisto_;
    double minJetPt_,maxJetEta_,minJetPtSel_,maxJetEtaSel_,minBJetPtSel_,maxBJetEtaSel_,minElectronPt_,minMuonPt_;
    //---- output TREE variables ------
    int run_,evt_,lumi_,nVtx_,njets_,nbjetsSel_,njetsSel_,nelectrons_,nmuons_,npu_;
    float rho_,sphericity_,aplanarity_;
    std::vector<bool> *triggerResult_;
    //---- electrons -------------------------
    std::vector<float> *elPt_,*elEta_,*elPhi_,*elE_,*elIso_,*elMva_;
    std::vector<int>   *elCh_;
    std::vector<bool>  *elEB_,*elEE_,*elIdL_,*elIdM_,*elIdT_;
    std::vector<LorentzVector> *elP4_;
    //---- muons -------------------------
    std::vector<float> *muPt_,*muEta_,*muPhi_,*muE_,*muIso_;
    std::vector<int>   *muId_,*muCh_;
    std::vector<LorentzVector> *muP4_;
    //---- jets -------------------------
    std::vector<bool>  *jetPuIdL_,*jetPuIdM_,*jetPuIdT_;
    std::vector<float> *jetPt_,*jetEta_,*jetPhi_,*jetE_,*jetChf_,*jetNhf_,*jetPhf_,*jetMuf_,*jetElf_,*jetQGL_;
    std::vector<float> *jetBtag_,*jetPuMva_; 
    std::vector<int>   *jetId_;
    std::vector<LorentzVector> *jetP4_;
    //---- di-electrons -------------------------
    float eePt_,eeEta_,eePhi_,eeE_,eeM_;
    //---- di-muons -------------------------
    float mmPt_,mmEta_,mmPhi_,mmE_,mmM_;
    //---- electron-muon -------------------------
    float emPt_,emEta_,emPhi_,emE_,emM_;
    //---- di-jets -------------------------
    float jjPt_,jjEta_,jjPhi_,jjE_,jjM_,jjDPhi_,jjDEta_;
    //---- MET -------------------------
    float met_,metPhi_;
    //---- lljj -----------------------
    float eejjPt_,eejjY_,eejjM_,mmjjPt_,mmjjY_,mmjjM_;
    //---- QGL tagger -----------------
    QGLCalculator *QGL_;
};

#endif
