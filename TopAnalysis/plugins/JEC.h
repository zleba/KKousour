#ifndef JECs_h
#define JECs_h
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "SMPJ/AnalysisFW/interface/QCDJet.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

using namespace edm;
using namespace std;



class JECs  {

    typedef reco::Particle::LorentzVector LorentzVector;

    public:
    //JECs(bool IsMCarlo, string GlobalTag, string JETTYPE, string jecUncSrc, vector<string> jecUncSrcNames);

    //static bool sort_pfjets(QCDPFJet j1, QCDPFJet j2) {
        //return j1.ptCor() > j2.ptCor();
    //}


    JetCorrectorParameters *L1Fast, *L2Relative, *L3Absolute, *L2L3Residual;
    vector<JetCorrectorParameters> vecL1Fast, vecL2Relative, vecL3Absolute, vecL2L3Residual;
    FactorizedJetCorrector *jecL1Fast, *jecL2Relative, *jecL3Absolute, *jecL2L3Residual;
    JetCorrectionUncertainty *mPFUnc;

    JetCorrectorParameters *par;
    JetCorrectionUncertainty *tmpUnc;
    std::vector<JetCorrectionUncertainty*> mPFUncSrc;
    std::vector<string> mJECUncSrcNames;


    virtual TLorentzVector JEC_CHScorrections(bool IsMCarlo, const pat::Jet &jetNow, double  pfRho, vector<string> jecUncSrcNames, double &CorFactor, double &unc)
    {

        // ---- declaration of the vector of jets ---- //
        //vector<QCDPFJet> mPFJetsCHS; mPFJetsCHS.clear();

        TLorentzVector UnCorrectedJet = TLorentzVector(jetNow.px(), jetNow.py(), jetNow.pz(), jetNow.energy() );

        double pfjetArea = jetNow.jetArea();//pfjetchs.area());
        //double pfRho  = 0;//Event->evtHdr().pfRho()

        // ----- looping over the number of jets in the event ---- //

            //TLorentzVector UnCorrectedJet = tmpJet; ///Jets are already uncorrected
            // ---- Evaluating the L1Fast correction factor ---- //
            jecL1Fast->setJetPt(UnCorrectedJet.Pt());
            jecL1Fast->setJetA(pfjetArea);
            jecL1Fast->setRho(pfRho);
            jecL1Fast->setJetEta(UnCorrectedJet.Eta());
            //std::cout<<std::endl<<pfjetchs.area()<<std::endl;

            double corFactorL1Fast = jecL1Fast->getCorrection();
            //cout<<"L1Fast Cor Factor = "<<corFactorL1Fast<<endl;
            TLorentzVector tmpJetL1FastCorrected =
                UnCorrectedJet * corFactorL1Fast; // ---- getting the jet corrected at L1Fast level ----- //

            // ---- Evaluating the L2Relative correction factor ---- //
            jecL2Relative->setJetPt(tmpJetL1FastCorrected.Pt());
            jecL2Relative->setJetEta(tmpJetL1FastCorrected.Eta());

            double corFactorL2Relative = jecL2Relative->getCorrection();
            //cout<<"L2Relative Cor Factor"<<corFactorL2Relative<<endl;
            TLorentzVector tmpJetL1FastL2RelativeCorrected =
                tmpJetL1FastCorrected * corFactorL2Relative; //  ---- getting the jet corrected at L1Fast*L2Relative level ----- //

            // ---- Evaluating the L3Absolute correction factor ---- //
            jecL3Absolute->setJetPt(tmpJetL1FastL2RelativeCorrected.Pt());
            jecL3Absolute->setJetEta(tmpJetL1FastL2RelativeCorrected.Eta());

            double corFactorL3Absolute = jecL3Absolute->getCorrection();
            //cout<<"L3Absolute Cor Factor"<<corFactorL3Absolute<<endl;
            TLorentzVector tmpJetL1FastL2RelativeL3AbsoluteCorrected =
                tmpJetL1FastL2RelativeCorrected * corFactorL3Absolute; // -- getting the jet corrected at L1Fast*L2Relative*L3Absolute level -- //


            // ---- Evaluating the L2L3Rsidual correction factor ---- //
            //TLorentzVector tmpJetL1FastL2RelativeL3AbsoluteCorrected = UnCorrectedJet;

            TLorentzVector tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected;
            double corFactorL2L3Residual(-1.);
            if(!IsMCarlo) {
                jecL2L3Residual->setJetPt(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Pt());
                jecL2L3Residual->setJetEta(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Eta());

                corFactorL2L3Residual = jecL2L3Residual->getCorrection();
                //cout<<"L2L3Rsidual Cor Factor"<<corFactorL2L3Residual<<endl;
                tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected =
                    tmpJetL1FastL2RelativeL3AbsoluteCorrected * corFactorL2L3Residual; //  -- getting the jet corrected at L1Fast*L2Relative*L3Absolute*L2L3Residual level -- /
            } // if(!IsMCarlo)


            auto getUnc =[&](double pt, double eta) {
                mPFUnc->setJetEta(eta);
                mPFUnc->setJetPt(pt);
                return mPFUnc->getUncertainty(true);
            };


            LorentzVector correctedJetP4; CorFactor = 0;
            if(!IsMCarlo) { // -- if data, then take the full L1FastL2RelativeL3AbsoluteL2L3Residual corrected jet - //
                correctedJetP4 = LorentzVector(tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.Px(),
                        tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.Py(),
                        tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.Pz(),
                        tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.E());

                unc = getUnc(correctedJetP4.pt(), correctedJetP4.eta());
                CorFactor = tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected.E()/UnCorrectedJet.E();
                return tmpJetL1FastL2RelativeL3AbsoluteL2L3ResidualCorrected;
            }
            else { // -- if mc, the take the L1FastL2RelativeL3Absolute corrected jet --//
                correctedJetP4 = LorentzVector(tmpJetL1FastL2RelativeL3AbsoluteCorrected.Px(),
                        tmpJetL1FastL2RelativeL3AbsoluteCorrected.Py(),
                        tmpJetL1FastL2RelativeL3AbsoluteCorrected.Pz(),
                        tmpJetL1FastL2RelativeL3AbsoluteCorrected.E());

                unc = getUnc(correctedJetP4.pt(), correctedJetP4.eta());
                CorFactor = tmpJetL1FastL2RelativeL3AbsoluteCorrected.E()/UnCorrectedJet.E();
                return tmpJetL1FastL2RelativeL3AbsoluteCorrected;
            }
            return TLorentzVector();

            // ----- Storing the JEC correction factor in each label at a vector --- //

            //JecFactors.push_back(corFactorL1Fast); JecFactors.push_back(corFactorL2Relative);
            //JecFactors.push_back(corFactorL3Absolute); JecFactors.push_back(corFactorL2L3Residual);
            //JecFactors.push_back(CorFactor);

            //pfjetchs.setP4(correctedJetP4); // -- replacing the old P4 by newly corrected P4 -- //


            ///JEC Uncertainties

            /*
            //
            vector<float> uncSrc(0);
            for(unsigned isrc=0;isrc<jecUncSrcNames.size();isrc++) {
                //for(unsigned isrc=0;isrc<2;isrc++) {
                mPFUncSrc[isrc]->setJetEta(pfjetchs.eta());
                mPFUncSrc[isrc]->setJetPt(pfjetchs.pt());
                float unc1 = mPFUncSrc[isrc]->getUncertainty(true);
                uncSrc.push_back(unc1);
            } // for(unsigned isrc=0;isrc<mPFJECUncSrcNames.size();isrc++)
            ///end of JEC Uncertainties

            pfjetchs.setCor(CorFactor);  // -- replacing the old corfactor by the new corfactor -- //
            pfjetchs.setJecLabels(JecFactors); // ---- setting each label of JEC  ------ //
            pfjetchs.setUnc(unc);
            pfjetchs.setUncSrc(uncSrc);
            mPFJetsCHS.push_back(pfjetchs);

            } // for(int iJet = 0; iJet < NJets; iJet++)

            // ----- pT sorting the PF jets ------- //
            sort(mPFJetsCHS.begin(),mPFJetsCHS.end(),sort_pfjets);
            Event->setPFJetsCHS(mPFJetsCHS);
            */
        }


        void Init(bool IsMCarlo, string GlobalTag, string JETTYPE, string jecUncSrc, vector<string> mJECUncSrcNames){

            string file_data_mc = "DATA";
            if(IsMCarlo)    file_data_mc = "MC";

            string path = "/afs/desy.de/user/z/zlebcr/cms/CMSSW_9_3_0/src/KKousour/TopAnalysis/data/JECtablas/";
            cout << "JECpath: " << path+"/"+GlobalTag+"_"+file_data_mc+"/"+GlobalTag+"_"+file_data_mc+"_L1FastJet_"+JETTYPE+".txt" << endl;
            //exit(1);
            L1Fast       = new JetCorrectorParameters(path+"/"+GlobalTag+"_"+file_data_mc+"/"+GlobalTag+"_"+file_data_mc+"_L1FastJet_"+JETTYPE+".txt");
            L2Relative   = new JetCorrectorParameters(path+"/"+GlobalTag+"_"+file_data_mc+"/"+GlobalTag+"_"+file_data_mc+"_L2Relative_"+JETTYPE+".txt");
            L3Absolute   = new JetCorrectorParameters(path+"/"+GlobalTag+"_"+file_data_mc+"/"+GlobalTag+"_"+file_data_mc+"_L3Absolute_"+JETTYPE+".txt");
            if(!IsMCarlo)
                L2L3Residual = new JetCorrectorParameters(path+"/"+GlobalTag+"_"+file_data_mc+"/"+GlobalTag+"_"+file_data_mc+"_L2L3Residual_"+JETTYPE+".txt");

            vecL1Fast.push_back(*L1Fast);
            vecL2Relative.push_back(*L2Relative);
            vecL3Absolute.push_back(*L3Absolute);
            if(!IsMCarlo)
                vecL2L3Residual.push_back(*L2L3Residual);


            jecL1Fast       = new FactorizedJetCorrector(vecL1Fast);
            jecL2Relative   = new FactorizedJetCorrector(vecL2Relative);
            jecL3Absolute   = new FactorizedJetCorrector(vecL3Absolute);
            if(!IsMCarlo)
                jecL2L3Residual = new FactorizedJetCorrector(vecL2L3Residual);

            ///Read Uncertainty txt files
            mPFUnc = new JetCorrectionUncertainty(path+"/"+GlobalTag+"_"+file_data_mc+"/"+GlobalTag+"_"+file_data_mc+"_Uncertainty_"+JETTYPE+".txt");    

            for(unsigned isrc=0;isrc<mJECUncSrcNames.size();isrc++) {
                cout<<mJECUncSrcNames[isrc]<<endl;

                par = new JetCorrectorParameters(jecUncSrc,mJECUncSrcNames[isrc]);
                tmpUnc = new JetCorrectionUncertainty(*par);

                mPFUncSrc.push_back(tmpUnc);
            }
        }

    };


#endif
