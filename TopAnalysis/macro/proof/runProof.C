void AddPythia(TChain *ch);
void AddMadgraph(TChain *ch);
#include <string>
using namespace std;

void runForFile(TString fName) {
    TChain * chain = new TChain("ak4PUPPI/events");
    //AddPythia(chain);
    //AddMadgraph(chain);
    //return;
    //chain->Add("/nfs/dust/cms/user/zlebcr/JEC/histos/runGpuppi.root"); 
    chain->Add(fName); 
    //bool on = nworkers>0;
    if (1) { 
        TProof::Open("workers=13");
        //TProof::Open(TString::Format("workers=%d",nworkers) );

        chain->SetProof(); 
    }

    chain->Process("jecFiller.C");//, "", 1000000, 3716000);


}


void runProof(int nMax, int nNow) { 

    TChain * chain = new TChain("ak4PUPPI/events");
    //AddPythia(chain);
    //AddMadgraph(chain);
    //return;
    //chain->Add("/nfs/dust/cms/user/zlebcr/JEC/histos/runGpuppi.root"); 
    chain->Add("/nfs/dust/cms/user/zlebcr/JEC/ntuples2/merged/jetsB.root"); 
    chain->Add("/nfs/dust/cms/user/zlebcr/JEC/ntuples2/merged/jetsC.root"); 
    //bool on = nworkers>0;

    int N = chain->GetEntries();

    if (0) { 
        TProof::Open("workers=12");//"workers=3");
        //TProof::Open(TString::Format("workers=%d",nworkers) );

        chain->SetProof(); 
    }



    int nStart = lround(nNow / double(nMax) * N);
    int nEnd   = lround((nNow+1) / double(nMax) * N);


    //TStopwatch w; w.Start();
    chain->Process("/afs/desy.de/user/z/zlebcr/cms/CMSSW_8_0_29/src/KKousour/TopAnalysis/macro/proof/jecFiller.C", "", nEnd - nStart + 1, nStart);
    //w.Print();

}

