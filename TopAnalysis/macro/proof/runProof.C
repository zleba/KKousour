void AddPythia(TChain *ch);
void AddMadgraph(TChain *ch);
#include <string>
using namespace std;

void runProof() { 

    TChain * chain = new TChain("ak4PUPPI/events");
    //AddPythia(chain);
    //AddMadgraph(chain);
    //return;
    //chain->Add("/nfs/dust/cms/user/zlebcr/JEC/histos/runGpuppi.root"); 
    chain->Add("/nfs/dust/cms/user/zlebcr/JEC/jets1.root"); 
    //bool on = nworkers>0;
    if (0) { 
        TProof::Open("workers=13");
        //TProof::Open(TString::Format("workers=%d",nworkers) );

        chain->SetProof(); 
    }

    TStopwatch w; w.Start();
    chain->Process("jecFiller.C");//, "", 1000000, 3716000);
    w.Print();

}

