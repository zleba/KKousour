void AddPythia(TChain *ch);
void AddMadgraph(TChain *ch);
#include <string>
#include <cmath>

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


void runProof(int nMax=1, int nNow=0) { 

    
    //gSystem->Load("../../plugins/QCDjet.h");
    gROOT->ProcessLine(".L ../../plugins/QCDjet.h");


    TChain * chain = new TChain("ak4/events");
    //AddPythia(chain);
    //AddMadgraph(chain);
    //return;
    //chain->Add("/nfs/dust/cms/user/zlebcr/JEC/histos/runGpuppi.root"); 
    chain->Add("/nfs/dust/cms/user/zlebcr/JEC/ntuplesNewFormat/merged/jetsB.root"); 
    //chain->Add("/nfs/dust/cms/user/zlebcr/JEC/ntuplesNewFormat/merged/jetsC.root"); 
    //bool on = nworkers>0;

    //chain->LoadTree(-1);
    //int N = chain->GetEntries();

    if (1) { 
        TProof::Open("workers=10");//"workers=3");
        //TString connect = gSystem->GetFromPipe("pod-info -c");
        //TProof::Open(connect);//"workers=3");
        //
        //TProof::Open(TString::Format("workers=%d",nworkers) );

        chain->SetProof(); 
    }



    //int nStart = lround(nNow / double(nMax) * N);
    //int nEnd   = lround((nNow+1) / double(nMax) * N);


    //TStopwatch w; w.Start();
    //chain->Process("matching.C", "", nEnd - nStart + 1, nStart);
    chain->Process("matching.C+");//,"", -1, 0);
    //w.Print();

}

