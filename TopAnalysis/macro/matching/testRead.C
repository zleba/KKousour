#include "match.h"

void testRead()
{
    TFile *f = TFile::Open("ahoj.root");
    //TTree *tr = f->Get("Tr");


    //gROOT->ProcessLine(".L match.h+");

    TTreeReader Reader("Tr", f);

    TTreeReaderValue<Jet> pvec = {Reader, "jet"};


    while(Reader.Next() ) {
        cout << pvec->pT << endl;
        cout << pvec->PUPPI.size() << endl;

    }


}
