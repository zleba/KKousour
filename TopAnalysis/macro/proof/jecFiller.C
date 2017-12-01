#define jecFiller_cxx
// The class definition in jecFiller.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("jecFiller.C")
// root> T->Process("jecFiller.C","some options")
// root> T->Process("jecFiller.C+")
//


#include "jecFiller.h"
#include <TH2.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>

void jecFiller::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void jecFiller::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   vector<double> ptEdges = {18,21,24,28,32,37,43,49,56,64,74,84,97,114,133,153,174,196,220,245,272,300,330,362,395,430,468,507,548,592,638,686,737,790,846,905,967,1032,1101,1172,1248,1327,1410,1497,1588, 1684,1784,1890,2000,2116,2238,2366,2500, 2640,2787,2941,3103,3273,3450,3637,3832,4037,4252,4477,4713,4961,5220,5492,5777,6076,6389,6717,7000 };

   hJetPt  = new TH1D("hJetPt", "hist", ptEdges.size()-1, ptEdges.data());
   fOutput->Add(hJetPt);

}

Bool_t jecFiller::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);

   for(unsigned i = 0; i < jetPt.GetSize(); ++i)
       hJetPt->Fill(jetPt[i]);



   return kTRUE;
}

void jecFiller::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void jecFiller::Terminate()
{
    hJetPt  = dynamic_cast<TH1D*>(fOutput->FindObject("hJetPt"));
    assert(hJetPt);
    TCanvas *can = new TCanvas("can", "can");
    hJetPt->Draw();
    can->SaveAs("radek.pdf");
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
