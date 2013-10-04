#include "FillHistograms.C"
#include "TString.h"

void FillHistogramsAll(TString SELECTION)
{
  FillHistograms("QCD100",SELECTION,true);
  FillHistograms("QCD250",SELECTION,true);
  FillHistograms("QCD500",SELECTION,true);
  FillHistograms("QCD1000",SELECTION,true);
  FillHistograms("ZJets",SELECTION,true);
  FillHistograms("WJets",SELECTION,true);
  FillHistograms("TTJets",SELECTION,true);
  FillHistograms("T_t-channel",SELECTION,true);
  FillHistograms("T_tW-channel",SELECTION,true);
  FillHistograms("T_s-channel",SELECTION,true);
  FillHistograms("Tbar_t-channel",SELECTION,true);
  FillHistograms("Tbar_tW-channel",SELECTION,true);
  FillHistograms("Tbar_s-channel",SELECTION,true);
  FillHistograms("VBFPowheg125",SELECTION,true);
  FillHistograms("GFPowheg125",SELECTION,true);
  if (SELECTION == "NOM") {
    FillHistograms("MultiJetA",SELECTION,false); 
    FillHistograms("BJetPlusXB",SELECTION,false);
    FillHistograms("BJetPlusXC",SELECTION,false);
    FillHistograms("BJetPlusXD",SELECTION,false);
  }
  if (SELECTION == "VBF") {
    FillHistograms("VBF1ParkedB",SELECTION,false);
    FillHistograms("VBF1ParkedC",SELECTION,false);
    FillHistograms("VBF1ParkedD",SELECTION,false);
  }
}
