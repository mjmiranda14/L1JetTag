#define L1JetTagPlotter_cxx
#include "L1JetTagPlotter.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void L1JetTagPlotter::Loop(TString outFileName)
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   nentries = 1000;  //Early Stopping
   std::cout<<"nentries="<<nentries<<std::endl;
   Long64_t nbytes = 0, nb = 0;

   h_pf_vz = new TH1F("pf_vz","pf_vz", 30,-15,15);
   h_pf_pt = new TH1F("pf_pt","pf_pt", 100,0,500);

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      for(unsigned i=0; i<pf->size(); i++){
       TLorentzVector v = (pf->at(i)).first;
       h_pf_pt->Fill(v.Pt());
      }

      for(unsigned i=0; i<pf_vz->size(); i++){
      h_pf_vz->Fill(pf_vz->at(i));
      }
   }
   TFile *outFile = 0;
   outFile = TFile::Open(outFileName, "RECREATE");
   outFile->cd();
   h_pf_vz->Write();
   h_pf_pt->Write();
   outFile->Close();
}
