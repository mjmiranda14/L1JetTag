#include "L1JetTagPlotter.h"
#include <TTree.h>
#include <TFile.h>
#include <TChain.h>
#include "TH1F.h"
#include "TH2F.h"
#include <fstream>

int main(){

TChain* chain = new TChain("ntuple0/objects");
TString inpath = "root://cmsxrootd.fnal.gov///";

//Fill Sample file Chain
std::fstream s;
TString sampleName = "stop";
TString outFileName = sampleName+".root";
s.open(sampleName+".txt", std::fstream::in);
int n = 0;
if (s.is_open()) {
            TString x;
            while (!s.eof()) {
                s>>x;
                chain->Add(inpath+x);
                std::cout<<inpath<<x<<std::endl;
                n++;
            }
}
s.close();

L1JetTagPlotter S(chain);

std::vector<float> *vz = S.vz;

S.Loop(outFileName);






//TFile *f = TFile::Open(sampleName+".root", "recreate");
//f->cd();
//   h_pf_pt->Write();
//   h_pf_vz->Write();
//f->Close();

std::cout<<"Done with Final Loop"<<std::endl;

return 0;
}
