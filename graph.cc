#include <iostream>

using namespace std;

void graph(){
  int runnum = 14;
  int nentry;
  double charge;

  double integral_ct, integral_si;
  
  TGraph *graph = new TGraph();
  TF1 *func = new TF1("func", "pol1", 1, 100);

  func -> FixParameter(0, 0);

  TString filename_ct = Form("./process/run_%05d/run%05d_ch1.root", runnum, runnum);
  TString filename_si = Form("./process/run_%05d/run%05d_ch3.root", runnum, runnum);

  TFile *inputfile_ct = new TFile(filename_ct, "read");
  TFile *inputfile_si = new TFile(filename_si, "read");

  TTree *inputtree_ct = (TTree*)inputfile_ct -> Get("tree");
  TTree *inputtree_si = (TTree*)inputfile_si -> Get("tree");

  inputtree_ct -> SetBranchAddress("integral", &integral_ct);
  inputtree_si -> SetBranchAddress("integral", &integral_si);

  nentry = inputtree_ct -> GetEntries();
  for(int ientry = 0; ientry < nentry; ientry++){
    inputtree_ct -> GetEntry(ientry);
    inputtree_si -> GetEntry(ientry);
    charge = integral_ct * 4 * pow(2, -13) / 50;
    integral_si = integral_si/pow(10, 6);

    if(charge < 0.005){
      continue;
    }

    graph -> SetPoint(ientry, integral_si, charge);
    
  }

   graph -> Draw("AP");  
   graph -> Fit("func");

}
