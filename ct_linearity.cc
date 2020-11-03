#include <iostream>
#include <vector>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TLegend.h"

using namespace std;

int main(int argc, char *argv[]){
  
  vector<TString> filelist_ct, filelist_si;
  vector<int> runnum;

  if (argc == 1){
    cout << "Error! Put the run numbers.\n";
    return 0;
  }
  
  for(int i = 1; i < argc; i++){
    cout << "analyze run " << atoi(argv[i]) << "\n";
    runnum.push_back(atoi(argv[i]));
    filelist_ct.push_back(Form("./process/run_%05d/run%05d_ch0.root", atoi(argv[i]), atoi(argv[i]) ));
    filelist_si.push_back(Form("./process/run_%05d/run%05d_ch1.root", atoi(argv[i]), atoi(argv[i]) ));
  }

  TString filename_ct, filename_si;
  TFile *inputfile_ct, *inputfile_si;
  TTree *inputtree_ct, *inputtree_si;
  int entry, nentry;
  double peak_ct, peak_si;
  double integral_ct, integral_si;
  int ct_range = 800;

  const int nfile = filelist_ct.size();
  TH2D *linearity_plot[nfile];
  TH2D *linearity_total = new TH2D("linearity_total", "; CT (integral); Si (integral)", ct_range, 0, ct_range, 10000, 0, 100000);
  TLegend *leg0 = new TLegend(0.7, 0.5, 0.9, 0.6);

  cout << "total number of files = " << nfile << "\n";
  
  for(int ifile = 0; ifile < nfile; ifile++){
    linearity_plot[ifile] = new TH2D(Form("linearity_plot_%d", ifile), "; CT (integral); Si (integral)", ct_range, 0, ct_range, 10000, 0, 100000);
    linearity_plot[ifile] -> SetMarkerStyle(8);
    linearity_plot[ifile] -> SetMarkerSize(0.5);
    linearity_plot[ifile] -> SetMarkerColor(ifile + 1);
    leg0 -> AddEntry(linearity_plot[ifile], Form("run%d", runnum[ifile]), "p");
    
    filename_ct = filelist_ct.at(ifile);
    filename_si = filelist_si.at(ifile);

    inputfile_ct = new TFile(filename_ct, "read");
    inputfile_si = new TFile(filename_si, "read");

    inputtree_ct = (TTree*)inputfile_ct->Get("tree");
    inputtree_si = (TTree*)inputfile_si->Get("tree");
    inputtree_ct->SetBranchAddress("entry", &entry);
    inputtree_ct->SetBranchAddress("peak", &peak_ct);
    inputtree_si->SetBranchAddress("peak", &peak_si);
    inputtree_ct->SetBranchAddress("integral", &integral_ct);
    inputtree_si->SetBranchAddress("integral", &integral_si);

    nentry = inputtree_ct -> GetEntries();
    for(int ientry = 0; ientry < nentry; ientry++){
      inputtree_ct->GetEntry(ientry);
      inputtree_si->GetEntry(ientry);

      //if(peak_bin > 680 && peak_bin < 700){ //if necessary
      linearity_plot[ifile] -> Fill(integral_ct, integral_si);
      linearity_total -> Fill(integral_ct, integral_si);
      //}
    }
    
  } // loop over files

  TApplication app("app", 0, 0, 0, 0);

  TF1 *func1 = new TF1("func1", "[0]*x+[1]", 0, ct_range);
  
  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 600);
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600);
  c0 -> cd();
  linearity_plot[0] -> Draw("P");
  for(int ifile = 1; ifile < nfile; ifile++){
    linearity_plot[ifile] -> Draw("P same");
  }
  leg0 -> Draw("same");

  c1 -> cd();
  linearity_total -> SetMarkerStyle(8);
  linearity_total -> SetMarkerSize(0.5);
  linearity_total -> SetMarkerColor(1);
  linearity_total -> Draw("P");
  linearity_total -> Fit("func1", "", "", 0, ct_range);
  
  c0->Update();
  c1->Update();
  app.Run();
  return 0;
}
