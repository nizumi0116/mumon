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
  
  vector<TString> filelist;
  vector<int> runnum;

  if (argc == 1){
    cout << "Error! Put the run numbers.\n";
    return 0;
  }
  
  for(int i = 1; i < argc; i++){
    cout << "analyze run " << atoi(argv[i]) << "\n";
    runnum.push_back(atoi(argv[i]));
    filelist.push_back(Form("./process/run_%05d/run_%05d.root", atoi(argv[i]), atoi(argv[i]) ));
  }

  int run_start = 0;
  int run_stop = 0;
  run_start = atoi(argv[1]);
  run_stop = atoi(argv[argc - 1]);
  
  const int nch = 3;
  
  TString filename;
  TFile *inputfile;
  TTree *inputtree;
  int entry[nch], nentry;
  double peak[nch];
  double integral[nch];
  int ct_range = 800;
  int si_range = 100000;

  const int nfile = filelist.size();
  TH2D *linearity_plot[nfile];
  TH2D *linearity_total = new TH2D("linearity_total", "; Si (integral); input (integral)", si_range/10, 0, si_range, ct_range, 0, ct_range);
  TLegend *leg0 = new TLegend(0.7, 0.5, 0.9, 0.7);

  cout << "total number of files = " << nfile << "\n";
  
  for(int ifile = 0; ifile < nfile; ifile++){
    linearity_plot[ifile] = new TH2D(Form("linearity_plot_%d", ifile), "; input (integral); Si (integral)", si_range/10, 0, si_range, ct_range, 0, ct_range);
    linearity_plot[ifile] -> SetMarkerStyle(8);
    linearity_plot[ifile] -> SetMarkerSize(0.5);
    linearity_plot[ifile] -> SetMarkerColor(ifile + 1);
    leg0 -> AddEntry(linearity_plot[ifile], Form("run%d", runnum[ifile]), "p");
    
    filename = filelist.at(ifile);
    inputfile = new TFile(filename, "read");

    inputtree = (TTree*)inputfile->Get("tree");
    inputtree -> SetBranchAddress("entry", entry);
    inputtree -> SetBranchAddress("peak", peak);
    inputtree -> SetBranchAddress("integral", integral);

    nentry = inputtree -> GetEntries();
    for(int ientry = 0; ientry < nentry; ientry++){
      inputtree->GetEntry(ientry);

      //if(peak_bin > 680 && peak_bin < 700){ //if necessary
      linearity_plot[ifile] -> Fill(integral[1], integral[0]);
      linearity_total -> Fill(integral[1], integral[0]);
      //}
    }
    
  } // loop over files

  TApplication app("app", 0, 0, 0, 0);

  TF1 *func1 = new TF1("func1", "[0]*x+[1]", 0, si_range);

  TFile ofn(Form("ct_linearity_%dto%d.root", run_start, run_stop), "RECREATE");
  ofn.cd();
  
  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 600);
  c0 -> Divide(2, 1);
  c0 -> cd(1);
  linearity_plot[0] -> Draw("P");
  linearity_plot[0] -> Write();
  for(int ifile = 1; ifile < nfile; ifile++){
    linearity_plot[ifile] -> Draw("P same");
    linearity_plot[ifile] -> Write();
  }
  leg0 -> Draw("same");

  c0 -> cd(2);
  linearity_total -> SetMarkerStyle(8);
  linearity_total -> SetMarkerSize(0.5);
  linearity_total -> SetMarkerColor(1);
  linearity_total -> Fit("func1", "", "", 0, si_range);
  linearity_total -> Draw("P");
  linearity_total -> Write();
  
  c0->Update();

  ofn.Close();

  cout << "\n... Done!!\n";

  app.Run();
  
  return 0;
}
