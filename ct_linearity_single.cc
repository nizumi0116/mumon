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
#include "TMath.h"

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
    filelist_ct.push_back(Form("./process/run_%05d/run_%05d.root", atoi(argv[i]), atoi(argv[i]) ));
    filelist_si.push_back(Form("./process/run_%05d/run_%05d.root", atoi(argv[i]), atoi(argv[i]) ));
  }

  int run_start = 0;
  int run_stop = 0;
  run_start = atoi(argv[1]);
  run_stop = atoi(argv[argc - 1]);
  
  TString filename_ct, filename_si;
  TFile *inputfile_ct, *inputfile_si;
  TTree *inputtree_ct, *inputtree_si;;
  int entry, nentry;
  double peak;
  double integral_ct, integral_si;
  int ct_range = 800;
  int si_range = 100000;
  int charge_range = 10;
  double charge; //nC

  const int nfile = filelist_ct.size();
  TH2D *linearity_plot[nfile];
  TH2D *charge_hist[nfile];
  TLegend *leg0 = new TLegend(0.7, 0.5, 0.9, 0.7);

  cout << "total number of files = " << nfile << "\n";
  
  for(int ifile = 0; ifile < nfile; ifile++){
    linearity_plot[ifile] = new TH2D(Form("linearity_plot_%d", ifile), "; Si (integral); input (integral)", si_range/10, 0, si_range, ct_range, 0, ct_range);
    charge_hist[ifile] = new TH2D(Form("charge_hist_%d", ifile), "; charge [nC]; Si (integral)", charge_range * 1000, 0, charge_range, si_range/10, 0, si_range);
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
    inputtree_ct -> SetBranchAddress("entry", &entry);
    inputtree_ct -> SetBranchAddress("peak", &peak);
    inputtree_ct -> SetBranchAddress("integral", &integral_ct);
    inputtree_si -> SetBranchAddress("integral", &integral_si);

    nentry = inputtree_ct -> GetEntries();
    for(int ientry = 0; ientry < nentry; ientry++){
      inputtree_ct->GetEntry(ientry);
      inputtree_si->GetEntry(ientry);

      //if(peak_bin > 680 && peak_bin < 700){ //if necessary
      linearity_plot[ifile] -> Fill(integral_si, integral_ct);
      charge = integral_ct * TMath::Power(10, -9) * 4 * pow(2, -13) / 50 * TMath::Power(10, 9);
      charge_hist[ifile] -> Fill(charge, integral_si);
      //}
    }
    
  } // loop over files

  TApplication app("app", 0, 0, 0, 0);

  TF1 *func1 = new TF1("func1", "[0]*x+[1]", 0, si_range);

  TFile ofn(Form("ct_linearity_%dto%d.root", run_start, run_stop), "RECREATE");
  ofn.cd();
  
  TF1 *f1;
  
  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 600);
  c0 -> Divide(2, 1);
  c0 -> cd(1);
  linearity_plot[0] -> Fit("pol1");
  f1 = linearity_plot[0] -> GetFunction("pol1");
  linearity_plot[0] -> Draw("P");
  linearity_plot[0] -> Write();
  cout << "ref. Si = ( " << f1 -> GetParameter(1) << " +/- " << f1 -> GetParError(1) << " ) CT + ( " << f1 -> GetParameter(0) << " +/- " << f1 -> GetParError(0) << " )\n";
  
  for(int ifile = 1; ifile < nfile; ifile++){
    linearity_plot[ifile] -> Fit("pol1");
    f1 = linearity_plot[ifile] -> GetFunction("pol1");
    linearity_plot[ifile] -> Draw("P same");
    cout << "ref. Si = ( " << f1 -> GetParameter(1) << " +/- " << f1 -> GetParError(1) << " ) + ( " << f1 -> GetParameter(0) << " +/- " << f1 -> GetParError(0) << " )\n";
    linearity_plot[ifile] -> Write();
  }
  leg0 -> Draw("same");

  c0 -> cd(2);
  
  charge_hist[0] -> Fit("pol1");
  f1 = charge_hist[0] -> GetFunction("pol1");
  charge_hist[0] -> Draw("P");
  charge_hist[0] -> Write();
  
  for(int ifile = 1; ifile < nfile; ifile++){
    charge_hist[ifile] -> Fit("pol1");
    f1 = charge_hist[ifile] -> GetFunction("pol1");
    charge_hist[ifile] -> Draw("P same");
    charge_hist[ifile] -> Write();
  }
  leg0 -> Draw("same");
  
  c0->Update();

  ofn.Close();

  cout << "\n... Done!!\n";

  app.Run();
  
  return 0;
}
