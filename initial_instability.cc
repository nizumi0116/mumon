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
#include "TLine.h"

using namespace std;

int main(int argc, char *argv[]){

  //////////////////////////
  //enter this parameter// 
  //////////////////////////
  double efficiency = 0.83; //EMT
  
  vector<TString> filelist;

  TLine l1;
  double ratio_max = 200;
  double ratio_max_si = 10;

  if (argc == 1){
    cout << "Error! Put the run numbers.";
    return 0;
  }

  for(int i = 1; i < argc; i++){
    cout << "analyze run " << atoi(argv[i]) << "\n";
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
  double irradiation = 0;
  double days = 0;
  vector<double> switching, switching_e;

  const int nfile = filelist.size();
  TH2D *history_plot[3];
  history_plot[0] = new TH2D("history_plot0", "; amount of irradiation [# of e / 10^11]; EMT/CT (integral)", 2000, 0, 200, 10000, 0, ratio_max);
  history_plot[1] = new TH2D("history_plot1", "; amount of irradiation [days in J-PARC @500 kW]; EMT/CT (integral)", 1500, 0, 15, 10000, 0, ratio_max);
  history_plot[2] = new TH2D("history_plot2", "; amount of irradiation [days in J-PARC @500 kW]; EMT/CT (integral)", 1500, 0, 15, 10000, 0, ratio_max);
  TH2D *history_plot_si[3];
  history_plot_si[0] = new TH2D("history_plot_si0", "; amount of irradiation [# of e / 10^11]; EMT/Si (integral)", 2000, 0, 200, 10000, 0, ratio_max_si);
  history_plot_si[1] = new TH2D("history_plot_si1", "; amount of irradiation [days in J-PARC @500 kW]; EMT/Si (integral)", 1500, 0, 15, 10000, 0, ratio_max_si);
  history_plot_si[2] = new TH2D("history_plot_si2", "; amount of irradiation [days in J-PARC @500 kW]; EMT/Si (integral)", 1500, 0, 15, 10000, 0, ratio_max_si);

  TH1D *h0, *h0_si;
  
  cout << "total number of files = " << nfile << "\n";

  for(int ifile = 0; ifile < nfile; ifile++){

    h0 = new TH1D("h0", "; EMT/CT", 10000, 0, ratio_max);
    h0_si = new TH1D("h0_si", "; EMT/Si", 10000, 0, ratio_max_si);
    
    filename = filelist.at(ifile);
    inputfile = new TFile(filename, "read");
    inputtree = (TTree*)inputfile->Get("tree");
    
    inputtree->SetBranchAddress("entry", entry);
    inputtree->SetBranchAddress("peak", peak);
    inputtree->SetBranchAddress("integral", integral);

    nentry = inputtree -> GetEntries();
    for(int ientry = 0; ientry < nentry; ientry++){
      inputtree->GetEntry(ientry);

      irradiation = irradiation + ( integral[0] * 4 * 1/pow(2, 13) / 50 / 1.6 * pow(10, 10) / pow(10, 11));
      days = irradiation * pow(10, 11) * efficiency / (2.8 * pow(10, 7) / 2.48 * 60 * 60 * 24);

      h0 -> Fill(integral[2]/integral[0]);
      //history_plot[0] -> Fill(irradiation, integral[2]/integral[0]);
      history_plot[1] -> Fill(days, integral[2]/integral[0]);
      cout << "days " << days  << " " << integral[2]/integral[0] << "\n"; 

      h0_si -> Fill(integral[2]/integral[1]);
      //history_plot_si[0] -> Fill(irradiation, integral[2]/integral[1]);
      history_plot_si[1] -> Fill(days, integral[2]/integral[1]);
    }

    switching_e.push_back(irradiation);
    switching.push_back(days);

    history_plot[2] -> Fill(days, h0->GetMean());
    history_plot_si[2] -> Fill(days, h0_si->GetMean());
    
    delete h0;
    delete h0_si;
    
  } // loop over files

  TApplication app("app", 0, 0, 0, 0);

  TFile ofn(Form("initial_instability_%dto%d.root", run_start, run_stop), "RECREATE");
  ofn.cd();
  
  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 600);
  c0 -> cd();
  //history_plot[0] -> Draw("colz");
  history_plot[2] -> SetMarkerStyle(4);
  history_plot[2] -> Draw("P");
  for (unsigned int i = 0; i < switching.size(); i++){
    l1.DrawLine(switching[i], 0, switching[i], ratio_max);
  }
  history_plot[2] -> Write();
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600);
  c1 -> cd();
  //history_plot_si[0] -> Draw("colz");
  history_plot_si[2] -> SetMarkerStyle(4);
  history_plot_si[2] -> Draw("P");
  for (unsigned int i = 0; i < switching.size(); i++){
    l1.DrawLine(switching[i], 0, switching[i], ratio_max_si);
  }
  history_plot_si[2] -> Write();

  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 600);
  c2 -> cd();
  history_plot[1] -> Draw("colz");
  for (unsigned int i = 0; i < switching.size(); i++){
    l1.DrawLine(switching[i], 0, switching[i], ratio_max);
  }
  history_plot[1] -> Write();
  
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 600);
  c3 -> cd();
  history_plot_si[1] -> Draw("colz");
  for (unsigned int i = 0; i < switching.size(); i++){
    l1.DrawLine(switching[i], 0, switching[i], ratio_max_si);
  }
  history_plot_si[1] -> Write();

  app.Run();

  ofn.Close();

  cout << "\n... Done!!\n";
  
  return 0;
  
}
