#include <iostream>
#include <fstream>
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
  //enter this parameter//// 
  //////////////////////////
  //double efficiency = 0.83; //EMT

  double efficiency;
  cout << "enter efficiency \n";
  cin >> efficiency;
  
  TString high_run = "./high_intensity_si.dat";
  TString low_run = "./low_intensity_si.dat";

  ifstream fin[2];
  fin[0].open(high_run);
  fin[1].open(low_run);

  vector<int> high_runv, low_runv;

  int run = 0;
  
  while( ! fin[0].eof() ){
    fin[0] >> run;
    high_runv.push_back(run);
  }
  high_runv.push_back(-1);

  while( ! fin[1].eof() ){
    fin[1] >> run;
    low_runv.push_back(run);
  }
  low_runv.push_back(-1);

  cout << "low " << low_runv.size() << "\n";
  cout << "high " << high_runv.size() << "\n";
  cout << "total number of files = " << low_runv.size() + high_runv.size() -2 << "\n";
  
  double ratio_max = 200; // Si/CT
  double ratio_max_si = 2; // Si/ref. Si

  int days_max = 500;

  TH2D *history_plot[4];
  //history_plot[4] = new TH2D("history_plot4", "; amount of irradiation [# of e / 10^11]; EMT/CT (integral)", 2000, 0, 200, 1000, 0, ratio_max);
  history_plot[0] = new TH2D("history_plot0", "; entries; Si/input (integral)", 500, 0, 5000, 1000, 0, ratio_max);
  history_plot[1] = new TH2D("history_plot1", "; amount of irradiation [days in J-PARC @500 kW]; Si/input (integral)", 1500, 0, days_max, 1000, 0, ratio_max);
  history_plot[2] = new TH2D("history_plot2", "; entries; Si/ref. Si (integral)", 500, 0, 5000, 1000, 0, ratio_max_si);
  history_plot[3] = new TH2D("history_plot3", "; amount of irradiation [days in J-PARC @500 kW]; Si/ref. Si (integral)", 1500, 0, days_max, 1000, 0, ratio_max_si);
  
  TH1D *h0[2];
  TLine l1;

  const int nch = 3;
  TString filename;
  TFile *inputfile;
  TTree *inputtree;
  int entry, nentry;
  double peak[nch];
  double integral[nch];
  double irradiation = 0;
  double days = 0;
  vector<double> switching, switching_e, switching_low;
  
  int high = 0;
  int low = 0;
  int type; // 0: low, 1: high 
  int runnum;
  int low_entry = 0;
  int run_start = 0;
  int run_stop = 0;
  
  for( unsigned int ifile = 0; ifile < low_runv.size() + high_runv.size() -2; ifile++){

    if( low_runv[low] == -1 ){
      runnum = high_runv[high];
      type = 1;
      high++;
    }
    else if( high_runv[high] == -1 ){
      runnum = low_runv[low];
      type = 0;
      low++;
    }
    else if ( low_runv[low] < high_runv[high] ){
      runnum = low_runv[low];
      type = 0;
      low++;
    }
    else{
      runnum = high_runv[high];
      type = 1;
      high++;
    }

    cout << "analyze run " << runnum << "\n"; 
    if( ifile == 0 ){
      run_start = runnum;
    }
    if( ifile == low_runv.size() + high_runv.size() -2 ){
      run_stop = runnum;
    }
    
    filename = Form("./process/run_%05d/run_%05d.root", runnum, runnum);
    
    inputfile = new TFile(filename, "read");    
    inputtree = (TTree*)inputfile -> Get("tree");
    
    inputtree -> SetBranchAddress("entry", &entry);
    inputtree -> SetBranchAddress("peak", peak);
    inputtree -> SetBranchAddress("integral", integral);

    if(type == 0){
      h0[0] = new TH1D("si/ct", "; Si/CT", 1000, 0, ratio_max);
      h0[1] = new TH1D("si/ref. si", "; Si/ref. Si", 1000, 0, ratio_max_si);
    }
    
    nentry = inputtree -> GetEntries();
    for(int ientry = 0; ientry < nentry; ientry++){
      inputtree->GetEntry(ientry);
      
      irradiation = irradiation + ( integral[0] * 4 * 1/pow(2, 13) / 50 / 1.6 * pow(10, 10) / pow(10, 11));
      days = irradiation * pow(10, 11) * efficiency / (2.8 * pow(10, 7) / 2.48 * 60 * 60 * 24);

      if(type == 0){
	h0[0] -> Fill(integral[1]/integral[0]);
	//history_plot[4] -> Fill(irradiation, integral[2]/integral[0]);
	history_plot[0] -> Fill(low_entry, integral[1]/integral[0]);

	h0[1] -> Fill(integral[1]/integral[2]);
	//history_plot[2] -> Fill(irradiation, integral[2]/integral[1]);
	history_plot[2] -> Fill(low_entry, integral[1]/integral[2]);
	low_entry++;
      }
    } // loop over entries

    if(type == 1){
      switching_e.push_back(irradiation);
      switching.push_back(days);
    }

    if(type == 0){
      history_plot[1] -> Fill(days, h0[0]->GetMean());
      history_plot[3] -> Fill(days, h0[1]->GetMean());

      switching_low.push_back(low_entry);
    }

    //delete *h0;
    
  } // loop over files

  TApplication app("app", 0, 0, 0, 0);

  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 600);
  c0 -> Divide(2, 2);
  c0 -> cd(1);
  //history_plot[4] -> Draw("colz");
  history_plot[1] -> SetMarkerStyle(4);
  history_plot[1] -> Draw("P");
  for (unsigned int i = 0; i < switching.size(); i++){
    l1.DrawLine(switching[i], 0, switching[i], ratio_max);
  }
  
  c0 -> cd(2);
  //history_plot[0] -> Draw("colz");
  history_plot[3] -> SetMarkerStyle(4);
  history_plot[3] -> Draw("P");
  for (unsigned int i = 0; i < switching.size(); i++){
    l1.DrawLine(switching[i], 0, switching[i], ratio_max_si);
  }  

  c0 -> cd(3);
  history_plot[0] -> Draw("colz");
  for (unsigned int i = 0; i < switching_low.size(); i++){
    l1.DrawLine(switching_low[i], 0, switching_low[i], ratio_max);
  }
  
  c0 -> cd(4);
  history_plot[2] -> Draw("colz");
  for (unsigned int i = 0; i < switching_low.size(); i++){
    l1.DrawLine(switching_low[i], 0, switching_low[i], ratio_max_si);
  }
 
  app.Run();

  TFile ofn(Form("history_%dto%d.root", run_start, run_stop), "RECREATE");
  ofn.cd();;
  for(int i = 0; i < 4; i++){
    history_plot[i] -> Write();
  }
  ofn.Close();

  cout << "\n... Done!!\n";
  
  return 0;
  
}
