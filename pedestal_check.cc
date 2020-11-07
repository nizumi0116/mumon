#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "TApplication.h"

#include "qualitycheck.h"

using namespace std;

int main(int argc, char *argv[]){

  int run = -1;
  int channel = -1;
  int type = -1;
  int daqid = -1;

  int p = -1;
  while((p = getopt(argc, argv, "r:c:t:d:")) != -1)
    {
      switch(p){
      case 'r':
	run = atoi(optarg);
	break;
      case 'c':
	channel = atoi(optarg);
	break;
      case 't':
	type = atoi(optarg);
	break;
      case 'd':
	daqid = atoi(optarg);
	break;
      }
    }

  if(run == -1 || channel == -1 || type == -1 || daqid == -1)
    {
      cout << "!!! usage !!!" << '\n';
      cout << "-r : input a run number" << '\n';
      cout << "-c : input a channel number" << "\n";
      cout << "-t : CT = 0, Si, EMT = 1" << "\n";
      cout << "-d : input a ID for daqpc (1, 2, 3)" << "\n";
      exit(0);
    }
 
  //const int data = 1024; //number of samples
  const int data = 2050; //number of samples //Si array = 2049
  const int range = 210000;

  int num = 0;
  int entry = 0;

  //for pedestal calculation
  double pedestal, pedestal_mean;
  int ped_bin;
  TF1 *f, *f_mean;
  double sigma, rms;
  TGraph *ped_check = new TGraph();
  TGraph *rms_check = new TGraph();
  TH1D *pedeHist;
  TH1D *pedmean;

  //TFile to save TTree
  TString ofn = Form("./process/run_000%d/run000%d_ch%d_pedestal.root", run, run, channel);
  TFile *fout = new TFile(ofn, "recreate");

  //define TTree
  TTree * tree = new TTree("tree", "");

  //define branches
  tree->Branch("entry", &entry, "entry/I");
  tree->Branch("ped_bin", &ped_bin, "ped_bin/I");
  tree->Branch("pedestal", &pedestal, "pedestal/D"); //pedestal
  tree->Branch("sigma", &sigma, "sigma/D");

  //input filename
  TString filename;
  if(type == 1){
    filename = Form("./rawdata/run_000%d/run000%d_ch%d.txt", run, run, channel);
  }
  if(type == 0){
    filename = Form("./process/run_000%d/run000%d_ch%d.dat", run, run, channel);
  }
  
  ifstream fin;

  int iter = 0;
  
  for(int pedbin = 0; pedbin < data; pedbin = pedbin + 100){
  
    fin.open(filename);
    pedmean = new TH1D(Form("h_mean%d", pedbin), "; pedestal", range, 0, range);
   
    for ( int itest = 0; itest < 100; itest++ ) {
    
    //while ( ! fin.eof() ) {
	//calculate the width of the waveform
	
	double dummy = 0; 
	double value = 0;
	
	pedeHist = new TH1D(Form("h_%d", num), "pedestal", range, 0, range); 
	
	if(type == 1){
	  for(int i = 0; i < 2; i++){
	    fin >> dummy;
	  }
	}
	
	for(int i = 0; i < pedbin; i++){
	  fin >> dummy;
	}      
	
	for(int i = pedbin; i < data; i++){
	  fin >> value;
	  pedeHist->Fill(value);	  
	}
	
	if(pedeHist->Fit("gaus","Q",0, range) == -1)
	  {
	    delete pedeHist;
	    cout << "entry " << num << '\n';
	    cout << "Fit failed. This event is skipped." << '\n';
	    continue;
	  }
	
	f = pedeHist->GetFunction("gaus");
	pedestal = f->GetParameter(1);
	sigma = f->GetParameter(2);
	delete pedeHist;
	
	if(pedcheck(pedestal))
	  {
	    pedmean -> Fill(pedestal);
	  }
	if(!pedcheck(pedestal))
	  {
	    pedestal = -1;
	    sigma = -1;	  
	  }      
	
	entry = num;
	num = num + 1;
	
	tree->Fill();
	
    }//loop over itest

    fin.close();
    
    pedmean->Fit("gaus", "" , 0, range);
    f_mean = pedmean -> GetFunction("gaus");
    pedestal_mean = f_mean -> GetParameter(1);
    ped_check -> SetPoint(iter, pedbin, pedestal_mean);
    rms = f_mean -> GetParameter(2);
    rms_check -> SetPoint(iter, pedbin, rms);
    delete pedmean;

    iter++;
    
    cout << "analyzed pedestal bin = " << pedbin << " pedestal: " << pedestal << " RMS: " << rms << "\n";
    
  } // loop over pedbin
  
  fout->cd();
  tree->Write();
  fout->Close();

  TApplication app("app", 0, 0, 0, 0);
  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 600);
  c0 ->Divide(2, 1);
  c0 -> cd(1);
  ped_check -> SetMarkerStyle(3);
  ped_check -> GetXaxis() -> SetTitle("bin number");
  ped_check -> GetYaxis() -> SetTitle("pedestal");
  ped_check -> Draw("AP");
  c0 -> cd(2);
  rms_check -> SetMarkerStyle(3);
  rms_check -> GetXaxis() -> SetTitle("bin number");
  rms_check -> GetYaxis() -> SetTitle("RMS");
  rms_check -> Draw("AP");
  
  cout << "done!" << '\n';

  app.Run();
}
