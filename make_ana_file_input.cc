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

#include "qualitycheck.h"

using namespace std;

int main(int argc, char *argv[]){

  int run = 0;

  int c = -1;
  while((c = getopt(argc, argv, "r:")) != -1)
    {
      switch(c){
      case 'r':
	run = atoi(optarg);
	break;
      }
    }

  if(run==0)
    {
      cout << "!!! usage !!!" << '\n';
      cout << "-r : input a run number" << '\n';
      exit(0);
    }
	
  int startbin = 600;  
  int endbin = 1400; 
  int pedbin = 1800;
  
  const int data = 2050; //number of samples
  const int range = 210000;

  int num = 0;
  int entry = 0;
  double array[data] = {};  //waveform

  //calculate the width of the waveform
  int sum_start = 0;
  int sum_stop = 0;
  int width = 0;

  //for pedestal calculation
  double pedestal;
  TF1 *f;
  double sigma, rms;
 
  //for charge calculation
  double integral;
  double peak;
  int peak_bin;
  double sum[data] = {};

  //TFile to save TTree
  TString ofn = Form("./process/run_000%d/run000%d_input.root", run, run);
  TFile *fout = new TFile(ofn, "recreate");

  //define TTree
  TTree * tree = new TTree("tree", "");

  //define branches
  tree->Branch("entry", &entry, "entry/I");
  tree->Branch("array", array, "array[2052]/D");  //waveform
  tree->Branch("pedestal", &pedestal, "pedestal/D"); //pedestal
  tree->Branch("rms", &rms, "rms/D");
  tree->Branch("peak", &peak, "peak/D");
  tree->Branch("peak_bin", &peak_bin, "peak_bin/I");
  tree->Branch("sum", sum, "sum[2052]/D");
  tree->Branch("integral", &integral, "integral/D");
  tree->Branch("width", &width, "width/I");
  tree->Branch("sum_start", &sum_start, "sum_start/I");
  tree->Branch("sum_stop", &sum_stop, "sum_stop/I");

  //input filename
  TString filename;
  filename = Form("./process/run_000%d/run000%d_ch0.dat", run, run);

  ifstream fin;

  fin.open(filename);

    //for ( int itest = 0; itest < 5; itest++ ) {

    
  while ( ! fin.eof() ) 
    {
	//calculate the width of the waveform
	
      double dummy = 0; 
      double value = 0;
      double min = 21000;
      
      TH1D* pedeHist = new TH1D(Form("h_%d", num), "pedestal", range, 0, range); 

      for(int i = 0; i < startbin; i++){
        fin >> dummy;
	array[i] = dummy;
      }      

      for(int i = startbin; i < endbin; i++){  //where signal is
	fin >> value;
	
	array[i] = value;
	
	  if(value < min){
	    min = value;
	    peak_bin = i;
	  }
      }

      for(int i = endbin; i < pedbin; i++){
	fin >> dummy;
	array[i] = dummy;
      }

      for(int i = pedbin; i < data; i++){
	fin >> value;
	
	array[i] = value;
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
	  peak = abs(min - pedestal);
	  rms = sigma;
	}
      if(!pedcheck(pedestal))
	{
	  pedestal = -1;
	  peak = -1;
	  integral = -1;
	  rms = -1;	  
	}      

      entry = num;
      num = num + 1;
      
      for(int i = peak_bin; i > startbin; i = i-1){
	sum_start = i;
	if(array[i] - pedestal > -peak*0.1){
	  break;
	}
      }

      for(int i = peak_bin; i < endbin; i++){
	sum_stop = i;
	if(array[i] - pedestal > -peak*0.1){
	  break;
	}
      }

      width = sum_stop - sum_start;
      
      double charge = 0;

      for(int i = 0; i < data; i++){
	array[i] = array[i] - pedestal;

	if(i > sum_start-11 && i <  sum_stop+11){
	  if(array[i] < 0){
	
	    charge = charge + abs(array[i]);
	    sum[i] = sum[i-1] + abs(array[i]);
	    
	  }
	  else{
	    sum[i] = sum[i-1];
	  }
	
	integral = charge;
	
	}
      } //loop over i
      
      if(num%100==0) cout << "processing : " << num << '\n';
      tree->Fill();

    }//while loop
  
  fout->cd();
  tree->Write();
  fout->Close();

  cout << "done!" << '\n';
 }

