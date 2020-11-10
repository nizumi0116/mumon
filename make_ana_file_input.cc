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
#include "TApplication.h"

#include "qualitycheck.h"

using namespace std;

int main(int argc, char *argv[]){

  TApplication app("app", 0, 0, 0, 0);
  
  int run = -1;
  int channel = -1;
  int type = -1;
  int daqid = -1;
  double factor = -1;

  int p = -1;
  while((p = getopt(argc, argv, "r:c:t:d:f:")) != -1)
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
      case 'f':
	factor = (double) atof(optarg);
	break;
      }
    }

  cout << "factor " << factor << "\n";

  if(run == -1 || channel == -1 || type == -1 || daqid == -1 || factor == -1)
    {
      cout << "!!! usage !!!" << '\n';
      cout << "-r : input a run number" << '\n';
      cout << "-c : input a channel number" << "\n";
      cout << "-t : CT = 0, Si = 1, EMT = 2" << "\n";
      cout << "-d : input a ID for daqpc (1, 2, 3)" << "\n";
      cout << "-f : enter factor if CT, f = 1\n";
      exit(0);
    }
  
  int startbin = 600;  
  int endbin = 1400; 
  int pedbin = 1400;
 
  //const int data = 1024; //number of samples
  const int data = 2050; //number of samples //Si array = 2049
  int range = 210000;

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
  //TString ofn = Form("./process/run_000%d/run000%d_ch16.root", run, run);
  TString ofn = Form("./process/run_000%d/run000%d_ch%d.root", run, run, channel);
  //TString ofn = Form("./process/position/%dmv_min_calc.root", run);
  TFile *fout = new TFile(ofn, "recreate");

  //define TTree
  TTree * tree = new TTree("tree", "");

  //define branches
  tree->Branch("entry", &entry, "entry/I");
  tree->Branch("array", array, "array[2050]/D");  //waveform
  tree->Branch("pedestal", &pedestal, "pedestal/D"); //pedestal
  tree->Branch("rms", &rms, "rms/D");
  tree->Branch("peak", &peak, "peak/D");
  tree->Branch("peak_bin", &peak_bin, "peak_bin/I");
  tree->Branch("sum", sum, "sum[2050]/D");
  tree->Branch("integral", &integral, "integral/D");
  tree->Branch("width", &width, "width/I");
  tree->Branch("sum_start", &sum_start, "sum_start/I");
  tree->Branch("sum_stop", &sum_stop, "sum_stop/I");

  //input filename
  TString filename;
  if (type == 0){
    filename = Form("./process/run_000%d/run000%d_ch%d.dat", run, run, channel);
  }
  else{
    filename = Form("./rawdata/run_000%d/run000%d_ch%d.txt", run, run, channel);
  }
  
  ifstream fin;

  fin.open(filename);

  TH1D *pedeHist;
  
    //for ( int itest = 0; itest < 5; itest++ ) {

  while ( ! fin.eof() ) //test 
  { 
	//calculate the width of the waveform
      
      double dummy = 0; 
      double value = 0;
      double min = 21000;
      double max = 0;

      int range_start = 0;
      int range_stop = 0;

      if(type == 0){
	pedeHist = new TH1D(Form("h_%d", num), "pedestal", 2000, -1000, 1000);
      }

      if(type != 0){
	pedeHist = new TH1D(Form("h_%d", num), "pedestal", 2000, 1000 * factor, 15000 * factor);
      }
	
      if(type == 1 || type == 2){
	for(int i = 0; i < 2; i++){
	  fin >> dummy;
	}
      }
      
      for(int i = 0; i < startbin; i++){
        fin >> dummy;
	dummy = dummy * factor;
	array[i] = dummy;

	if(type != 0){
	  pedeHist->Fill(dummy);
	}
	
      }
	
      for(int i = startbin; i < endbin; i++){  //where signal is
	fin >> value;
	value = value * factor;
	array[i] = value;
	
	if(type == 0  || type == 2){
	  if(value < min){
	    min = value;
	    peak_bin = i;
	  }
	}
	         	   
	if(type == 1 && value > max){
	  max = value;
	  peak_bin = i;	  
	}
	
      }
      
      for(int i = pedbin; i < data; i++){
	fin >> value;
	value = value * factor;
	array[i] = value;

	if(type == 0){
	  pedeHist->Fill(value);
	}
	
      }

      if(type == 0){
	range_start = -100;
	range_stop = 100;
      }
      
      if(type != 0){
	range_start = 1000 * factor;
	range_stop = 15000 * factor;
      }
      
      if(pedeHist->Fit("gaus", "Q", "", range_start, range_stop) == -1)
	{
	  delete pedeHist;
	  cout << "entry " << num << '\n';
	  cout << "Fit failed. This event is skipped." << '\n';
	  continue; //test
	}
      

      f = pedeHist -> GetFunction("gaus");
      pedestal = f ->GetParameter(1);
      sigma = f -> GetParameter(2);
      delete pedeHist; //test

      ///pedeHist -> Draw(); // test
      //app.Run(); //test

      
      if(type == 0 || type == 2){
	peak = abs(min - pedestal);
      }
      if(type == 1){
	peak = abs(max - pedestal);
      }
      rms = sigma;          
  
      /*
      if(!pedcheck(pedestal))
	{
	  pedestal = -1;
	  peak = -1;
	  integral = -1;
	  rms = -1;
	}
      */
      
      entry = num;
      num = num + 1;
      
      for(int i = peak_bin; i > startbin; i = i - 1){
	sum_start = i;
	if(type == 0 || type == 2){
	  if(array[i] - pedestal > -peak*0.1){
	    break;
	  }
	}
	if(type == 1 && array[i] - pedestal < peak*0.1){
	  break;
	}
      }

      for(int i = peak_bin; i < endbin; i++){
	sum_stop = i;
	if(type == 0 || type == 2){
	  if(array[i] - pedestal > -peak*0.1){
	    break;
	  }
	}
	if(type == 1 && array[i] - pedestal < peak*0.1){
	  break;
	}
      }

      width = sum_stop - sum_start;
      
      double charge = 0;

      for(int i = 0; i < data; i++){
	array[i] = array[i] - pedestal;
	
	if(i > sum_start-11 && i <  sum_stop+11){
	  if(type == 0 || type == 2){
	    if(array[i] < 0){
	      charge = charge + abs(array[i]);
	      sum[i] = sum[i-1] + abs(array[i]);
	    }
	  }
	 
	  if(type == 1 && array[i] > 0){
	    charge = charge + abs(array[i]);
	    sum[i] = sum[i-1] + abs(array[i]);
	  }
	  else{
	    sum[i] = sum[i-1];
	  }
	
	integral = charge;

	}
	
	}//loop over i //test
      
      if(num%100==0) cout << "processing : " << num << '\n';
      tree->Fill();
      
    }//while loop
  
  fout->cd();
  tree->Write();
  fout->Close();
  
  cout << "done!" << '\n';
  
}

