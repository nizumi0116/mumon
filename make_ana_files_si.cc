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

  int runstart = -1;
  int runend = -1;
  int daqid = -1;

  int c = -1;
  while((c = getopt(argc, argv, "s:e:d:")) != -1)
    {
      switch(c){
      case 's':
	runstart = atoi(optarg);
	break;
      case 'e':
	runend = atoi(optarg);
	break;
      case 'd':
	daqid = atoi(optarg);
	break;
      }
    }

  if(daqid==-1 || runstart==-1 || runend==-1)
    {
      cout << "!!! usage !!!" << '\n';
      cout << "-s : input a start run number" << '\n';
      cout << "-e : input a end run number" << '\n';
      cout << "-d : input a ID for daqpc (1, 2, 3)" << '\n';
      exit(0);
    }

  const int nrun = runend - runstart + 1;
  int startbin = 600;
  int endbin = 1400;
  
  const int data = 2050; //number of samples
  const int nch = 3;
  const int max = 2050;
  const int range = 17000;

  int ifile=0;
  int channel[nch];
  int num[nch] = {};
  int entry[nch] = {};
  double array[nch][max] = {};  //waveform

  //calculate the width
  int sum_start[nch];
  int sum_stop[nch];
  int width[nch];

  //for pedestal calculation
  double pedestal[nch];
  TF1 *f;
  double sigma[nch], rms[nch];
 
  //for charge calculation
  double integral[nch];
  double peak[nch];
  int peak_bin[nch];

  //TFile to save TTree
  TString ofn;
  TFile *fout[nrun];

  cout << "make " << nrun << " files" << '\n';

  for(int irun=runstart; irun<=runend; irun++)
    {

      for(int inum=0; inum<3; inum++) //initialize
	{
	  channel[inum] = 0;
	  num[inum] = 0;
	  entry[inum] = 0;	  
	}

      TTree * tree;
      cout << "make rootfile for run" << irun << '\n';
      
      //define branches
      ofn = Form("./process/run_000%d/run_000%d.root", irun, irun);
      fout[ifile] = new TFile(ofn, "recreate");
      tree = new TTree("tree", "tree");
      tree->Branch("channel", channel, "channel[3]/I");
      tree->Branch("entry", &entry, "entry/I");
      tree->Branch("array", array, "array[3][2050]/D");  //waveform
      tree->Branch("pedestal", pedestal, "pedestal[3]/D"); //pedestal
      tree->Branch("rms", rms, "rms[3]/D");
      tree->Branch("peak", peak, "peak[3]/D");
      tree->Branch("peak_bin", peak_bin, "peak[3]/I"); 
      tree->Branch("integral", integral, "integral[3]/D");
      tree->Branch("width", width, "width[3]/I");
      tree->Branch("sum_start", sum_start, "sum_start[3]/I");
      tree->Branch("sum_stop", sum_stop, "sum_stop[3]/I");
      
      //input filename
      TString filename[nch];
      filename[0] = Form("./rawdata/run_000%d/run000%d_ch0.dat", irun, irun);
      for(int ich = 1; ich < nch; ich++){
	filename[ich] = Form("./rawdata/run_000%d/run000%d_ch%i.txt", irun, irun, ich);
      }

      ifstream fin[nch];      
      for(int ich = 0; ich < nch; ich++){
	fin[ich].open(filename[ich]);
      }
      
      //for ( int itest = 0; itest < 5; itest++ ) {      
      
      while ( ! fin[ 0 ].eof() 
	      && ! fin[ 1 ].eof() 
	      && ! fin[ 2 ].eof() ) 
	{
	  
	  for(int ich = 0; ich < nch; ich++){
	    double dummy = 0.; 
	    double value = 0.;
	    double charge = 0;
	    double min = 17000;
	    double max = 0;
	    
	    TH1D* pedeHist = new TH1D(Form("h_ch%d_%d", ich, num[ich]), "pedestal", range, 0, range);

	    for(int i = 0; i < 2; i++){
	      fin[ich] >> dummy;
	    }
	    
	    for(int i = 0; i < startbin; i++){
	      fin[ich] >> dummy;
	      array[ich][i] = dummy;
	    }
	    
	    for(int i = startbin; i < endbin; i++){  //where signal is
	      fin[ich] >> value;
	      array[ich][i] = value;	      

	      if(ich == 0){
		if(value < min){
		  min = value;
		  peak_bin[ich] = i;
		}
	      }

	      if(ich == 1 || ich == 2){
		if(value > max){
		  max = value;
		  peak_bin[ich] = i;
		}
	      }
	      
	    }
	    
	    for(int i = endbin; i < data; i++){
	      fin[ich] >> value;
	      array[ich][i] = value;
	      pedeHist->Fill(value);
	    }

	    if(pedeHist->Fit("gaus","Q",0, range) == -1) 
	      {
		delete pedeHist;
		cout << "entry " << num[0] << '\n';
		cout << "Fit failed. This event is skipped." << '\n';
		continue;
	      }
	    f = pedeHist->GetFunction("gaus");
	    pedestal[ich] = f->GetParameter(1);
	    sigma[ich] = f->GetParameter(2);
	    delete pedeHist;
	    
	    if(pedcheck(pedestal[ich]))
	      {
		if(ich == 0){
		  peak[ich] = abs(min - pedestal[ich]);
		}
		
		if(ich == 1 || ich == 2){
		  peak[ich] = abs(max - pedestal[ich]);
		}

		rms[ich] = sigma[ich];
		
	      }

	    if(!pedcheck(pedestal[ich]))
	      {
		pedestal[ich] = -1;
		peak[ich] = -1;
		integral[ich] = -1;
		rms[ich] = -1;	  
	      }
	    
	    entry[ich] = num[ich];
	    num[ich] = num[ich] + 1;
	    channel[ich] = ich;

	    for(int i = peak_bin[ich]; i > startbin; i = i - 1){
	      sum_start[ich] = i;
	      if(ich == 0){
		if(array[ich][i] - pedestal[ich] > -peak[ich] * 0.1){
		  break;
		}
	      }

	      if(ich == 1 || ich == 2){
		if(array[ich][i] - pedestal[ich] < peak[ich] * 0.1){
		  break;
		}
	      }
	      
	    }

	    for(int i = peak_bin[ich]; i < endbin; i++){
	      sum_stop[ich] = i;
	      if(ich == 0){
		if(array[ich][i] - pedestal[ich] > -peak[ich] * 0.1){
		  break;
		}
	      }

	      if(ich == 1 || ich == 2){
		if(array[ich][i] - pedestal[ich] < peak[ich] * 0.1){
		  break;
		}
	      }
	      
	    }

	    width[ich] = sum_stop[ich] - sum_start[ich];

	    for(int i = 0; i < data; i++){
	      array[ich][i] = array[ich][i] - pedestal[ich];

	      if(i > sum_start[ich]-11 && i < sum_stop[ich]+11){
		if(ich == 0 && array[ich][i] < 0){
		    charge = charge + abs(array[ich][i]);
		  }
		
		if(ich == 1 || ich == 2){
		  if(array[ich][i] > 0){
		    charge = charge + abs(array[ich][i]);
		  }
		}

		integral[ich] = charge;
	      }
	    }
	    
	  } //loop ich
	  
	  if(num[0]%100==0) cout << "processing : " << num[0] << '\n';
	  tree->Fill();
	}//while loop
      
      fout[ifile]->cd();
      tree->Write();
      fout[ifile]->Close();

      //tree->Delete();
      fin[0].close();
      fin[1].close();
      fin[2].close();

      ifile++;
    } //loop irun    
  cout << "done!" << '\n';
}

