#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TF2.h"
#include "TString.h"
#include "TTree.h"
#include "TGraph.h"
#include "math.h"
#include "TApplication.h"

using namespace std;

const int CH = 32;
const int SAMPLE = 2049;

int main(int argc, char *argv[]) {
  //void si_array(){

  TApplication app("app", 0, 0, 0, 0);
  TCanvas *c0 = new TCanvas("c0", "c0", 500, 500);
  
  //TH1D *pedeHist = new TH1D("h","pedestal:ADCcounts:events",500,2000,2500);
  TH1D **pedeHist; pedeHist = new TH1D*[CH];
  for(int p=0;p<CH;p++) {
    pedeHist[p] = new TH1D(Form("h_%d",p),"pedestal:ADCcounts:events",1000,2000,3000);
  }
  TH1D **Pedestal; Pedestal = new TH1D*[CH];
  for(int q=0;q<CH;q++) {
    Pedestal[q] = new TH1D(Form("h__%d",q),"pedestal:ADCcounts:events",1000,2000,3000);
  }
  
  TString txtname; TString txtname1;

  TH2D* map_all = new TH2D("map_all", "map_all", 8, -12, 12, 8, -12, 12);
  TF2* g2 = new TF2("g2", "[0]*exp(- pow((x-[1]),2) / (2*pow([2],2)) ) * exp(- pow((y-[3]),2) / (2*pow([4],2)) )", -10, 10, -10, 10); 
  //"[0]*gaus(x, [1], [2])*gaus(y, [3], [4])", -10, 10, -10, 10);
  g2->SetParameters(1, 0, 1, 0, 1);

  int ADC_R[CH][SAMPLE]; double pedestal_R[CH]; double RMS_R[CH]; int integral_R[CH]; int ientry_R; 
  int ADC_L[CH][SAMPLE]; double pedestal_L[CH]; double RMS_L[CH]; int integral_L[CH]; int ientry_L; 
  TTree *profile_R = new TTree("profile_R","profile_R");
  profile_R->Branch("ADCcount", ADC_R,      "ADCcount[32][2049]/I");
  profile_R->Branch("pedestal", pedestal_R, "pedestal[32]/D");
  profile_R->Branch("RMS"     , RMS_R,      "RMS[32]/D");
  profile_R->Branch("integral", integral_R, "integral[32]/I");
  profile_R->Branch("ientry",   &ientry_R,  "ientry/I");
    
  double xArray_R[32] = {9,3,6,3,9,6,0,0,-3,-6,-9,-12,-6,-12,-9,-3,-12,-6,-6,-9,-12,-9,-3,-3,0,0,3,6,9,3,9,6}; // mm (view of down stream)
  double xArray_L[32] = {-12,-6,-9,-6,-12,-9,-3,-3,0,3,6,9,3,9,6,0,9,3,3,6,9,6,0,0,-3,-3,-6,-9,-12,-6,-12,-9}; // mm edited 12/28/2019
  double yArray_R[32] = {0,0,3,6,6,9,3,9,6,9,6,9,3,3,0,0,-3,-3,-9,-6,-9,-12,-6,-12,-3,-9,-12,-9,-12,-6,-6,-3}; // mm
  double yArray_L[32] = {0,0,3,6,6,9,3,9,6,9,6,9,3,3,0,0,-3,-3,-9,-6,-9,-12,-6,-12,-3,-9,-12,-9,-12,-6,-6,-3}; //edited 12/28/2019

  ifstream fin[CH]; ifstream fin1[CH];

  int runNum; int runNum1;
  cout << "select DAQPC==1== run number >" << flush; cin >> runNum;
  cout << "select DAQPC==2== run number >" << flush; cin >> runNum1;

  for(int i=0;i<CH;i++) {
    txtname.Form("../profile/run000%d/run000%d_ch%d.txt", runNum, runNum, i);
    fin[i].open(txtname); 
    txtname1.Form("../profile/run000%d/run000%d_ch%d.txt", runNum1, runNum1, i);
    fin1[i].open(txtname1);
  }

  int Sum; int x; int cont=0; int pede; double pedeRMS=0;

  //daqpc1
  ientry_R = 0; 
  while(!fin[0].eof()) {
    for(int channel=0; channel<CH; channel++) {
      Sum=0;pede=0;pedeRMS=0;
      for(int i=0;i<2;i++) {
	double dummy;
	fin[channel] >> dummy;
      }
      for(int isamp=0; isamp<SAMPLE; isamp++) {
	fin[channel] >> x;
	ADC_R[channel][isamp] = x;
	Sum += ADC_R[channel][isamp];
	if(isamp<100) {
	  pede =  pede + ADC_R[channel][isamp];
	  pedeRMS = pedeRMS + pow(ADC_R[channel][isamp],2.0)/100;
	}
      }
      integral_R[channel] = Sum; 
      pedestal_R[channel] = pede/100; //dataout(pedestal[channel]);
      RMS_R[channel] = sqrt(pedeRMS); //dataout(sqrt(pedeRMS));//dataout(RMS[channel]);
    }
    profile_R->Fill(); cont++; ientry_R++;
    //if(ientry==2000)break; //debug
  }

  TString ofn;
  ofn.Form("../profile/profile_daqpc1_run000%d.root", runNum);
  //TFile er(ofn);
  //rofile->Write();
  //er.Close();
  profile_R->SaveAs(ofn);
  //profile_R->Write();
  //tree1->Reset();

  TFile fi(ofn);
  int ADcount_R[CH][SAMPLE]; double peds_R[CH]; double rms_R[CH]; int integ_R[CH]; int ient_R;
  TTree *tre = (TTree*)fi.Get("profile_R");
  tre->SetBranchAddress("ADCcount", ADcount_R);
  tre->SetBranchAddress("pedestal", peds_R);
  tre->SetBranchAddress("RMS"     , rms_R);
  tre->SetBranchAddress("integral", integ_R);
  tre->SetBranchAddress("ientry"  , &ient_R);

  int nevt = tre->GetEntries();
  cout << "nevt = " << nevt << "\n";

  double realInt[CH][nevt];
  TH2D *map = new TH2D("map","map",8,-12,12,8,-12,12);
  TH2D **siArrayHist_R; siArrayHist_R = new TH2D*[nevt];
  for(int i=0;i<nevt;i++) {
    siArrayHist_R[i]=new TH2D(Form("h_%d",i),"si array;x [mm];y [mm]",8,-12,12,8,-12,12);
  }

  for(int evt=0; evt<nevt; evt++) {
    tre->GetEntry(evt);
    for(int ich=0; ich<CH; ich++) {
      realInt[ich][evt] = integ_R[ich] - peds_R[ich]*SAMPLE;
      //dataout(realInt[ich][evt]);
      if(realInt[ich][evt] <0 ){
	realInt[ich][evt]=0;
      }
    }
    for(int nch=0;nch<CH;nch++) {
      for(int s=0;s<realInt[nch][evt]/nevt;s++) {
	siArrayHist_R[evt]->Fill(xArray_R[nch],yArray_R[nch]);
	map->Add(siArrayHist_R[evt],1);
	map_all->Add(siArrayHist_R[evt], 1);
	siArrayHist_R[evt]->Reset();
	}
    }
  }
  TH1D *projectX_R = map->ProjectionX("pjx",-12,12);
  TH1D *projectY_R = map->ProjectionY("pjy",-12,12);            

  TString fileName; fileName.Form("../profile/hoge_%d_%d.root", runNum, runNum1);
  TFile f1(fileName, "RECREATE");
  f1.cd();
  //siArrayHist_R->Write();
  map->Write();
  projectX_R->Write();
  projectY_R->Write();
  f1.Close();

  fi.Close();

  TTree *profile_L = new TTree("profile_L","profile_L");
  profile_L->Branch("ADCcount", ADC_L,      "ADCcount[32][2049]/I");
  profile_L->Branch("pedestal", pedestal_L, "pedestal[32]/D");
  profile_L->Branch("RMS"     , RMS_L,      "RMS[32]/D");
  profile_L->Branch("integral", integral_L, "integral[32]/I");
      profile_L->Branch("ientry",   &ientry_L,  "ientry/I");

      // daqpc2
      ientry_L = 0; 
      while(!fin1[0].eof()) {
	for(int channel=0; channel<CH; channel++) {
	  Sum=0;pede=0;pedeRMS=0;
	  for(int i=0;i<2;i++) {
	    double dummy;
	    fin1[channel]>>dummy;
	  }
	  for(int isamp=0; isamp<SAMPLE; isamp++) {
	    fin1[channel]>>x;
	    ADC_L[channel][isamp] = x;
	    Sum += ADC_L[channel][isamp];
	    if(isamp<100) {
	      pede =  pede + ADC_L[channel][isamp];
	      pedeRMS = pedeRMS + pow(ADC_L[channel][isamp],2.0)/100;
	    }
	  }
	  integral_L[channel] = Sum; 
	  pedestal_L[channel] = pede/100;
	  RMS_L[channel] = sqrt(pedeRMS); //dataout(RMS_L[channel]);
	}
	profile_L->Fill(); cont++; ientry_L++;
      }
        
      ofn.Form("../profile/profile_daqpc2_run000%d.root", runNum1);
      //TFile iu(ofn);
      //profile_L->Write();
      //iu.Close();
      profile_L->SaveAs(ofn);
      //profile_L->Write();
      
      TFile fi1(ofn);
      int ADcount_L[CH][SAMPLE]; double peds_L[CH]; double rms_L[CH]; int integ_L[CH]; int ient_L;
      TTree *tee = (TTree*)fi1.Get("profile_L");
      tee->SetBranchAddress("ADCcount", ADcount_L);
      tee->SetBranchAddress("pedestal", peds_L);
      tee->SetBranchAddress("RMS"     , rms_L);
      tee->SetBranchAddress("integral", integ_L);
      tee->SetBranchAddress("ientry"  , &ient_L);

      int nevt1 = tee->GetEntries();
      cout << "nevt1 = " << nevt1 << "\n";
      double realInt1[CH][nevt1];
      TH2D *map1 = new TH2D("map1","map1",8,-12,12,8,-12,12);
      TH2D **siArrayHist_L; siArrayHist_L = new TH2D*[nevt1];
      for(int i=0;i<nevt1;i++) {
	siArrayHist_L[i]=new TH2D(Form("h_L%d",i),"si array;x [mm];y [mm]",8,-12,12,8,-12,12);
      }

      for(int evt=0; evt<nevt1; evt++) {
	tee->GetEntry(evt);
	for(int ich=0; ich<CH; ich++) {
	  realInt1[ich][evt] = integ_L[ich] - peds_L[ich]*SAMPLE;
	  if(realInt1[ich][evt]<0) realInt1[ich][evt]=0;
	}
	for(int nch=0;nch<CH;nch++) {
	  for(int s=0;s<realInt1[nch][evt]/nevt1;s++) {
	    siArrayHist_L[evt]->Fill(xArray_L[nch],yArray_L[nch]);
	    map1->Add(siArrayHist_L[evt],1);
	    map_all->Add(siArrayHist_L[evt], 1);
	    siArrayHist_L[evt]->Reset();
	  }
	}
      }

      //map_all->Fit(g2);

      //map->Add(map1,1); //problem
      TH1D *projectX_L = map1->ProjectionX("pjx1",-12,12);
      TH1D *projectY_L = map1->ProjectionY("pjy1",-12,12);            
    
      //TString fileName = Form("output_%d_%d.root", runNum, runNum1);
      TFile f2(fileName, "UPDATE");
      f2.cd();
      map1->Write();
      map_all->Write();
      c0->cd();
      map_all->Draw("colz");
      c0->Update();
      app.Run();
      //g2->Draw("same");
      //g2->Write();
      //siArrayHist_R->Write();
      //map->Write();
      //projectX_R->Write();
      //projectY_R->Write();
      projectX_L->Write();
      projectY_L->Write();
      f2.Clear();
}
