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

#include <TStyle.h>
//#define test
#ifdef  test
#include "TROOT.h"
#endif

using namespace std;

const int CH = 32;
const int SAMPLE = 2049;
const TString OPTDIR = "./";

int main(int argc, char *argv[]) {

  #ifdef test
  gROOT->SetBatch();
  #endif

  TApplication app("app", 0, 0, 0, 0);

  /*
  TH1D **pedeHist; pedeHist = new TH1D*[CH];
  TH1D **Pedestal; Pedestal = new TH1D*[CH];
  for (int i = 0; i < CH; i++) {
    pedeHist[i] = new TH1D(Form("h_%d",i), "pedestal:ADCcounts:events", 1000, 2000, 3000);
    Pedestal[i] = new TH1D(Form("h__%d",i),"pedestal:ADCcounts:events", 1000, 2000, 3000);
  }
  */

  TH2D **map_all = 0;
  map_all = new TH2D*[2];
  for (int i = 0; i < 2; i++) map_all[i] = new TH2D(Form("map_all_%d", i), "si array map;x[mm];y[mm]", 8, -12, 12, 8, -12, 12);
  //TH2D* map_all = new TH2D("map_all", "si array map;x[mm];y[mm]", 8, -12, 12, 8, -12, 12);

  // 2D Gaussian Distribution
  TF2* g2 = new TF2("g2", "[0]*exp(- pow((x-[1]),2) / (2*pow([2],2)) ) * exp(- pow((y-[3]),2) / (2*pow([4],2)) )", -10, 10, -10, 10);
  g2->SetParameters(1, 0, 1, 0, 1);
  //

  if (argc != 3) {
    cout << "\nusage : " << argv[0] << "DAQPC1 run number & DAQPC2 run number\n";
    exit(0);
  }
  int runNumberDAQ1 = atoi(argv[1]);
  int runNumberDAQ2 = atoi(argv[2]);
  //int runNumberDAQ1 = 32, runNumberDAQ2 = 26; // test
  //cout << "select DAQPC==1== run number > " << flush;
  //cin  >> runNumberDAQ1;
  //cout << "select DAQPC==2== run number > " << flush;
  //cin  >> runNumberDAQ2;

  ifstream inputDAQ1[CH], inputDAQ2[CH];
  TString txtName;
  for (int i = 0; i < CH; i++) {
    txtName.Form("%srun000%d/run000%d_ch%d.txt", OPTDIR.Data(), runNumberDAQ1, runNumberDAQ1, i); // apply your own path
    inputDAQ1[i].open(txtName);
    txtName.Form("%srun000%d/run000%d_ch%d.txt", OPTDIR.Data(), runNumberDAQ2, runNumberDAQ2, i);
    inputDAQ2[i].open(txtName);
  }

  double xArray_R[CH] = {9,3,6,3,9,6,0,0,-3,-6,-9,-12,-6,-12,-9,-3,-12,-6,-6,-9,-12,-9,-3,-3,0,0,3,6,9,3,9,6}; // mm (view of down stream)
  double xArray_L[CH] = {-12,-6,-9,-6,-12,-9,-3,-3,0,3,6,9,3,9,6,0,9,3,3,6,9,6,0,0,-3,-3,-6,-9,-12,-6,-12,-9}; // mm edited 12/28/2019
  double yArray_R[CH] = {0,0,3,6,6,9,3,9,6,9,6,9,3,3,0,0,-3,-3,-9,-6,-9,-12,-6,-12,-3,-9,-12,-9,-12,-6,-6,-3}; // mm
  double yArray_L[CH] = {0,0,3,6,6,9,3,9,6,9,6,9,3,3,0,0,-3,-3,-9,-6,-9,-12,-6,-12,-3,-9,-12,-9,-12,-6,-6,-3}; // mm edited 12/28/2019

  int ADC_R[CH][SAMPLE]; double pedestal_R[CH]; double RMS_R[CH]; int integral_R[CH]; int ientry_R;
  int ADC_L[CH][SAMPLE]; double pedestal_L[CH]; double RMS_L[CH]; int integral_L[CH]; int ientry_L;

  TString optFileName = Form("%sprofile_daqpc1_run000%d.root", OPTDIR.Data(), runNumberDAQ1);
  TFile *fout_R = new TFile(optFileName, "RECREATE");
  TTree *profile_R = new TTree("profile_R","profile_R");
  profile_R->Branch("ADCcount", ADC_R,      "ADCcount[32][2049]/I");
  profile_R->Branch("pedestal", pedestal_R, "pedestal[32]/D");
  profile_R->Branch("RMS"     , RMS_R,      "RMS[32]/D");
  profile_R->Branch("integral", integral_R, "integral[32]/I");
  profile_R->Branch("ientry",   &ientry_R,  "ientry/I");

  int sumADC = 0;
  int x;
  int pedestalValue = 0;
  double pedestalRMS = 0;

  //daqpc1
  ientry_R = 0;
  while (!inputDAQ1[0].eof()) {
    for (int ich = 0; ich < CH; ich++) {
      sumADC = 0;
      pedestalValue = 0;
      pedestalRMS   = 0;
      for (int i = 0; i < 2; i++) {
	double dummy;
	inputDAQ1[ich] >> dummy;
      }
      for (int isample = 0; isample < SAMPLE; isample++) {
	inputDAQ1[ich] >> x;
	ADC_R[ich][isample] = x;
	sumADC += ADC_R[ich][isample];
	if (isample < 100) {
	  pedestalValue = pedestalValue + ADC_R[ich][isample];
	  pedestalRMS   = pedestalRMS + pow(ADC_R[ich][isample], 2.0) / 100;
	}
      }
      integral_R[ich] = sumADC;
      pedestal_R[ich] = pedestalValue / 100;
      RMS_R     [ich] = sqrt(pedestalRMS);
    }
    profile_R->Fill();
    ientry_R++;
  }

  // Write to TTree & Close
  profile_R->Write();
  fout_R->Close();
  //

  int ADcount_R[CH][SAMPLE]; double peds_R[CH]; double rms_R[CH]; int integ_R[CH]; int ient_R;
  TFile iptFile_R(optFileName);
  TTree *tree_R = (TTree*)iptFile_R.Get("profile_R");
  tree_R->SetBranchAddress("ADCcount", ADcount_R);
  tree_R->SetBranchAddress("pedestal", peds_R);
  tree_R->SetBranchAddress("RMS"     , rms_R);
  tree_R->SetBranchAddress("integral", integ_R);
  tree_R->SetBranchAddress("ientry"  , &ient_R);

  int nevt_R = tree_R->GetEntries();
  cout << "\nDAQPC1 event = " << nevt_R << "\n\n";

  double realInt[CH][nevt_R];
  TH2D *map_R = new TH2D("map_R", "map_R", 8, -12, 12, 8, -12, 12);
  TH2D **siArrayHist_R;
  siArrayHist_R = new TH2D*[nevt_R];
  for (int i = 0; i < nevt_R; i++) siArrayHist_R[i]=new TH2D(Form("h_%d",i),"si array;x [mm];y [mm]",8,-12,12,8,-12,12);

  for (int evt = 0; evt < nevt_R; evt++) {
    tree_R->GetEntry(evt);
    for (int ich = 0; ich < CH; ich++) {
      realInt[ich][evt] = integ_R[ich] - peds_R[ich]*SAMPLE;
      if (realInt[ich][evt] < 0) realInt[ich][evt] = 0;
    }
    for (int nch = 0; nch < CH; nch++) {
      for (int s = 0; s < realInt[nch][evt]/nevt_R; s++) { // what does it do?
	siArrayHist_R[evt]->Fill(xArray_R[nch], yArray_R[nch]);
	map_R->Add(siArrayHist_R[evt], 1);
	for (int i = 0; i < 2; i++) map_all[i]->Add(siArrayHist_R[evt], 1);
	siArrayHist_R[evt]->Reset();
      }
    }
  }

  TH1D *projectX_R = map_R->ProjectionX("pjx_R", -12, 12);
  TH1D *projectY_R = map_R->ProjectionY("pjy_R", -12, 12);

  optFileName.Form("%sprofile_daqpc2_run000%d.root", OPTDIR.Data(), runNumberDAQ2);
  TFile *fout_L = new TFile(optFileName, "RECREATE");
  TTree *profile_L = new TTree("profile_L","profile_L");
  profile_L->Branch("ADCcount", ADC_L,      "ADCcount[32][2049]/I");
  profile_L->Branch("pedestal", pedestal_L, "pedestal[32]/D");
  profile_L->Branch("RMS",      RMS_L,      "RMS[32]/D");
  profile_L->Branch("integral", integral_L, "integral[32]/I");
  profile_L->Branch("ientry",   &ientry_L,  "ientry/I");

  // daqpc2
  ientry_L = 0;
  while (!inputDAQ2[0].eof()) {
    for (int ich = 0; ich < CH; ich++) {
      sumADC = 0;
      pedestalValue = 0;
      pedestalRMS   = 0;
      for (int i = 0; i < 2; i++) {
	double dummy;
	inputDAQ2[ich] >> dummy;
      }
      for (int isample = 0; isample < SAMPLE; isample++) {
	inputDAQ2[ich] >> x;
	ADC_L[ich][isample] = x;
	sumADC += ADC_L[ich][isample];
	if (isample < 100) {
	  pedestalValue = pedestalValue + ADC_L[ich][isample];
	  pedestalRMS   = pedestalRMS + pow(ADC_L[ich][isample], 2.0) / 100;
	}
      }
      integral_L[ich] = sumADC;
      pedestal_L[ich] = pedestalValue / 100;
      RMS_L     [ich] = sqrt(pedestalRMS);
    }
    profile_L->Fill();
    ientry_L++;
  }

  // Write to TTree & Close
  profile_L->Write();
  fout_L->Close();
  //

  TFile iptFile_L(optFileName);
  int ADcount_L[CH][SAMPLE]; double peds_L[CH]; double rms_L[CH]; int integ_L[CH]; int ient_L;
  TTree *tree_L = (TTree*)iptFile_L.Get("profile_L");
  tree_L->SetBranchAddress("ADCcount", ADcount_L);
  tree_L->SetBranchAddress("pedestal", peds_L);
  tree_L->SetBranchAddress("RMS"     , rms_L);
  tree_L->SetBranchAddress("integral", integ_L);
  tree_L->SetBranchAddress("ientry"  , &ient_L);

  int nevt_L = tree_L->GetEntries();
  cout << "DAQPC2 event = " << nevt_L << "\n\n";
  double realInt1[CH][nevt_L];
  TH2D *map_L = new TH2D("map_L", "map_L", 8, -12, 12, 8, -12, 12);
  TH2D **siArrayHist_L; siArrayHist_L = new TH2D*[nevt_L];
  for (int i = 0; i < nevt_L; i++) siArrayHist_L[i] = new TH2D(Form("h_L%d",i),"si array;x [mm];y [mm]",8,-12,12,8,-12,12);

  for (int evt = 0; evt < nevt_L; evt++) {
    tree_L->GetEntry(evt);
    for (int ich = 0; ich < CH; ich++) {
      realInt1[ich][evt] = integ_L[ich] - peds_L[ich]*SAMPLE;
      if (realInt1[ich][evt] < 0) realInt1[ich][evt]=0;
    }
    for (int nch = 0; nch < CH; nch++) {
      for (int s = 0; s < realInt1[nch][evt]/nevt_L; s++) {
	siArrayHist_L[evt]->Fill(xArray_L[nch], yArray_L[nch]);
	map_L->Add(siArrayHist_L[evt], 1);
	for (int i = 0; i < 2; i++) map_all[i]->Add(siArrayHist_L[evt], 1);
	siArrayHist_L[evt]->Reset();
      }
    }
  }

  TH1D *projectX_L = map_L->ProjectionX("pjx_L", -12, 12);
  TH1D *projectY_L = map_L->ProjectionY("pjy_L", -12, 12);
  TH1D *projectX_all = map_all[0]->ProjectionX("pjx_all", -12, 12);
  TH1D *projectY_all = map_all[0]->ProjectionY("pjy_all", -12, 12);

  optFileName.Form("%soutput_%d_%d.root", OPTDIR.Data(), runNumberDAQ1, runNumberDAQ2);
  TFile optFile(optFileName, "RECREATE");
  optFile.cd();
  map_R->Write();
  projectX_R->Write();
  projectY_R->Write();
  map_L->Write();
  projectX_L->Write();
  projectY_L->Write();
  map_all[0]->Write();
  projectX_all->Write();
  projectY_all->Write();
  optFile.Close();

  // for offline
  TCanvas *c1 = new TCanvas("c1", "c1");
  c1->cd();
  map_all[0]->Draw("colz");
  map_all[0]->Fit(g2);
  TCanvas *c2 = new TCanvas("c2", "c2");
  c2->cd();
  map_all[1]->Draw("lego2z");

  cout << "\n=============================================\n";
  cout << "           Si Array Profile Result           \n";
  cout << "=============================================\n";
  cout << "\t   x    = " << g2->GetParameter(1) << " mm\n";
  cout << "\t   y    = " << g2->GetParameter(3) << " mm\n";
  cout << "\tsigma_x = " << g2->GetParameter(2) << " mm\n";
  cout << "\tsigma_y = " << g2->GetParameter(4) << " mm\n";
  cout << "=============================================\n\n";
  //

  // make backup data as pdf file
  TString fileName = Form("%stest.pdf", OPTDIR.Data());
  TCanvas *c5 = new TCanvas(fileName, fileName);
  c5->Print(fileName + "[", "pdf");
  map_R->Draw("colz");
  c5->Print(fileName, "pdf");
  projectX_R->Draw();
  c5->Print(fileName, "pdf");
  projectY_R->Draw();
  c5->Print(fileName, "pdf");
  map_L->Draw("colz");
  c5->Print(fileName, "pdf");
  projectX_L->Draw();
  c5->Print(fileName, "pdf");
  projectY_L->Draw();
  c5->Print(fileName, "pdf");
  map_all[0]->Draw("colz");
  map_all[0]->Fit(g2, "Q");
  //gStyle->SetOptStat(1000000001);
  gStyle->SetOptFit(1);
  c5->Print(fileName, "pdf");
  map_all[1]->Draw("lego2z");
  c5->Print(fileName, "pdf");
  projectX_all->Draw();
  c5->Print(fileName, "pdf");
  projectY_all->Draw();
  c5->Print(fileName, "pdf");
  c5->Print(fileName + "]", "pdf");
  delete c5;
  //

  cout << "\n... Back to terminal > Type 'ctl + C'\n";
  app.Run();
}
