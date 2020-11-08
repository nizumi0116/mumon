#include <iostream>
#include <unistd.h>
#include <vector>
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TLegend.h"

using namespace std;

const int CH = 3;

int main(int argc, char *argv[]){

  if (argc == 1){
    cout << "Error! Put the run numbers.";
    return 0;
  }

  // l (0) : low intensity RUN
  // t (1) : 0.3 (three) nC / spill RUN
  // o (2) ;  1  (one)   nC / spill RUN
  int c = -1, totalRunNum[3] = {-1,-1,-1};
  while ((c = getopt(argc, argv, "l:t:o:")) != -1) {
    switch (c) {
    case 'l':
      totalRunNum[0] = atoi(optarg);
      break;
    case 't':
      totalRunNum[1] = atoi(optarg);
      break;
    case 'o':
      totalRunNum[2] = atoi(optarg);
      break;
    }
  }

  vector<vector<TString>> filelist(CH);
  vector<vector<int>> runnum(CH);

  // Loading kind of RUN files
  const int startArgc = 8;
  int loop = startArgc;
  TString anaRun[3] = {"Low Intensity", "0.3nC/spill", "1.0nC/spill"};
  cout << "\n====================================\n";
  cout << "        Analysis Information     \n";
  cout << "====================================\n\n";
  for (int i = 0; i < CH; i++) {
    int stopNum = loop + totalRunNum[i];
    cout << "    " << anaRun[i] << " <number of RUN> : " << totalRunNum[i] << endl;
    while (loop < stopNum) {
      cout << "    RUN number : " << atoi(argv[loop-1]) << endl;
      filelist[i].push_back(Form("~/mumon/process/run_%05d.root", atoi(argv[loop-1])));
      runnum[i].push_back(atoi(argv[loop-1]));
      loop++;
    }
    cout << "\n";
  }
  cout << "====================================\n\n";
  //

  TString filename;
  int *entry, nentry;
  double peak[CH];
  double integral[CH];
  const int ct_range = 800;
  const int si_range = 100000;

  TH2D ***linearity_plot = 0, ***linearity_si = 0;
  linearity_plot = new TH2D**[CH];
  linearity_si   = new TH2D**[CH];
  TH2D *linearity_total    = new TH2D("linearity_total", "; input (integral); EMT (integral)", ct_range, 0, ct_range, 1000, 0, 100000);
  TH2D *linearity_total_si = new TH2D("linearity_total_si", "; Si (integral); EMT (integral)", si_range, 0, si_range, 1000, 0, 100000);
  linearity_total->SetMarkerStyle(8);
  linearity_total->SetMarkerSize(0.5);
  linearity_total->SetMarkerColor(1);
  linearity_total_si->SetMarkerStyle(8);
  linearity_total_si->SetMarkerSize(0.5);
  linearity_total_si->SetMarkerColor(1);
  TLegend **leg_plot = 0, **leg_si = 0;
  leg_plot = new TLegend*[CH];
  leg_si   = new TLegend*[CH];

  loop = 0;
  for (int ich = 0; ich < CH; ich++) {
    linearity_plot[ich] = new TH2D*[totalRunNum[ich]];
    linearity_si  [ich] = new TH2D*[totalRunNum[ich]];

    leg_plot[ich] = new TLegend(0.7, 0.5, 0.9, 0.6);
    leg_si  [ich] = new TLegend(0.7, 0.5, 0.9, 0.6);

    for (int inum = 0; inum < totalRunNum[ich]; inum++) {
      linearity_plot[ich][inum] = new TH2D(Form("linearity_plot_ch_%d_%d", ich, inum), Form("%s;input (integral);EMT (integral)", anaRun[ich].Data()), ct_range, 0, ct_range, 1000, 0, 100000);
      linearity_si  [ich][inum] = new TH2D(Form("linearity_si_ch_%d_%d",   ich, inum), Form("%s;Si (integral);EMT (integral)", anaRun[ich].Data()),    si_range, 0, si_range, 1000, 0, 100000);
      linearity_plot[ich][inum]->SetMarkerStyle(8);
      linearity_plot[ich][inum]->SetMarkerSize(0.5);
      linearity_plot[ich][inum]->SetMarkerColor(ich + inum + 1);
      linearity_si  [ich][inum]->SetMarkerStyle(8);
      linearity_si  [ich][inum]->SetMarkerSize(0.5);
      linearity_si  [ich][inum]->SetMarkerColor(ich + inum + 1);
      leg_plot[ich]->AddEntry(linearity_plot[ich][inum], Form("run%d", runnum[ich][inum]), "p");
      leg_si  [ich]->AddEntry(linearity_si  [ich][inum], Form("run%d", runnum[ich][inum]), "p");
      filename = filelist[ich].at(inum);

      TFile *inputfile = new TFile(filename, "read");
      TTree *inputtree = (TTree*)inputfile->Get("tree");
      inputtree -> SetBranchAddress("entry", &entry);
      inputtree -> SetBranchAddress("peak", peak);
      inputtree -> SetBranchAddress("integral", integral);

      nentry = inputtree -> GetEntries();
      for(int ientry = 0; ientry < nentry; ientry++){
	inputtree->GetEntry(ientry);

	//if(peak_bin > 680 && peak_bin < 700){ //if necessary
	linearity_plot[ich][inum]->Fill(integral[0], integral[2]);
	linearity_si  [ich][inum]->Fill(integral[1], integral[2]);
	linearity_total    ->Fill(integral[0], integral[2]);
	linearity_total_si ->Fill(integral[1], integral[2]);
	//}
      }
      loop++;
      delete inputfile;
    }
  } // loop over files

  TApplication app("app", 0, 0, 0, 0);

  TF1 *func1 = new TF1("func1", "[0]*x+[1]", 0, ct_range);

  //TFile ofn(Form("emt_linearity_%dto%d.root", run_start, run_stop), "RECREATE");
  TFile ofn(Form("abc.root"), "RECREATE");
  ofn.cd();

  for (int ich = 0; ich < CH; ich++) {
    for (int inum = 0; inum < totalRunNum[ich]; inum++) {
      linearity_plot[ich][inum]->Write();
      linearity_si  [ich][inum]->Write();
    }
  }
  linearity_total->Write();
  linearity_total_si->Write();
  ofn.Close();

  TCanvas *c1 = new TCanvas("c1", "c1", 1500, 1500);
  c1->Divide(3, 2);
  int loopRunNum[CH] = {4, 5, 6};;
  for (int ich = 0; ich < CH; ich++) {
    if (totalRunNum[ich] == 0) continue;
    c1->cd(ich+1);
    linearity_plot[ich][0]->Draw("P");
    leg_plot[ich]->Draw("SAME");
    c1->cd(loopRunNum[ich]);
    linearity_si[ich][0]->Draw("P");
    leg_si[ich]->Draw("SAME");
    for (int inum = 1; inum < totalRunNum[ich]; inum++) {
      c1->cd(ich+1);
      linearity_plot[ich][inum]->Draw("SAME P");
      c1->cd(loopRunNum[ich]);
      linearity_si[ich][inum]->Draw("SAME P");
    }
    loop++;
  }
  c1->Print("sample.pdf");

  // if you want to get all run information
  /*
  TString canvasName("sample1.pdf");
  TCanvas *c2 = new TCanvas(canvasName, canvasName);
  c2->Print(canvasName + "[", "pdf");
  linearity_total->Draw("P");
  linearity_total->Fit(func1, "Q", "", 0, ct_range);
  c2->Print(canvasName, "pdf");
  linearity_total_si->Draw("P");
  linearity_total_si -> Fit(func1, "Q", "", 0, si_range);
  c2->Print(canvasName, "pdf");
  c2->Print(canvasName + "]", "pdf");
  delete c2;
  */
  //

  c1->Update();
  cout << "\n ... Back to terminal > type 'ctl + C'\n";

  app.Run();

  /*
  int run_start = 0;
  int run_stop  = 0;
  run_start = atoi(argv[1]);
  run_stop  = atoi(argv[argc - 1]);
  
  const int nch = 3;
  
  const int nfile = filelist.size();
  TH2D *linearity_plot[nfile];
  TH2D *linearity_si[nfile];
  TH2D *linearity_total = new TH2D("linearity_total", "; input (integral); EMT (integral)", ct_range, 0, ct_range, 1000, 0, 100000);
  TH2D *linearity_total_si = new TH2D("linearity_total_si", "; Si (integral); EMT (integral)", si_range, 0, si_range, 1000, 0, 100000);
  TLegend *leg0 = new TLegend(0.7, 0.5, 0.9, 0.6);
  TLegend *leg1 = new TLegend(0.7, 0.5, 0.9, 0.6);
  
  cout << "total number of files = " << nfile << "\n";
  
  for(int ifile = 0; ifile < nfile; ifile++){
    linearity_plot[ifile] = new TH2D(Form("linearity_plot_%d", ifile), "; input (integral); EMT (integral)", ct_range, 0, ct_range, 1000, 0, 100000);
    linearity_si[ifile] = new TH2D(Form("linearity_si_%d", ifile), "; Si (integral); EMT (integral)", si_range, 0, si_range, 1000, 0, 100000);
    linearity_plot[ifile] -> SetMarkerStyle(8);
    linearity_plot[ifile] -> SetMarkerSize(0.5);
    linearity_plot[ifile] -> SetMarkerColor(ifile + 1);
    
    linearity_si[ifile] -> SetMarkerStyle(8);
    linearity_si[ifile] -> SetMarkerSize(0.5);
    linearity_si[ifile] -> SetMarkerColor(ifile + 1);
    leg0 -> AddEntry(linearity_plot[ifile], Form("run%d", runnum[ifile]), "p");
    leg1 -> AddEntry(linearity_si[ifile], Form("run%d", runnum[ifile]), "p");
    
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
      linearity_plot[ifile] -> Fill(integral[0], integral[2]);
      linearity_si[ifile] -> Fill(integral[1], integral[2]);
      linearity_total -> Fill(integral[0], integral[2]);
      linearity_total_si -> Fill(integral[1], integral[2]);
      //}
    }
    
  } // loop over files

  TApplication app("app", 0, 0, 0, 0);
  
  TF1 *func1 = new TF1("func1", "[0]*x+[1]", 0, ct_range);
  
  TFile ofn(Form("emt_linearity_%dto%d.root", run_start, run_stop), "RECREATE");
  ofn.cd();
  
  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 600);
  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600);
  TCanvas *c2 = new TCanvas("c2", "c2", 1000, 600);
  TCanvas *c3 = new TCanvas("c3", "c3", 1000, 600);
  c0 -> cd();
  linearity_plot[0] -> Draw("P");
  linearity_plot[0] -> Write();
  leg0 -> Draw("same");
  for(int ifile = 1; ifile < nfile; ifile++){
    linearity_plot[ifile] -> Draw("P same");
    linearity_plot[ifile] -> Write();
  }
  
  c1 -> cd();
  linearity_si[0] -> Draw("P");
  linearity_si[0] -> Write();
  leg1 -> Draw("same");
  for(int ifile = 1; ifile < nfile; ifile++){
  linearity_si[ifile] -> Draw("P same");
  linearity_si[ifile] -> Write();
  }
  
  c2 -> cd();
  linearity_total -> SetMarkerStyle(8);
  linearity_total -> SetMarkerSize(0.5);
  linearity_total -> SetMarkerColor(1);
  linearity_total -> Draw("P");
  linearity_total -> Fit("func1", "", "", 0, ct_range);
  linearity_total -> Write();
  
  c3 -> cd();
  linearity_total_si -> SetMarkerStyle(8);
  linearity_total_si -> SetMarkerSize(0.5);
  linearity_total_si -> SetMarkerColor(1);
  linearity_total_si -> Draw("P");
  linearity_total_si -> Fit("func1", "", "", 0, si_range);
  linearity_total_si -> Write();
  
  c0->Update();
  c1->Update();
  c2->Update();
  c3->Update();
  
  ofn.Close();
  
  cout << "\n... Done!!\n";

  //app.Run();
  */
}
