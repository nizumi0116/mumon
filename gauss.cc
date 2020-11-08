#include <iostream>
#include "TMath.h"

using namespace std;

TCanvas* c0 = new TCanvas("c0");

void gauss(){
  double sigma_x = 0.17;
  double sigma_y = 0.17;
  double num = 9;
  double min, max;
  //double gauss1, gauss2, sum;

  TGraph *eff_g = new TGraph();
  TGraph *per_si_g = new TGraph();
  TGraph *per_emt_g = new TGraph();
  
  //TCanvas* c0 = new TCanvas("c0", "", 600, 600);
  //TCanvas* c1 = new TCanvas("c1", "", 600, 600);

  TH2D* hist = new TH2D("hist", "; x; y; z", 600, 0, 600, 600, 0, 600);
  
  TF2* gauss = new TF2("gauss", "[0]*( exp(- pow((x - [1]), 2) / (2 * pow([2], 2)) + (- pow((y - [3]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [5]), 2) / (2 * pow([2], 2)) + (- pow((y - [6]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [7]), 2) / (2 * pow([2], 2)) + (- pow((y - [8]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [9]), 2) / (2 * pow([2], 2)) + (- pow((y - [10]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [11]), 2) / (2 * pow([2], 2)) + (- pow((y - [12]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [13]), 2) / (2 * pow([2], 2)) + (- pow((y - [14]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [15]), 2) / (2 * pow([2], 2)) + (- pow((y - [16]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [17]), 2) / (2 * pow([2], 2)) + (- pow((y - [18]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [19]), 2) / (2 * pow([2], 2)) + (- pow((y - [20]), 2) / (2 * pow([4], 2)) )))", -0.5, 0.5, -0.5, 0.5);
  
  //4 points
  //TF2* gauss = new TF2("gauss", "[0]*( exp(- pow((x - [5]), 2) / (2 * pow([2], 2)) + (- pow((y - [6]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [7]), 2) / (2 * pow([2], 2)) + (- pow((y - [8]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [9]), 2) / (2 * pow([2], 2)) + (- pow((y - [10]), 2) / (2 * pow([4], 2)) )) + exp(- pow((x - [11]), 2) / (2 * pow([2], 2)) + (- pow((y - [12]), 2) / (2 * pow([4], 2)) )) )", -1., 1., -1., 1.);

  //TF2* gauss = new TF2("gauss", "[0]*( exp(- pow((x - [1]), 2) / (2 * pow([2], 2)) + (- pow((y - [3]), 2) / (2 * pow([4], 2)) )) )", -1., 1., -1., 1.);

  int iter = 0;
  
  for(double position = 0.2; position < 0.8; position = position + 0.02){
    
    gauss->FixParameter(0, 1/(2*TMath::Pi()*sigma_x*sigma_y*num));
    gauss->FixParameter(1, 0.0);
    gauss->FixParameter(2, sigma_x);
    gauss->FixParameter(3, -0.0);
    gauss->FixParameter(4, sigma_y);
    
    gauss->FixParameter(5, position);
    gauss->FixParameter(6, position);
    
    gauss->FixParameter(7, -position);
    gauss->FixParameter(8, position);
    
    gauss->FixParameter(9, position);
    gauss->FixParameter(10, -position);
    
    gauss->FixParameter(11, -position);
    gauss->FixParameter(12, -position);
    
    gauss->FixParameter(13, 0);
    gauss->FixParameter(14, -position);
    
    gauss->FixParameter(15, 0);
    gauss->FixParameter(16, position);
    
    gauss->FixParameter(17, -position);
    gauss->FixParameter(18, 0);
    
    gauss->FixParameter(19, position);
    gauss->FixParameter(20, 0);
    
    //gauss->GetHistogram()->GetXaxis()->SetRangeUser(0.0, 0.5);
    //gauss->GetHistogram()->GetZaxis()->SetRangeUser(0.0, 2.0);
    gauss->SetRange(-0.5, -0.5, 0.0, 0.5, 0.5, 2.0);
    gauss->SetMinimum(0.0);
    gauss->SetMaximum(1.2);
    //gauss->GetXaxis()->SetTitle("x [cm]");
    //gauss->GetXaxis()->SetTitle("y [cm]");
    
    c0->cd();
    //gauss->Draw("lego2");
    //c0->Update();
    //c1->cd();
    //gauss->Draw("col");
    //c1->Update();
    
    min = gauss->GetMinimum();
    max = gauss->GetMaximum();
    //cout << "min = " << min << " max = " << max << " percentage " << (max - min)/ max * 100 << "\n";
    
    double integral_si;
    double integral_emt = 0;
    integral_si = gauss->Integral(-0.5, 0.5, -0.5, 0.5);
    
    for(double x = -1.0; x < 1.0; x = x + 0.1){
      for(double y = -1.0; y < 1.0; y = y + 0.1){
	if(x * x + y * y <= 0.16){
	  integral_emt = integral_emt + gauss->Eval(x, y) * 0.01;
	}
      }
    }

    cout << "sigma " << sigma_x << "\n";
    cout << "position " << position << "\n";
    cout << "max " << gauss->GetMaximum() << "\n";
    cout << "integral si " << integral_si << "\n";
    cout << "integral EMT " << integral_emt << "\n";
    double average = (max + integral_si) / 2;
    double emt_calc = integral_emt/(0.4 * 0.4 * 3.14);
    double average_emt = (max + emt_calc ) / 2;
    cout << "Si percentage " << (average - integral_si) / average * 100 << "\n";
    double per_si = (average - integral_si) / average * 100;
    cout << "emt_calc " << emt_calc << "\n";
    cout << "EMT percentage " << (average_emt - emt_calc) / average_emt * 100 << "\n\n";
    double per_emt = (average_emt - emt_calc) / average_emt * 100;

    eff_g -> SetPoint(iter, position, integral_si);
    per_si_g -> SetPoint(iter, position, per_si);
    per_emt_g -> SetPoint(iter, position, per_emt);
    iter++;
  } //loop over position

  TCanvas *c0 = new TCanvas("c0", "c0", 1000, 600);
  c0 -> Divide(1, 2);
  c0 -> cd(1);
  per_si_g -> SetMarkerStyle(3);
  per_si_g -> GetXaxis() -> SetTitle("position [cm]");
  per_si_g -> GetYaxis() -> SetTitle("Si variation");
  per_si_g -> Draw("AP");
  c0 -> cd(2);
  per_emt_g -> SetMarkerStyle(3);
  per_emt_g -> GetXaxis() -> SetTitle("position [cm]");
  per_emt_g -> GetYaxis() -> SetTitle("EMT variation");
  per_emt_g -> Draw("AP");

  TCanvas *c1 = new TCanvas("c1", "c1", 1000, 600);
  eff_g -> SetMarkerStyle(3);
  eff_g -> GetXaxis() -> SetTitle("position [cm]");
  eff_g -> GetYaxis() -> SetTitle("Si efficiency");
  eff_g -> Draw("AP");
}
