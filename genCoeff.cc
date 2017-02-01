// Standard C and C++ libraries
#include         <vector>
#include         <iostream>
#include         <algorithm>
#include         <functional>
#include         <fstream>
#include         <sstream>
#include         <stdlib.h>
#include         <math.h>
#include         <string.h>
#include         <time.h>
#include         <assert.h>

// Pretty much all the ROOT libraries I have ever used.
#include         <TROOT.h>
#include         <TSystem.h>
#include         <TMath.h>
#include         <TF1.h>
#include         <TGaxis.h>
#include         <TGraph.h>
#include         <TGraphErrors.h>
#include         <TCanvas.h>
#include         <TApplication.h>
#include         <TH1.h>
#include         <TProfile.h>
#include         <TObjArray.h>
#include         <TStyle.h>
#include         <TMarker.h>
#include         <TPaveStats.h>
#include         <TPaveText.h>
#include         <TFile.h>
#include         <TLegend.h>
#include         <TLegendEntry.h>
#include         <TH2F.h>
#include         <TRandom.h>
#include         <TTree.h>
#include         <TChain.h>
#include         <TObjArray.h>
#include         <TFractionFitter.h>
#include         <TLatex.h>
#include         <TMatrixD.h>
#include	 <TRandom3.h>

using            namespace std;

// Fundamental constants that get used
const double m_e = 511.00;                                              ///< electron mass, keV/c^2

#define		PARAM_FILE_NAME		"params_2010.txt"

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
double SampleGaus(double mean, double sigma);

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    cout << "Improper format. Needs: (executable) (either 'save' or 'plot)" << endl;
    return 0;
  }
  string option = argv[1];

  gRandom->SetSeed(0);

  double Ce139 = 131.9; 	// KeV central value. +/- 0.1
  double Sn113 = 368.1; 	// KeV. +/- 0.1
  double Bi207_low = 502.5; 	// KeV +/- 0.3
  double Bi207_high = 994.8; 	// KeV +/- 0.5

  double Ce_offset = 0.6;
  double Ce_error = 1.3;	// error taken from error envelope at Ce139 point.
  double Sn_offset = -2.1;	// All values taken from page 99, Michael Mendenhall thesis.
  double Sn_error = 2.3;
  double Bi2_offset = -0.8;
  double Bi2_error = 3.2;
  double Bi1_offset = 2.4;
  double Bi1_error = 4.2;

  TH1D* hCe = new TH1D("Ce139", "Ce139", 100, 100, 200);

  // Run 4 sets of random samplings of Gaussian functions. Save the doubles returned.
  // Create a TF1 4th order polynomial using those doubles.
  // Also print to file those coefficients.
  vector<double> root1, root2, root3, root4;
  for(int i = 0; i < 1000; i++)
  {
    double randVal = SampleGaus(Ce139 + Ce_offset, Ce_error);
    root1.push_back(randVal);
//    hCe->Fill(randVal);
  }

/*  if(option == "save")
  {
    ofstream outfile;
    outfile.open(PARAM_FILE_NAME, ios::app);
    outfile 	<< 1 << "\t"
		<< SampleGaus(Ce139 + Ce_offset, Ce_error) << "\t"
		<< SampleGaus(Sn113 + Sn_offset, Sn_error) << "\t"
		<< SampleGaus(Bi207_low + Bi2_offset, Bi2_error) << "\t"
		<< SampleGaus(Bi207_high + Bi1_offset, Bi1_error) << "\n";

    outfile.close();
  }
  else */if(option == "plot")
  {
  // plot everything and visualize
  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT->SetStyle("Plain");     //on my computer this sets background to white, finally!

  TF1* gausFunction = new TF1("name", "gaus(0)", 0, 850);
  gausFunction -> SetParName(0, "norm");
  gausFunction -> SetParName(1, "mean");
  gausFunction -> SetParName(2, "sigma");
  gausFunction -> SetParameter(0, 100);
  gausFunction -> SetParameter(1, Ce139 + Ce_offset);
  gausFunction -> SetParameter(2, Ce_error);

  for(int i = 0; i < 1; i++)
  {
    hCe->FillRandom("name", 1000);
  }
  PlotHist(C, 1, 1, hCe, "Ce139", "");

  gausFunction->Draw("SAME");

  // prints the canvas with a dynamic TString name of the name of the file
  C -> Print(Form("%s.pdf", "test_genCoeff"));
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  }

  // Create the linear polynomial. Get coefficients by fitting to 4 points we have.

  // Subtract that linear polynomial from the 4th order polynomial.
  // Plot it.




  return 0;
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Counts (N)");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
//    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 3)
  {
    hPlot -> SetFillColor(29);
//    hPlot -> SetFillStyle(3005);
    hPlot -> SetFillStyle(3001);
  }

  hPlot -> Draw(command);
  C -> Update();
}

double SampleGaus(double mean, double sigma)
{
  TF1* gausFunction = new TF1("gaussian", "gaus(0)", 0, 850);
  gausFunction -> SetParName(0, "norm");
  gausFunction -> SetParName(1, "mean");
  gausFunction -> SetParName(2, "sigma");
  gausFunction -> SetParameter(0, 1);
  gausFunction -> SetParameter(1, mean);
  gausFunction -> SetParameter(2, sigma);

  return gausFunction->GetRandom();
}
