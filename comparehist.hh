#ifndef	COMPAREHIST_HH
#define	COMPAREHIST_HH
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
#include 	 <TChain.h>
#include	 <TObjArray.h>
#include	 <TFractionFitter.h>
#include	 <TLatex.h>
#include 	 <TMatrixD.h>
using            namespace std;

// forward declarations of functions used
void PlotHist(TCanvas *C, int styleIndex, int canvaxIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors* gPlot, TString title, TString command);
TChain* MakeTChain(TString baseName, TString treeName, int fileNumMin, int fileNumMax);
TH1D* ExtractHistFromChain(TString varName, TString cutsUsed, TChain* chain,
			   TString name, TString title, int nbBins, double minX, double maxX);
double CalculateChiSquared(TH1D* hdat, TH1D* hthBeta, TH1D* hthFierz, double frac0, double frac1, double xBinMin, double xBinMax);

// global TF1's to be accessed by the structs in order to have proper scope
// This is poorly coded and needs to be here for fsum() to access it.
/*TF1 *f1, *f2;

double fsum(double *x, double *par)
{
    return par[0]*(f1->EvalPar(x,par) + par[1]*f2->EvalPar(x,par));
}

struct betaSpectralIndexSpectrum
{
  betaSpectralIndexSpectrum(int si): SI(si) {}

  double operator() (double *x, double *p)  const
  {
    double betaSpectrumValue = 0;

    return betaSpectrumValue;
  }

  int SI;
};
*/
#endif


