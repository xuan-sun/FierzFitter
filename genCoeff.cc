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
#define		INPUT_EQ2ETRUE_PARAMS	"/home/xuansun/Documents/ParallelAnalyzer/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat"

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString command);
void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command);
vector < vector < vector <double> > > GetEQ2EtrueParams(string geometry);
double CalculateErecon(double totalEvis, vector < vector < vector <double> > > tempEQ2Etrue, int type, int side);
TF1* ErrorEnvelope_2010();


TApplication plot_program("FADC_readin",0,0,0,0);

struct EreconFunction
{
  EreconFunction(vector < vector < vector <double> > > calibrationCoeffs, int type, int side, TF1 *f_var):
  COEFFS(calibrationCoeffs), TYPE(type), SIDE(side), XAXIS(f_var) {}

  double operator() (double *x, double *p) const
  {
    return CalculateErecon(XAXIS->EvalPar(x,p), COEFFS, TYPE, SIDE);
//    return CalculateErecon(*x, COEFFS, TYPE, SIDE);
  }

  vector < vector < vector <double> > > COEFFS;
  int TYPE;
  int SIDE;
  TF1* XAXIS;
};

struct TwiddleFunctionErecon
{
  TwiddleFunctionErecon(TF1* base, TF1* varied): BASE(base), VARIED(varied) {}

  double operator() (double *x, double *p) const
  {
    return VARIED -> EvalPar(x, p) - BASE -> EvalPar(x, p);
  }

  TF1* BASE;
  TF1* VARIED;
};

int main(int argc, char *argv[])
{
/*  if(argc != 2)
  {
    cout << "Improper format. Needs: (executable) (either 'save' or 'plot)" << endl;
    return 0;
  }
  // Takes in initial argument and ensures the seed is different for randomizing in ROOT.
  string option = argv[1];
*/  gRandom->SetSeed(0);		// Makes sure that each call to GetRandom() is different

  // Generate twiddle polynomials in EQ space.
  vector <TF1*> twiddles_East;
  vector <TF1*> twiddles_West;

  TF1* pure_East = new TF1("pureE", "x", 0, 800);
  TF1* pure_West = new TF1("pureW", "x", 0, 800);

//  if(option == "save")
//  {
    ofstream outfile;
    outfile.open(PARAM_FILE_NAME, ios::app);

    int counter = 0;

    for(double a = -1; a <= 1; a = a + 0.5)
    {
      for(double b = -0.1; b <= 0.1; b = b + 0.05)
      {
        for(double c = -1e-5; c <= 1e-5; c = c + 5e-6)
        {
          for(double d = -1e-7; d <= 1e-7; d = d + 5e-8)
          {
/*	    outfile 	<< counter << "\t"
			<< a << "\t"
			<< b << "\t"
			<< c << "\t"
			<< d << "\t"
			<< a << "\t"
			<< b << "\t"
			<< c << "\t"
			<< d << "\n";
*/
	    twiddles_East.push_back(new TF1(Form("polyE_%i", counter), "[0] + (1+[1])*x + [2]*x*x + [3]*x*x*x", 0, 800));
	    twiddles_East.back() -> SetParameter(0, a);
	    twiddles_East.back() -> SetParameter(1, b);
	    twiddles_East.back() -> SetParameter(2, c);
	    twiddles_East.back() -> SetParameter(3, d);

            twiddles_West.push_back(new TF1(Form("polyW_%i", counter), "[0] + (1+[1])*x + [2]*x*x + [3]*x*x*x", 0, 800));
            twiddles_West.back() -> SetParameter(0, a);
            twiddles_West.back() -> SetParameter(1, b);
            twiddles_West.back() -> SetParameter(2, c);
            twiddles_West.back() -> SetParameter(3, d);

            counter++;
          }
        }
      }
    }
    outfile.close();
//  }

  cout << "\nNumber of twiddle coefficients generated: " << counter << "\n" << endl;

  // Load the converter to get Erecon from a single EQ value.
  vector < vector < vector <double> > > EQ2Etrue = GetEQ2EtrueParams("2010");

  // Calculate Erecon
  int typeIndex = 0;	// Since we are doing calibration curves on East/West scintillators,
			// it is equivalent to looking at a spectrum of type 0's.
  int sideIndex = 0;	// testing just the Erecons due to East side calibration

  TF1 *Erecon0 = new TF1("Erecon0", EreconFunction(EQ2Etrue, typeIndex, sideIndex, pure_East), 0.1, 800, 0);
  vector <TF1*> Erecon_Twiddles_East;
  vector <TF1*> delta_Erecon_East;
  for(int i = 0; i < counter; i++)
  {
    Erecon_Twiddles_East.push_back(new TF1(Form("Erecon_twiddles_%i", i),
					EreconFunction(EQ2Etrue, typeIndex, sideIndex, twiddles_East[i]), 0.1, 800, 0));
  }


  // Create arrays of Erecon0 and EreconError so we can scatter plot them (hence have error as function of Erecon0).
  vector <double> Evis_axis;
  vector <double> Erecon0_values;
  vector < vector <double> > delta_Erecon_values;	// first index ranges over twiddles
							// second index ranges over step values in Evis space
  double Evis_min = 1;		// both of these values are KeV.
  double Evis_max = 800;
  double Evis_step = 1;
  int nbPoints = 0;
  for(int i = Evis_min; i <= Evis_max; i = i + Evis_step)
  {
    Evis_axis.push_back(i);	// note: i is in whatever units Evis is in.
    Erecon0_values.push_back(Erecon0 -> Eval(i));
    nbPoints++;
  }

  // Get our error envelope for 2010 defined by Michael Mendenhall's thesis.
  TF1* errEnv = ErrorEnvelope_2010();

  // Plot all the twiddle functions and error envelope
  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT->SetStyle("Plain");

/*
  for(int i = 0; i < counter; i++)
  {
    delta_Erecon_East.push_back(new TF1(Form("Delta_Erecon_%i", i),
					TwiddleFunctionErecon(Erecon0, Erecon_Twiddles_East[i]), 0.1, 800, 0));
    vector <double> temp;
    for(int j = Evis_min; j <= Evis_max; j = j + Evis_step)
    {
      temp.push_back(delta_Erecon_East[i] -> Eval(j));
    }

    delta_Erecon_values.push_back(temp);
  }

  vector <TGraph*> graphs;
  for(int i = 0; i < counter; i++)
  {
    graphs.push_back(new TGraph(nbPoints, &(Erecon0_values[0]), &(delta_Erecon_values[i][0])));
  }

  for(int i = 0; i < 625; i++)
  {
    if(i == 0)
	PlotGraph(C, i, 1, graphs[i], "AL");
    else
	PlotGraph(C, i, 1, graphs[i], "SAME");

//    PlotFunc(C, 4, 1, delta_Erecon_East[i], "SAME");
  }

  // Save our plot and print it out as a pdf.
  C -> Print("output_genCoeff.pdf");
*/
  plot_program.Run();
  cout << "-------------- End of Program ---------------" << endl;
  return 0;
}

void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command)
{
  C -> cd(canvasIndex);

  fPlot->SetLineColor(styleIndex % 50);	// only 50 colors in set line color.
  fPlot->GetYaxis()->SetRangeUser(-40, 40);
  fPlot->GetYaxis()->SetTitle("Erecon error");
  fPlot->GetXaxis()->SetTitle("Evis");

  fPlot->Draw(command);
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString command)
{
  C -> cd(canvasIndex);
  gPlot->SetLineColor(styleIndex);
  gPlot->GetYaxis()->SetRangeUser(-40, 40);
  gPlot->Draw(command);
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

//Get the conversion from EQ2Etrue
vector < vector < vector <double> > > GetEQ2EtrueParams(string geometry)
{
  ifstream infile;
  if (geometry=="2010") infile.open(INPUT_EQ2ETRUE_PARAMS);
//  else if (geometry=="2011/2012") infile.open("../simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat");
//  else if (geometry=="2012/2013") infile.open("../simulation_comparison/EQ2EtrueConversion/2012-2013_EQ2EtrueFitParams.dat");
  else {
    cout << "Bad geometry passed to getEQ2EtrueParams\n";
    exit(0);
  }
  vector < vector < vector < double > > > params;
  params.resize(2,vector < vector < double > > (3, vector < double > (6,0.)));

  char holdType[10];
  int side=0, type=0;
  while (infile >> holdType >> params[side][type][0] >> params[side][type][1] >> params[side][type][2] >> params[side][type][3] >> params[side][type][4] >> params[side][type][5])
  {
    std::cout << holdType << " " << params[side][type][0] << " " << params[side][type][1] << " " << params[side][type][2] << " " << params[side][type][3] << " " << params[side][type][4] << " " << params[side][type][5] << std::endl;
    type+=1;
    if (type==3) {type=0; side=1;}
  }
  return params;
}

double CalculateErecon(double totalEvis, vector < vector < vector <double> > > tempEQ2Etrue, int type, int side)
{
  return tempEQ2Etrue[side][type][0]
	+tempEQ2Etrue[side][type][1]*totalEvis
	+tempEQ2Etrue[side][type][2]/(totalEvis+tempEQ2Etrue[side][type][3])
	+tempEQ2Etrue[side][type][4]/((totalEvis+tempEQ2Etrue[side][type][5])*(totalEvis+tempEQ2Etrue[side][type][5]));;
}

TF1* ErrorEnvelope_2010()
{
  TF1* fEnv = new TF1("2010_error_envelope", "(x <= 200)*2.5 + (x > 200 && x <= 500)*(2.5 + 0.0125*(x-200)) + (x>500)*6.25", 0, 1000);

  return fEnv;
}
