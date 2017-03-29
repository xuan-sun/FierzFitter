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

// Input and output names and paths used in the code.
// My pseudo version of environment variables.
#define		PARAM_FILE_NAME		"params_2010_EastWest_2nd_order.txt"
#define		INPUT_EQ2ETRUE_PARAMS	"/home/xuansun/Documents/ParallelAnalyzer/simulation_comparison/EQ2EtrueConversion/2011-2012_EQ2EtrueFitParams.dat"

// Plotting functions.
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraph *gPlot, TString command);
void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command);

// Michael Brown's SimulationProcessor.cpp uses this to convert to Erecon
vector < vector < vector <double> > > GetEQ2EtrueParams(string geometry);

// Does the Erecon calculation using calibration information.
double CalculateErecon(double totalEvis, vector < vector < vector <double> > > tempEQ2Etrue, int type, int side);

// Michael Mendenhall's 2010 error envelope.
TF1* ErrorEnvelope_2010(double factor);

// Perform a single twiddle so we can loop over it in main(), check against a save condition.
// Return whether or not the thrown polynomial passed the save condition.
bool PerformVariation(double a, double b, double c, double d, int numPassed,
                      vector < vector < vector <double> > > EQ2Etrue, TRandom3 *factor,
		      int passNumber, int sideIndex);

// Used for visualization, keeps the graph on screen.
TApplication plot_program("FADC_readin",0,0,0,0);

// Counters to allow us to properly sample beyond 1 sigma distribution.
int num1sigma = 0;
int num2sigma = 0;

// coefficients we are saving on first pass. Speeds up program
vector < vector <double> > goodTwiddles;

// Takes an Evis function and converts it to an Erecon function.
struct EreconFunction
{
  EreconFunction(vector < vector < vector <double> > > calibrationCoeffs, int type, int side, TF1 *f_var):
  COEFFS(calibrationCoeffs), TYPE(type), SIDE(side), XAXIS(f_var) {}

  double operator() (double *x, double *p) const
  {
    return CalculateErecon(XAXIS->EvalPar(x,p), COEFFS, TYPE, SIDE);
  }

  vector < vector < vector <double> > > COEFFS;
  int TYPE;
  int SIDE;
  TF1* XAXIS;
};

// Creates a twiddle function in Erecon space. Mostly this just does a subtraction.
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

// ----------------------------------------------------------------- //
// -------------------- Start of program code. --------------------- //
// ----------------------------------------------------------------- //

int main(int argc, char *argv[])
{
  // Ensures the seed is different for randomizing in ROOT.
  TRandom3* engine = new TRandom3(0);
  gRandom->SetSeed(0);

  // Start the plotting stuff so we can loop and use "SAME" as much as possible.
  TCanvas *C = new TCanvas("canvas", "canvas");
//  C -> Divide(5, 2);
  C -> cd(1);
  gROOT->SetStyle("Plain");
  TF1* line = new TF1("origin", "0", 0, 1000);
  line -> GetYaxis() -> SetRangeUser(-0.25, 0.25);
  line -> GetYaxis() -> SetTitle("E_{recon} Error Fraction");
  line -> GetXaxis() -> SetTitle("E_{recon} (keV)");
  line -> SetTitle("Fractional Spectrum Variations");
  line -> SetLineStyle(1);
  line -> SetLineWidth(1);
  line -> Draw();


  // Load the converter to get Erecon from a single EQ value.
  cout << "Using following calibration for 2010 geometry to convert Evis to Erecon..." << endl;
  vector < vector < vector <double> > > converter = GetEQ2EtrueParams("2010");

  int counter = 0;
  int index = 0;
  double ae, be, ce, de, aw, bw, cw, dw;

  string buf1;
  ifstream infile1;
  cout << "The file being opened is: " << PARAM_FILE_NAME << endl;
  infile1.open(PARAM_FILE_NAME);
  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << PARAM_FILE_NAME << endl;

  while(true)
  {
    getline(infile1, buf1);
    istringstream bufstream1(buf1);

    if(!infile1.eof())
    {
      bufstream1 >> index
                >> ae
                >> be
                >> ce
                >> de
                >> aw
                >> bw
                >> cw
                >> dw;
    }

    if(infile1.eof() == true)
    {
      break;
    }

    // at the level of each line read-in, i.e. each loop, perform twiddle and plot the residual
    if(index == 2)	// this if statement ensures we only plot a fixed number, not all
    {
									// 0 = East side calib, start with that
      bool save = PerformVariation(ae, be, ce, de, index, converter, engine, 2, 0);

      if(save == true)
      {
	counter++;
      }
    }
  }

  int n = 10;
  double x[] = {70, 150, 200, 250, 300, 350, 400, 450, 500, 550};
  double y[] = {-0.23, -0.009, };

  TGraph* mpm_frac = new TGraph()



  // Save our plot and print it out as a pdf.
  C -> Print("output_fracResiduals.pdf");
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}

bool PerformVariation(double a, double b, double c, double d, int numPassed,
		      vector < vector < vector <double> > > EQ2Etrue, TRandom3* factor,
		      int passNumber, int sideIndex)
{
  // booleans needed for first pass to count successes
  bool firstBand = true;
  bool secondBand = true;

  bool saveCondition = true;
  bool throwCondition = false;

  double xMin = 0.1;	// For all polynomial ranges, in Evis units.
  double xMax = 900;

  // Generate twiddle polynomials in EQ space (equivalently, Evis space).
  TF1* pure_Evis = new TF1("pureE", "x", xMin, xMax);
  TF1* twiddle_Evis = new TF1("polyE", "[0] + (1+[1])*x + [2]*x*x + [3]*x*x*x", xMin, xMax);
  twiddle_Evis -> SetParameter(0, a);
  twiddle_Evis -> SetParameter(1, b);
  twiddle_Evis -> SetParameter(2, c);
  twiddle_Evis -> SetParameter(3, d);

  // Calculate Erecon
  int typeIndex = 0;

  // Create our twiddled and untwiddled functions in Erecon space.
  TF1 *Erecon0 = new TF1("Erecon0", EreconFunction(EQ2Etrue, typeIndex, sideIndex, pure_Evis), xMin, xMax, 0);
  TF1* Erecon_Twiddle = new TF1("Erecon_twiddle", EreconFunction(EQ2Etrue, typeIndex, sideIndex, twiddle_Evis), xMin, xMax, 0);
  TF1* delta_Erecon = new TF1("Delta_Erecon", TwiddleFunctionErecon(Erecon0, Erecon_Twiddle), xMin, xMax, 0);

  // Create arrays of Erecon0 and EreconError so we can scatter plot them (hence have error as function of Erecon0).
  vector <double> Evis_axis;
  vector <double> Erecon0_values;
  vector <double> delta_Erecon_values;
  vector <double> frac_delta_Erecon_values;
  double Evis_min = 1;
  double Evis_max = xMax;
  double Evis_step = 1;
  int nbPoints = 0;
  for(int i = Evis_min; i <= Evis_max; i = i + Evis_step)
  {
    Evis_axis.push_back(i);     // note: i is in whatever units Evis is in.
    Erecon0_values.push_back(Erecon0 -> Eval(i));
    delta_Erecon_values.push_back(delta_Erecon -> Eval(i));
    frac_delta_Erecon_values.push_back((delta_Erecon -> Eval(i)) / (Erecon0 -> Eval(i)));
    nbPoints++;
  }

  // Create our scatter plot as a TGraph.
  TGraph* graph = new TGraph(nbPoints, &(Erecon0_values[0]), &(frac_delta_Erecon_values[0]));

  // Get our error envelopes so we can check polynomial values against them.
  TF1* errEnv1 = ErrorEnvelope_2010(1);
  TF1* errEnv2 = ErrorEnvelope_2010(2);

  // Check our polynomial (the scatter plot) against a save condition.
  double x, y;
  for(int i = 1; i <= graph->GetN(); i++)
  {
    graph->GetPoint(i, x, y);

    // if, at any point, we are over 2 sigma away, exit and don't save and don't throw a number.
    if(abs(y) > errEnv2->Eval(x))
    {
      if(passNumber == 1)
      {
        firstBand = false;
        secondBand = false;
        break;
      }
      else if(passNumber == 2)
      {
        saveCondition = false;
        throwCondition = false;
        break;
      }
    }
    // if we are ever between 1 and 2 sigma and never outside 2 sigma, set flag to throw number to true
    else if(abs(y) > errEnv1->Eval(x) && abs(y) < errEnv2->Eval(x))
    {
      if(passNumber == 1)
      {
        firstBand = false;
        secondBand = true;
      }
      if(passNumber == 2)
      {
        throwCondition = true;
      }
    }
  }

  if(passNumber == 1)
  {
    if(firstBand == true)	// i.e. if we complete the search and never trip either flag, we're in first band
    {
      secondBand = false;
    }
    if(firstBand == true)	// now do the counting.
    {
      num1sigma++;
    }
    else if(secondBand == true)
    {
      num2sigma++;
    }
    if(firstBand == true || secondBand == true)
    {
      vector <double> temp;
      temp.push_back(a);
      temp.push_back(b);
      temp.push_back(c);
      temp.push_back(d);
      goodTwiddles.push_back(temp);
    }

  }
  else if(passNumber == 2)
  {
    // if our condition is tripped, means our curve lies (not exclusively from below) between 1 and 2 sigma
    if(throwCondition == true)
    {
      double percentageToSave = (double)num1sigma/num2sigma;
      if(factor->Rndm() < (1 - percentageToSave))
      {
        saveCondition = false;
      }
    }
  }

  if(passNumber == 2)
  {
    if(saveCondition == true)
    {
      // Plotting stuff
      graph->SetLineColor(numPassed % 50);
      graph->Draw("SAMEP");
    }
    // memory management. Delete the left-over pointers. Absolutely necessary or program doesn't run.
    else if(saveCondition == false)
    {
      delete graph;
    }
  }
  else if(passNumber == 1)
  {
    delete graph;
  }
  delete pure_Evis;
  delete twiddle_Evis;
  delete Erecon0;
  delete Erecon_Twiddle;
  delete delta_Erecon;
  delete errEnv2;
  delete errEnv1;

  return saveCondition;
}

void PlotFunc(TCanvas *C, int styleIndex, int canvasIndex, TF1 *fPlot, TString command)
{
  C -> cd(canvasIndex);

  fPlot->SetLineColor(styleIndex % 50); // only 50 colors in set line color.
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
    std::cout << holdType << " " << params[side][type][0] << " " << params[side][type][1]
	      << " " << params[side][type][2] << " " << params[side][type][3]
	      << " " << params[side][type][4] << " " << params[side][type][5] << std::endl;
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

TF1* ErrorEnvelope_2010(double factor)
{
  TF1* fEnv = new TF1("2010_error_envelope", Form("%f*((x <= 200)*2.5 + (x > 200 && x <= 500)*(2.5 + 0.0125*(x-200)) + (x>500)*6.25)", factor), 0, 1000);

  return fEnv;
}
