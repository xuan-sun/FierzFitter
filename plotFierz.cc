#include "comparehist.hh"

#define		HIST_IMAGE_PRINTOUT_NAME	"b_0_notwiddles_statistical error"
#define		INPUT_DATA_FILE			"AnalyzedTextFiles/b_0_SimProcessed_noTwiddles_100keV-650keV_firstPass_forStatError.txt"
#define		INPUT_PARAM_FILE		"/mnt/Data/xuansun/analyzed_files/matchingParams_2010_0.dat"

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvaxIndex, TH1D *hPlot, TString title, TString command);
void FillArrays(TString fileName, TString paramFile, TH1D *hist);

struct event
{
  // from data file
  double b;
  double avg_mOverE;
  double calcErr;
  int binMin;
  int binMax;
  int entriesFitted;
  char rootFile[60];
  int fileNum;
  double chisquared;
  int ndf;
  double chisquaredperdof;
  // from inital param file
  int index;
  double ea;
  double eb;
  double ec;
  double ed;
  double wa;
  double wb;
  double wc;
  double wd;
};


int main()
{
  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH1D *h1 = new TH1D("myhist", "myhist", 50, -0.05, 0.05);

  FillArrays(INPUT_DATA_FILE, INPUT_PARAM_FILE, h1);

  gStyle->SetOptStat("mr");
  gStyle->SetOptTitle(0);

  PlotHist(C, 1, 1, h1, "Extracted b values.", "");

  //prints the canvas with a dynamic TString name of the name of the file
  C -> Print(Form("%s.png", HIST_IMAGE_PRINTOUT_NAME));
  cout << "-------------- End of Program ---------------" << endl;
  plot_program.Run();

  return 0;
}
void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Extracted b");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("N(b)");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
  }

  hPlot -> Draw(command);
}

void FillArrays(TString fileName, TString paramFile, TH1D* hist)
{

  event evt;
  int counter = 0;

  //opens the file that I name in DATA_FILE_IN
  string buf1;
  ifstream infile1;
  cout << "The file being opened is: " << fileName << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName << endl;
/*
  string buf2;
  ifstream infile2;
  cout << "The file being opened is: " << paramFile << endl;
  infile2.open(paramFile);
  if(!infile2.is_open())
    cout << "Problem opening " << paramFile << endl;
*/
  while(true)
  {
    getline(infile1, buf1);
    istringstream bufstream1(buf1);

/*    getline(infile2, buf2);
    istringstream bufstream2(buf2);
*/
    if(!infile1.eof() /*&& !infile2.eof()*/)
    {
      bufstream1 >> evt.b
		>> evt.avg_mOverE
		>> evt.calcErr
		>> evt.binMin
		>> evt.binMax
		>> evt.entriesFitted
		>> evt.rootFile
		>> evt.fileNum
                >> evt.chisquared
		>> evt.ndf
		>> evt.chisquaredperdof;
/*      bufstream2 >> evt.index
		>> evt.ea
		>> evt.eb
		>> evt.ec
		>> evt.ed
		>> evt.wa
		>> evt.wb
		>> evt.wc
		>> evt.wd;
*/
//      if(evt.ed == 0)
      {
	counter++;
        hist -> Fill(evt.b);
      }
    }

    if(infile1.eof() == true /*|| infile2.eof() == true*/)
    {
      break;
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}
