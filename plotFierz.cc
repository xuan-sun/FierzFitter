#include	"comparehist.hh"

#define		HIST_IMAGE_PRINTOUT_NAME	"test_plotFierz"
#define		INPUT_DATA_FILE			"AnalyzedTextFiles/Fierz_Analysis_b_1_SimAnalyzed_100-650KeV.txt"

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void PlotHist(TCanvas *C, int styleIndex, int canvaxIndex, TH1D *hPlot, TString title, TString command);
void FillArrays(TString fileName, TH1D *hist);

struct event
{
  double b;
  double avg_mOverE;
  double calcErr;
  double frac0Val;
  double frac0Err;
  double frac1Val;
  double frac1Err;
  int binMin;
  int binMax;
  int entriesFitted;
  char rootFile[25];
  double chisquared;
  int ndf;
  double chisquaredperdof;
};


int main()
{
  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT -> SetStyle("Plain");	//on my computer this sets background to white, finally!

  TH1D *h1 = new TH1D("myhist", "myhist", 100, 0.5, 1.5);

  FillArrays(INPUT_DATA_FILE, h1);

  PlotHist(C, 1, 1, h1, "Extracted b values.", "");

  //prints the canvas with a dynamic TString name of the name of the file
  C -> Print(Form("%s.pdf", HIST_IMAGE_PRINTOUT_NAME));
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

void FillArrays(TString fileName, TH1D* hist)
{

  event evt;

  //opens the file that I name in DATA_FILE_IN
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  //a check to make sure the file is open
  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);

    if(!bufstream.eof())
    {
      bufstream >> evt.b
		>> evt.avg_mOverE
		>> evt.calcErr
		>> evt.frac0Val
		>> evt.frac0Err
		>> evt.frac1Val
		>> evt.frac1Err
		>> evt.binMin
		>> evt.binMax
		>> evt.entriesFitted
		>> evt.rootFile
                >> evt.chisquared
		>> evt.ndf
		>> evt.chisquaredperdof;
        hist -> Fill(evt.b);
    }
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}
