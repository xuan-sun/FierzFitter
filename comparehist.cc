#include	"comparehist.hh"

#define		HIST_IMAGE_PRINTOUT_NAME	"Test_comparehist"

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

int main()
{
  // Create a TChain
  TString treeName = Form("SimAnalyzed");
  TChain *MCTheoryChainBeta = MakeTChain("Data/SimBeta/SimAnalyzed_2010_Beta_b_0_paramSet_42", treeName, 1);
  TChain *MCTheoryChainFierz = MakeTChain("Data/SimFierz/SimAnalyzed_2010_Beta_b_inf_paramSet_42", treeName, 1);
  TChain *dataChain = MakeTChain("Data/Simb_1/SimAnalyzed_2010_Beta_FierzIs1_paramSet_42", treeName, 1);

  // Get the Erecon histogram out with appropriate cuts
  TString variableName = Form("Erecon");
  TString cutsUsed = Form("type != 4 && side !=2");
  TH1D* dataHist = ExtractHistFromChain(variableName, cutsUsed, dataChain,
				      "myHist", "Test of comparehist code", 100, 0, 1000);
  TH1D* mcTheoryHistBeta = ExtractHistFromChain(variableName, cutsUsed, MCTheoryChainBeta,
                                      "mcBeta", "Test of comparehist code", 100, 0, 1000);
  TH1D* mcTheoryHistFierz = ExtractHistFromChain(variableName, cutsUsed, MCTheoryChainFierz,
                                      "mcFierz", "Test of comparehist code", 100, 0, 1000);

  // Create a TFractionFitter and do the fit.
  TObjArray *MCTheory = new TObjArray(2);
  MCTheory -> Add(mcTheoryHistBeta);
  MCTheory -> Add(mcTheoryHistFierz);
  TFractionFitter* fit = new TFractionFitter(dataHist, MCTheory);	// initialise
  fit -> SetRangeX(5, 65);	// Set range in bin numbers
  int status = fit->Fit();	// perform the fit
  TH1D* resultHist = (TH1D*)fit->GetPlot();	// extract the plot from the fit.

  // Get valuable numbers and add them to the legend
  double chisquared = fit->GetChisquare();
  double ndf = fit->GetNDF();
  double frac0Val, frac0Err, frac1Val, frac1Err;
  fit->GetResult(0, frac0Val, frac0Err);
  fit->GetResult(1, frac1Val, frac1Err);

  // Add everything to legend.


  // plot everything and visualize
  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT->SetStyle("Plain");	//on my computer this sets background to white, finally!
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("em");
  PlotHist(C, 1, 1, dataHist, "Test of non-zero Fierz fit.", "");
  PlotHist(C, 1, 2, resultHist, "", "SAME");

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
  hPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Counts (N)");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(0, 0.000004);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
//    hPlot -> SetFillStyle(3004);
    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
//    hPlot -> SetFillStyle(3005);
    hPlot -> SetFillStyle(3001);
  }
  if(styleIndex == 3)
  {
    hPlot -> SetFillColor(29);
//    hPlot -> SetFillStyle(3005);
    hPlot -> SetFillStyle(3001);
  }

  hPlot -> Draw(command);
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors* gPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  gPlot -> SetTitle(title);
  gPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  gPlot -> GetXaxis() -> SetRangeUser(0, 1000);
  gPlot -> GetXaxis() -> CenterTitle();
  gPlot -> GetYaxis() -> SetTitle("Hits (N)");
  gPlot -> GetYaxis() -> CenterTitle();

  if(styleIndex == 1)
  {
    gPlot -> SetMarkerColor(46);
  }
  if(styleIndex == 2)
  {
    gPlot -> SetMarkerColor(38);
  }
  if(styleIndex == 3)
  {
    gPlot -> SetMarkerColor(30);
  }
  gPlot -> SetMarkerStyle(7);
  gPlot -> SetMarkerSize(1);

  // add a flat line at y = 0
  TLine *line = new TLine(0, 0, 1000, 0);

  gPlot -> Draw(command);
  line -> Draw("SAME");

}

TH1D* ExtractHistFromChain(TString varName, TString cutsUsed, TChain* chain,
                           TString name, TString title, int nbBins, double minX, double maxX)
{
  TH1D* hist = new TH1D(name.Data(), title.Data(), nbBins, minX, maxX);

  chain -> Draw(Form("%s >> %s",varName.Data(),name.Data()), cutsUsed.Data());

  return hist;
}

TChain* MakeTChain(TString baseName, TString treeName, int nbFiles)
{
  TChain* chain = new TChain(treeName.Data());

  for(int i = 0; i < nbFiles; i++)
  {
    chain -> AddFile(Form("%s_%i.root", baseName.Data(), i));
  }

  cout << "Loaded trees from files identified by the template: " << baseName.Data() << "_#.root" << endl;

  return chain;
}