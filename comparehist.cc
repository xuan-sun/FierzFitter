#include	"comparehist.hh"

#define		HIST_IMAGE_PRINTOUT_NAME	"Test_comparehist"

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

int main()
{
  // Create a TChain
  TString treeName = Form("SimAnalyzed");
  TChain *MCTheoryChainBeta = MakeTChain("Data/20mill_FierzAndBeta/SimAnalyzed_2010_Beta_paramSet_42", treeName, 0, 20);
  TChain *MCTheoryChainFierz = MakeTChain("Data/20mill_FierzAndBeta/SimAnalyzed_2010_Beta_fierz_paramSet_42", treeName, 0, 20);
  TChain *dataChain = MakeTChain("Data/Sim_b_1/8mill_beta_b_1/SimAnalyzed_2010_Beta_paramSet_42", treeName, 0, 8);

  // Get the Erecon histogram out with appropriate cuts
  TString variableName = Form("Erecon");
  TString cutsUsed = Form("type == 0 && side != 2");
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
  TFractionFitter* fit = new TFractionFitter(dataHist, MCTheory, "V");	// initialise
//  TVirtualFitter* vfit = fit->GetFitter();
  int fitMin = 10;
  int fitMax = 60;
  fit -> SetRangeX(fitMin, fitMax);	// Set range in bin numbers
  int status = fit->Fit();
  if(status != 0)
  {
    cout << "Fit straight up didn't work. Getting out now." << endl;
    return 0;
  }
  TH1D* resultHist = (TH1D*)fit->GetPlot();	// extract the plot from the fit.

  double avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, fitMin, fitMax);
  cout << "Our average value and hence scaling factor is: " << avg_mE << endl;

  // This works and returns the elements of the covariance matrix. Indexed from 0.
//  cout  << "Testing TVirtualFitter. What is the return of GetCovarianceMatrixElement? "
//	<< vfit->GetCovarianceMatrixElement(0, 0) << endl;

  // Get valuable numbers for later
  double chisquared = fit->GetChisquare();
  int ndf = fit->GetNDF();
  double frac0Val, frac0Err, frac1Val, frac1Err;
  fit->GetResult(0, frac0Val, frac0Err);
  fit->GetResult(1, frac1Val, frac1Err);

  double bErr = Fierz_b_Error(frac0Val, frac0Err, frac1Val, frac1Err, avg_mE, -0.989);
  cout << "b error without the m/E term: " << bErr << endl;
  cout << "Also b value without m/E term: " << frac1Val/(frac0Val*avg_mE) << endl;
  cout << "To compare, the limitation from 100KeV and up is: " << 10.1/sqrt(dataHist->GetEntries()) << endl;

  // plot everything and visualize
  TCanvas *C = new TCanvas("canvas", "canvas");
  gROOT->SetStyle("Plain");	//on my computer this sets background to white, finally!
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat("en");
  gStyle->SetStatH(0.45);
  gStyle->SetStatW(0.45);
  PlotHist(C, 1, 1, dataHist, "Test of non-zero Fierz fit.", "");
  PlotHist(C, 2, 1, resultHist, "", "SAME");
//  PlotHist(C, 1, 1, mcTheoryHistBeta, "", "");
//  PlotHist(C, 2, 1, mcTheoryHistFierz, "", "SAME");


  // update the legend to include valuable variables.
  TPaveStats *ps = (TPaveStats*)C->GetPrimitive("stats");
  ps->SetName("mystats");
  TList *listOfLines = ps->GetListOfLines();
  TLatex *myText1 = new TLatex(0,0,Form("#Chi^{2} = %f", chisquared));
  listOfLines->Add(myText1);
  TLatex *myText2 = new TLatex(0,0,Form("NDF = %d", ndf));
  listOfLines->Add(myText2);
  TLatex *myText3 = new TLatex(0,0,Form("#frac{#Chi^{2}}{NDF} = %f", chisquared/ndf));
  listOfLines->Add(myText3);
  TLatex *myText4 = new TLatex(0,0,Form("Beta = %f #pm %f", frac0Val, frac0Err));
  listOfLines->Add(myText4);
  TLatex *myText5 = new TLatex(0,0,Form("Fierz = %f #pm %f", frac1Val, frac1Err));
  listOfLines->Add(myText5);
  // the following line is needed to avoid that the automatic redrawing of stats
  dataHist->SetStats(0);

  // prints the canvas with a dynamic TString name of the name of the file
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

TChain* MakeTChain(TString baseName, TString treeName, int fileNumMin, int fileNumMax)
{
  TChain* chain = new TChain(treeName.Data());

  for(int i = fileNumMin; i < fileNumMax; i++)
  {
    chain -> AddFile(Form("%s_%i.root", baseName.Data(), i));
  }

  cout << "Loaded trees from files identified by the template: " << baseName.Data() << "_#.root" << endl;

  return chain;
}

double CalculateChiSquared(TH1D* hdat, TH1D* hthBeta, TH1D* hthFierz, double frac0, double frac1, double xBinMin, double xBinMax)
{
  TH1D* hmc = new TH1D("MCHist", "MCHist", 100, 0, 1000);
  hmc -> Add(hthBeta, hthFierz, frac0, frac1);
  double norm = hdat->GetEntries() / hmc -> GetEntries();
  hmc -> Scale(norm);

  double chisquare = 0;
  for(int i = xBinMin; i <= xBinMax; i++)
  {
    chisquare = chisquare + ((hdat->GetBinContent(i) - hmc->GetBinContent(i))*
			     (hdat->GetBinContent(i) - hmc->GetBinContent(i)))/
			    hdat->GetBinContent(i);
  }

  return chisquare;
}

double CalculateAveragemOverE(TH1D* gammaSM, int binMin, int binMax)
{
  double num = 0;
  double denom = 0;

  for(int i = binMin; i < binMax; i++)
  {
    num = num + (m_e*gammaSM->GetBinContent(i)) / (gammaSM->GetBinCenter(i) + m_e);
    denom = denom + gammaSM->GetBinContent(i);
  }

  return num/denom;
}

double Fierz_b_Error(double f0v, double f0e, double f1v, double f1e, double avgInverseW, double corrCoeff)
{
  double errb = 0;

  errb = (f1v/(f0v*avgInverseW))
	 * sqrt((f0e/f0v)*(f0e/f0v) + (f1e/f1v)*(f1e/f1v) - (2*f0e*f1e*corrCoeff) / (f0v*f1v));

  return errb;
}
