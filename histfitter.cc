#include	"comparehist.hh"

#define		HIST_IMAGE_PRINTOUT_NAME	"Test_master_histfitter"
#define		OUTPUT_ANALYSIS_FILE		"Fierz_Analysis_b_0.txt"

int main()
{
  TString treeName = Form("Evts");
  TChain *MCTheoryChainBeta = MakeTChain("/home/xuansun/Documents/Analysis_Code/ucna_g4_2.1/UCN/UK_EventGen_2016/Evts_Files/b_0_200mill/Evts", treeName, 0, 100);
  TChain *MCTheoryChainFierz = MakeTChain("/home/xuansun/Documents/Analysis_Code/ucna_g4_2.1/UCN/UK_EventGen_2016/Evts_Files/b_inf_100mill/Evts", treeName, 0, 100);
  TChain *dataChain = MakeTChain("/home/xuansun/Documents/Analysis_Code/ucna_g4_2.1/UCN/UK_EventGen_2016/Evts_Files/b_0_200mill/Evts", treeName
				, 299, 300);

  TString variableName = Form("KE");
  TString cutsUsed = Form("");
  TH1D* mcTheoryHistBeta = ExtractHistFromChain(variableName, cutsUsed, MCTheoryChainBeta,
                                      "mcBeta", "Beta", 100, 0, 1000);
  TH1D* mcTheoryHistFierz = ExtractHistFromChain(variableName, cutsUsed, MCTheoryChainFierz,
                                      "mcFierz", "Fierz", 100, 0, 1000);
  TH1D* dataHist = ExtractHistFromChain(variableName, cutsUsed, dataChain, "myHist", "Data", 100, 0, 1000);

  TObjArray *MCTheory = new TObjArray(2);
  MCTheory -> Add(mcTheoryHistBeta);
  MCTheory -> Add(mcTheoryHistFierz);
  TFractionFitter* fit = new TFractionFitter(dataHist, MCTheory, "Q");  // initialise

  int fitMin = 10;
  int fitMax = 85;
  fit -> SetRangeX(fitMin, fitMax);
  int status = fit->Fit();
  if(status != 0)
  {
    cout << "Fit straight up didn't work. Getting out now." << endl;
    return 0;
  }

  TH1D* resultHist = (TH1D*)fit->GetPlot();	// extract the plot from the fit.
  int entries = 0;
  for(int i = fitMin; i < fitMax; i++)
  {
    entries = entries + resultHist->GetBinContent(i);
  }

  double avg_mE = CalculateAveragemOverE(mcTheoryHistBeta, fitMin, fitMax);

  // Get valuable numbers for later
  double chisquared = fit->GetChisquare();
  int ndf = fit->GetNDF();
  double frac0Val, frac0Err, frac1Val, frac1Err;
  fit->GetResult(0, frac0Val, frac0Err);
  fit->GetResult(1, frac1Val, frac1Err);

  ofstream outfile;
  outfile.open(OUTPUT_ANALYSIS_FILE, ios::app);
  outfile << frac1Val/(frac0Val*avg_mE) << "\t"
	  << 10.1/sqrt(dataHist->GetEntries()) << "\t"
          << fitMin << "\t" << fitMax << "\t"
	  << entries << "\t"
          << "Evts_" << 299 << ".root" << "\t"
	  << chisquared << "\t"
	  << ndf << "\t"
	  << chisquared/ndf << "\n";
  outfile.close();

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

  cout << "Completed storing histogram of variable " << varName.Data() << " into histogram " << name.Data() << endl;

  return hist;
}

TChain* MakeTChain(TString baseName, TString treeName, int fileNumMin, int fileNumMax)
{
  TChain* chain = new TChain(treeName.Data());

  for(int i = fileNumMin; i < fileNumMax; i++)
  {
    chain -> AddFile(Form("%s_%i.root", baseName.Data(), i));
  }

  cout << "Loaded trees from files identified by the template: " << baseName.Data() << "_#.root \n"
       << "Lower index = " << fileNumMin << " and upper index = " << fileNumMax << endl;

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

double Fierz_b_Error(double f0v, double f0e, double f1v, double f1e, double avgInverseW
		    , double cov00, double cov01, double cov10, double cov11)
{
  double errb = 0;

  errb = (f1v/(f0v*avgInverseW))
	 * sqrt((f0e/f0v)*(f0e/f0v) + (f1e/f1v)*(f1e/f1v) - (2*cov01*avgInverseW*f0v*avgInverseW*f1v));


  return errb;
}
