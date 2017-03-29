#include	"comparehist.hh"

#define		INPUT_DATA_FILE			"AnalyzedTextFiles/b_0_SimProcessed_EastWest_allTwiddles_100keV-650keV_secondPass.txt"
#define		INPUT_PARAM_FILE		"/mnt/Data/xuansun/analyzed_files_statDependent/matchingParams_2010_0.dat"

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


void loadtrees()
{
  TFile file("fitter_tree_secondPass.root", "RECREATE");
  TTree* fitTree = new TTree("fitTree", "tree for examining fit results");
  FillArrays(INPUT_DATA_FILE, INPUT_PARAM_FILE, fitTree);
  fitTree->Write();


  int nEntries = fitTree- > GetEntries();

  double chisquared_used = 0;
  int ndf_used = 0;
  double ed_used = 0;
  double wd_used = 0;

  fitTree -> SetBranchAddress("chisquared", &chisquared_used);
  fitTree -> SetBranchAddress("ndf", &ndf_used);
  fitTree -> SetBranchAddress("East_d", &ed_used);
  fitTree -> SetBranchAddress("West_d", &wd_used);

  double totalProbWithCuts = 0;

  for(int i = 0; i < nEntries; i++)
  {
    fitTree->GetEntry(i);
    // These are the cuts we are applying
    if(ed_used == 0 && wd_used == 0)
    {
      totalProbWithCuts = totalProbWithCuts + TMath::Prob(chisquared_used, ndf_used);
    }
  }

  gStyle->SetOptStat("mr");
  gStyle->SetOptTitle(0);

//  fitTree -> Draw("b", Form("(TMath::Prob(chisquared, ndf) / %f)*(East_d == 0 && West_d == 0)", totalProbWithCuts));
  fitTree -> Draw("b", "(1/sqrt(chisquaredperdof)) * (East_d == 0 && West_d == 0)");
}

void FillArrays(TString fileName, TString paramFile, TTree* tree)
{
  event evt;

  // set branches
  tree -> Branch("b", &evt.b, "evt.b/D");
  tree -> Branch("mOverE", &evt.avg_mOverE, "evt.ave_mOverE/D");
  tree -> Branch("bError", &evt.calcErr, "evt.bError/D");
  tree -> Branch("binMin", &evt.binMin, "evt.binMin/I");
  tree -> Branch("binMax", &evt.binMax, "evt.binmax/I");
  tree -> Branch("entries", &evt.entriesFitted, "evt.entriesFitted/I");
  tree -> Branch("paramNumber", &evt.fileNum, "evt.fileNum/I");
  tree -> Branch("chisquared", &evt.chisquared, "evt.chisquared/D");
  tree -> Branch("ndf", &evt.ndf, "evt.ndf/I");
  tree -> Branch("chisquaredperdof", &evt.chisquaredperdof, "evt.chisquaredperdof/D");
  tree -> Branch("paramIndex", &evt.index, "evt.index/I");
  tree -> Branch("East_a", &evt.ea, "evt.ea/D");
  tree -> Branch("East_b", &evt.eb, "evt.eb/D");
  tree -> Branch("East_c", &evt.ec, "evt.ec/D");
  tree -> Branch("East_d", &evt.ed, "evt.ed/D");
  tree -> Branch("West_a", &evt.wa, "evt.wa/D");
  tree -> Branch("West_b", &evt.wb, "evt.wb/D");
  tree -> Branch("West_c", &evt.wc, "evt.wc/D");
  tree -> Branch("West_d", &evt.wd, "evt.wd/D");

  //opens the file that I name in DATA_FILE_IN
  string buf1;
  ifstream infile1;
  cout << "The file being opened is: " << fileName << endl;
  infile1.open(fileName);

  //a check to make sure the file is open
  if(!infile1.is_open())
    cout << "Problem opening " << fileName << endl;

  string buf2;
  ifstream infile2;
  cout << "The file being opened is: " << paramFile << endl;
  infile2.open(paramFile);
  if(!infile2.is_open())
    cout << "Problem opening " << paramFile << endl;

  while(true)
  {
    getline(infile1, buf1);
    istringstream bufstream1(buf1);

    getline(infile2, buf2);
    istringstream bufstream2(buf2);

    if(!infile1.eof() && !infile2.eof())
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
      bufstream2 >> evt.index
		>> evt.ea
		>> evt.eb
		>> evt.ec
		>> evt.ed
		>> evt.wa
		>> evt.wb
		>> evt.wc
		>> evt.wd;

    }

    if(infile1.eof() == true || infile2.eof() == true)
    {
      break;
    }

    tree->Fill();
  }

  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;

}
