plotpolynomials()
{
  FillArrays("params_2010.txt", "SAME");
}

void FillArrays(TString fileName, TString command)
{
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  //a check to make sure the file is open
  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  int counter = 0;

  int paramNb = 0;
  double aE, bE, cE, dE, aW, bW, cW, dW;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);

    if(!bufstream.eof())
    {
      bufstream >> paramNb
                >> aE
                >> bE
                >> cE
                >> dE
                >> aW
                >> bW
                >> cW
                >> dW;

      TF1* polyE = new TF1(Form("polyE_%i", counter), "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0, 800);
      polyE->SetParameter(0, aE);
      polyE->SetParameter(1, bE);
      polyE->SetParameter(2, cE);
      polyE->SetParameter(3, dE);
      polyE->SetLineColor(counter);

      if(counter == 0)
      {
        polyE->Draw();
      }
      else
        polyE->Draw(command);
/*
      TF1* polyW = new TF1(Form("polyW_%i", counter), "[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0, 1000);
      polyW->SetParameter(0, aW);
      polyW->SetParameter(1, bW);
      polyW->SetParameter(2, cW);
      polyW->SetParameter(3, dW);
      polyW->SetLineColor(counter);
      polyW->Draw(command);
*/
    }

//    if(counter == 100) { break; }

    counter++;
  }


  cout << "Data from " << fileName << " has been filled into all arrays successfully." << endl;
}


