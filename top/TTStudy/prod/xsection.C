void xsection(Double_t rsi, Char_t *ofile = "bases.root")
{
// Macro example for bases calculation.
  gROOT->Reset();
  TFile file(ofile,"RECREATE");

  jsf = new JSFSteer();    // required to read parameter from jsf.conf
  
  bases = new TTBases();
  bases->SetRoots(rsi); // rsi [GeV]

//  bases->SetNoOfSample(10000);
//  bases->SetIteration1( 0.2, 10);
//  bases->SetIteration2( 0.1, 10);
  bases->Bases();
  bases->Bh_plot();
  bases->Userout();
  bases->Write();
#if 1
  Double_t rs  = bases->GetRoots();
  Double_t sg  = bases->GetEstimate();
  Double_t dsg = bases->GetError();
  ofstream outf("xsection.dat",std::ios::app);
  outf << rs << " " << sg << " " << dsg << endl;
#endif

  file.Write();
}
