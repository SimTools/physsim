{
// Macro example for bases calculation.

  gROOT->Reset();
  TFile file("bases.root","RECREATE");

  jsf = new JSFSteer();    // required to read parameter from jsf.conf

  bases = new XN1XN2Bases();

//  bases->SetNoOfSample(5000);
//  bases->SetIteration1( 0.2, 1);
//  bases->SetIteration2( 0.1, 1);
  bases->Bases();
  bases->Bh_plot();
  bases->Userout();

  bases->Write();
  file->Write();

}



