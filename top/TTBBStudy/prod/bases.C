{
  gROOT->Reset();
  TFile file("bases.root","RECREATE");

  jsf = new JSFSteer();    // required to read parameter from jsf.conf
  
  bases = new TTBBBases();

  bases->Bases();
  bases->Bh_plot();
  bases->Userout();
  bases->Write();

  file->Write();
}
