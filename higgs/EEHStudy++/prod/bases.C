TFile *file;

Int_t bases()
{ 
  // Macro example for bases calculation.

  gROOT->Reset();

  // Allocate a file to output bases result.
  file = new TFile("bases.root","RECREATE");

  jsf = new JSFSteer();      // Parameters are obtained from jsf.conf

  bs  = new EEHBases();
  bs->Bases();

  file->cd();
  bs->Bh_plot();
  printf(" Scalls is %d\n",bs->GetScalls());
  bs->Userout();
  bs->Write();
  file->Write(); 
}
