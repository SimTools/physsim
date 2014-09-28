TFile *file;

Int_t bases()
{ 
  // Macro example for bases calculation.

  gROOT->Reset();

  // Allocate a file to output bases result.
  file = new TFile("bases.root","RECREATE");

  jsf = new JSFSteer();      // Parameters are obtained from jsf.conf

  bs  = new ENWBases();
  bs->Bases();

  bs->Bh_plot();
  file->cd();
  printf(" Scalls is %d\n",bs->GetScalls());
  bs->Userout();
  bs->Write();
  file->Write(); 
}
