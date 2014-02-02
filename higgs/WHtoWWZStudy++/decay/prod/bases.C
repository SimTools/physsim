TFile *file;

Int_t bases()
{ 
  // Macro example for bases calculation.

  gROOT->Reset();

  // Allocate a file to output bases result.
  file = new TFile("bases.root","RECREATE");

  jsf = new JSFSteer();      // Parameters are obtained from jsf.conf

  Double_t mx   = 150.;
  Double_t fhwz =   1.;
  Double_t fhwa =   0.;
  bs  = new XBoson(mx,fhwz,fhwa);
}
