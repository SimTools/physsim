{
// Macro example for bases calculation.

  gROOT->Reset();
  TFile file("bases.root","RECREATE");
  
  if( strncmp(gSystem->HostName(),"ccjlc",5)  != 0 ) {
    if( strncmp(gSystem->Getenv("OSTYPE"),"hpux",4) ==0 ) {
      gSystem->Load("SFSFSpring.sl");
    }
    else {
      gSystem->Load("SFSFSpring.so");
   }
  }

  jsf = new JSFSteer();    // required to read parameter from jsf.conf
  bases = new SFSFBases();

//  bases->SetNoOfSample(10000);
//  bases->SetIteration1( 0.2, 10);
//  bases->SetIteration2( 0.1, 10);
  bases->Bases();
  bases->Bh_plot();
  bases->Userout();

  bases->Write();
  file->Write();

}

