{
// Macro example for bases calculation.

  gROOT->Reset();
  TFile file("bases.root","RECREATE");

  jsf = new JSFSteer();    // required to read parameter from jsf.conf
  
  if( strncmp(gSystem->HostName(),"ccjlc",5)  != 0 ) {
    if( strncmp(gSystem->Getenv("OSTYPE"),"hpux",4) ==0 ) {
      gSystem->Load("NNWWSpring.sl");
    }
    else {
      gSystem->Load("NNWWSpring.so");
   }
  }
  bases = new NNWWBases();


//  bases->SetNoOfSample(10000);
//  bases->SetIteration1( 0.2, 10);
//  bases->SetIteration2( 0.1, 10);
  bases->Bases();
  bases->Bh_plot();
  bases->Userout();
  bases->Write();

  file->Write();

}



