{
// Macro example for bases calculation.

  gROOT->Reset();
  TFile file("bases.root","RECREATE");

  jsf = new JSFSteer();    // required to read parameter from jsf.conf
  
  if( strncmp(gSystem->HostName(),"ccjlc",5)  != 0 ) {
    if( strncmp(gSystem->Getenv("OSTYPE"),"hpux",4) ==0 ) {
      gSystem->Load("XCXCSpring.sl");
    }
    else {
      gSystem->Load("XCXCSpring.so");
   }
  }
  bases = new XCXCBases();

//  bases->SetNoOfSample(5000);
//  bases->SetIteration1( 0.2, 1);
//  bases->SetIteration2( 0.1, 1);
  bases->Bases();
  bases->Bh_plot();
  bases->Userout();

  bases->Write();
  file->Write();

}



