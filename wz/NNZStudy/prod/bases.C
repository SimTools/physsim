{
// Macro example for bases calculation.

  gROOT->Reset();
  TFile file("bases.root","RECREATE");

  jsf = new JSFSteer();    // required to read parameter from jsf.conf
  
  if( strncmp(gSystem->HostName(),"ccjlc",5)  != 0 ) {
    if( strncmp(gSystem->Getenv("OSTYPE"),"hpux",4) ==0 ) {
      gSystem->Load("NNZSpring.sl");
    }
    else {
      gSystem->Load("NNZSpring.so");
   }
  }
  bases = new NNZBases();

//  bases->SetNCALL(5000);
  bases->fPrintInfo=kTRUE;
  bases->fPrintHist=kTRUE;
//  bases->SetITMX1(1);
//  bases->SetITMX2(1);

  bases->DoBases();
  bases->Write();

  file->Write();

}



