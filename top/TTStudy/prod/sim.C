{ gROOT->Reset();
  TFile file("ttsim.root","RECREATE");  // Output file
 
  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
 
  if( strncmp(gSystem->HostName(),"ccjlc",5)  != 0 ) {
    if( strncmp(gSystem->Getenv("OSTYPE"),"hpux",4) ==0 ) {
      gSystem->Load("TTSpring.sl");
    }
    else {
      gSystem->Load("TTSpring.so");
   }
  }
 
  spring = new TTSpring();
  spring->ReadBases("bases.root");
 
  printf(" Roots is %g\n",((TTBases*)spring->GetBases())->GetRoots());

  hdr=new JSFHadronizer();
  sim=new JSFQuickSim();
  
//  full->SetMakeBranch(kFALSE);   // suppress output of EventBuf 
//  hdr->SetMakeBranch(kFALSE);    // suppress output of EventBuf 
//  sim->SetMakeBranch(kFALSE);    // suppress output of EventBuf

  Int_t maxevt=2000;      // Number of events.
  jsf->Initialize();

  jsf->BeginRun(30);      // Set run number to 30.
  for(Int_t ev=1;ev<=maxevt;ev++){
    printf(" start event %d\n",ev);
    if( !jsf->Process(ev) ) break;
    jsf->FillTree();
    jsf->Clear();
  }

  jsf->Terminate();
  file->Write();

}


