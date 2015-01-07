{
  gROOT->Reset();
  TFile file("nnzhsim.root","RECREATE");  // Output file
 
  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
 
  spring = new NNZHSpring();
  spring->ReadBases("bases.root");
 
  printf(" Roots is %g\n",((NNZHBases*)spring->GetBases())->GetRoots());

  hdr=new JSFHadronizer();
  sim=new JSFQuickSim();
  //  stdhep = new JSFWriteStdHep();
  
//  full->SetMakeBranch(kFALSE);   // suppress output of EventBuf 
//  hdr->SetMakeBranch(kFALSE);    // suppress output of EventBuf 
//  sim->SetMakeBranch(kFALSE);    // suppress output of EventBuf

  Int_t maxevt=200;      // Number of events.
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
