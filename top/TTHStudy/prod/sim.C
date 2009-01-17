{
  gROOT->Reset();
  TFile file("tthsim.root","RECREATE");  // Output file
 
  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
 
  spring = new TTHSpring();
  char *bsfile = jsf->Env()->GetValue("JSFGUI.Spring.BasesFile","bases.root");
  spring->ReadBases(bsfile);
 
  printf(" Roots is %g\n",((TTHBases*)spring->GetBases())->GetRoots());

  hdr=new JSFHadronizer();
  sim=new JSFQuickSim();
  
//  full->SetMakeBranch(kFALSE);   // suppress output of EventBuf 
//  hdr->SetMakeBranch(kFALSE);    // suppress output of EventBuf 
//  sim->SetMakeBranch(kFALSE);    // suppress output of EventBuf

  Int_t maxevt=50000;      // Number of events.
//  Int_t maxevt=2000;      // Number of events.
//  Int_t maxevt=200;      // Number of events.
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
