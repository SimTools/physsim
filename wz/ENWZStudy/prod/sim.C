{
  gROOT->Reset();
  TFile file("enwzsim.root","RECREATE");  // Output file
 
  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
 
  spring = new ENWZSpring();
  spring->ReadBases("bases.root");
 
  printf(" Roots is %g\n",((ENWZBases*)spring->GetBases())->GetRoots());

  hdr=new JSFHadronizer();
  sim=new JSFQuickSim();
  
//  full->SetMakeBranch(kFALSE);   // suppress output of EventBuf 
//  hdr->SetMakeBranch(kFALSE);    // suppress output of EventBuf 
//  sim->SetMakeBranch(kFALSE);    // suppress output of EventBuf

  Int_t maxevt=200;      // Number of events.
  jsf->Initialize();

  jsf->BeginRun(30);      // Set run number to 30.
  ev = 1;
  while (1) {
    printf(" start event %d\n",ev);
    if (jsf->Process(ev)) {
      if (ev >= maxevt) break;
      ev++;
    } else continue;
    jsf->FillTree();
    jsf->Clear();
  }
  jsf->Terminate();
  file->Write();
}
