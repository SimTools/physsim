{ gROOT->Reset();
  TFile file("ttsim.root","RECREATE");  // Output file
 
  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
 
  spring = new TTSpring("TTSpring");
  spring->ReadBases("bases.root");
 
  printf(" Roots is %g\n",((TTBases*)spring->GetBases())->GetRoots());

  hdr=new JSFHadronizer();
  sim=new JSFQuickSim();
  
//  TTree::SetBranchStyle(0);
//  full->SetMakeBranch(kFALSE);   // suppress output of EventBuf 
//  hdr->SetMakeBranch(kFALSE);    // suppress output of EventBuf 
//  sim->SetMakeBranch(kFALSE);    // suppress output of EventBuf

  //Int_t maxevt=50000;      // Number of events.
  Int_t maxevt=200;      // Number of events.
  jsf->Initialize();

  jsf->BeginRun(30);      // Set run number to 30.

  Int_t ev = 1;
  while (1) {
    printf(" start event %d ",ev);

    if(jsf->Process(ev)) {
      printf("Processed event %d ",ev);

      jsf->FillTree();
      jsf->Clear();

      printf(" End event %d \n",ev);
      if (ev >= maxevt) break;
      ev++;
    }
  }
  jsf->Terminate();
  file->Write();
}
