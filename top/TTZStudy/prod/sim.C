{
  gROOT->Reset();
  TFile file("ttzsim01.root","RECREATE");  // Output file
 
  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
 
  spring = new TTZSpring();
  spring->ReadBases("bases.root");
 
  printf(" Roots is %g\n",((TTZBases*)spring->GetBases())->GetRoots());

  hdr=new JSFHadronizer();
  sim=new JSFQuickSim();
  
//  full->SetMakeBranch(kFALSE);   // suppress output of EventBuf 
//  hdr->SetMakeBranch(kFALSE);    // suppress output of EventBuf 
//  sim->SetMakeBranch(kFALSE);    // suppress output of EventBuf

//  Int_t maxevt=10000;      // Number of events.
  Int_t maxevt=700;      // Number of events.
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
