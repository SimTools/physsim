{ 
  gROOT->Reset();
  TFile file("tbwsim02.root","RECREATE");  // Output file
 
  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
  spring = new TBWSpring("TBWSpring");
#if 0
  hdr=new JSFHadronizer();
  sim=new JSFQuickSim();
 
  TFile flast("tbwsim.root","READ");
  jsf->GetLastRunInfo(&flast);
  flast.Close();
#else
  TFile flast("tbwsim.root","READ");
  jsf->GetLastRunInfo(&flast);
  flast.Close();

  hdr=new JSFHadronizer();
  sim=new JSFQuickSim();
#endif

  spring->ReadBases("bases.root");

  Int_t maxevt=10000;      // Number of events.
  jsf->Initialize();


  printf(" Roots is %g\n",((TBWBases*)spring->GetBases())->GetRoots());
  
//  TBWree::SetBranchStyle(0);
//  full->SetMakeBranch(kFALSE);   // suppress output of EventBuf 
//  hdr->SetMakeBranch(kFALSE);    // suppress output of EventBuf 
//  sim->SetMakeBranch(kFALSE);    // suppress output of EventBuf

  jsf->BeginRun(31);      // Set run number to 30.
  for(Int_t ev=1;ev<=maxevt;ev++){
    printf(" start event %d\n",ev);
    if( !jsf->Process(ev) ) break;
    jsf->FillTree();
    jsf->Clear();
  }
  jsf->Terminate();
  file->Write();
}
