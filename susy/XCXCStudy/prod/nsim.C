TFile *file;

void nsim()
{
  gROOT->Reset();
  file = new TFile("xcxcsim.run2.root","RECREATE");  // Output file

  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
  spring = new XCXCSpring();
  hdr    = new JSFHadronizer();
  sim    = new JSFQuickSim();

  Int_t maxevt = 200;      // Number of events.
  jsf->Initialize();

  //  Get seeds of prevous run from a file, jsf.root.
  //  This part must be executed prior to the begin run.
  TFile *flast = new TFile("xcxcsim.root","READ");
  jsf->GetLastRunInfo(flast);       // Get seeds of last run. 
  flast->Close(); 

  spring->ReadBases("bases.root");  // Bases must be initialized after 
                                    // GetLastRunInfo  
  printf(" Roots is %g\n",((XCXCBases*)spring->GetBases())->GetRoots());

  // Begin run
  jsf->BeginRun(31);      // Set run number to 31.
  for(Int_t ev=1;ev<=maxevt;ev++){
    printf(" start event %d\n",ev);
    if( !jsf->Process(ev) ) break;
    jsf->FillTree();
    jsf->Clear();
  }
  jsf->Terminate();
  file->Write();
}

