//***************************************************************
//*  sim.C
//*  An example script to run Spring, Hadronizer and QuickSimulator.
//*  Bases data, bases.root, is prepared by bases.C script in the
//*  same directory.
//*  
//*  24-July-1999 Akiya Miyamoto 
//*  31-July-1999 Akiya Miyamoto Use dynamic loading on system other than 
//*                              ccjlc.         
//*$Id$
//***************************************************************

TFile *file;  // file must be global to browse it at the end.

Int_t sim()
{
  gROOT->Reset();

  file=new TFile("zhhsim.root","RECREATE");  // Output file

  // Define modules.  JSFLCFULL must be declared to use JSFHadronizer and 
  // JSFQuickSim, which uses LCLIB libraries.  Each modules are executed
  // in the order according to the order of definition.
  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
  spring = new ZHHSpring(); 
  pythia = new JSFHadronizer();
  //  sim    = new JSFQuickSim();
  stdhep = new JSFWriteStdHep();

  Int_t maxevt=50000;      // Number of events
  //  Int_t maxevt=1000;        // Number of events
  jsf->Initialize();

  spring->ReadBases("bases.root");  // ReadBases must be called before BeginRun
  spring->SetPrintInfo(kTRUE);
  spring->SetPrintHist(kTRUE);

  //  jsf->BeginRun(30);      // Set run number to 30.
  jsf->BeginRun(20130513);      // Set run number to 30.

  //  ------------------------------------------------------------
  //  Event loop.
  //  ------------------------------------------------------------
  Int_t ev = 1;
  while (1) {
    printf(" start event %d ",ev);
    if (jsf->Process(ev)) {
      printf("Processed event %d ",ev);

      jsf->FillTree();
      jsf->Clear();

      printf(" End event %d \n",ev);
      if (ev >= maxevt) break;
      ev++;
    }
  }
  //  ------------------------------------------------------------
  jsf->Terminate();
  file->Write();
}
