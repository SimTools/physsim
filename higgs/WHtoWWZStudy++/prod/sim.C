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

  file=new TFile("whsim.root","RECREATE");  // Output file

  // Define modules.  JSFLCFULL must be declared to use JSFHadronizer and 
  // JSFQuickSim, which uses LCLIB libraries.  Each modules are executed
  // in the order according to the order of definition.
  jsf    = new JSFSteer();
  full   = new JSFLCFULL();
  spring = new WHSpring(); 
  pythia = new JSFHadronizer();
#if 1
  sim    = new JSFQuickSim();
#else
  stdhep=new JSFWriteStdHep();
#endif

  //Int_t maxevt=10000;      // Number of event 
  Int_t maxevt=100;     // Number of event is 10.
  jsf->Initialize();

  //spring->GetBases()->Bases();
  spring->ReadBases("bases.root");  // ReadBases must be called before BeginRun
  //spring->Bases()->SetSeed(12345);  // Set seed for Spring
  spring->SetPrintInfo(kTRUE);
  spring->SetPrintHist(kTRUE);

  jsf->BeginRun(30);      // Set run number to 30.

  //  ------------------------------------------------------------
  //  Event loop.
  //  ------------------------------------------------------------

  for(Int_t ev=1;ev<=maxevt;ev++){
    printf(" start event %d ",ev);

    if( !jsf->Process(ev)) break;  

    printf("Processed event %d ",ev);

    jsf->FillTree();
    jsf->Clear();

    printf(" End event %d \n",ev);
  }
  
  //  ------------------------------------------------------------

  jsf->Terminate();
  file->Write();
}
