Int_t maxevt = 5000;
Int_t freq   = 10;

int anl4J()
{  
  TFile *file;
  TFile *fin;
  JSFSteer *jsf;

  Char_t *outputfile="jsf.root";  // A file to output histograms
  Char_t *inputfile="../../../../wz/WWStudy/prod/wwsim.root"; 
  // Char_t *inputfile="wwsim.root";	// Input simulator file.

  if( strncmp(gSystem->HostName(),"ccjlc",5)  != 0 ) {
    if( strncmp(gSystem->Getenv("OSTYPE"),"hpux",4) ==0 ) {
      gSystem->Load("libAnlib.sl");
      gSystem->Load("../../../../wz/WWStudy/prod/WWSpring.sl");
      gSystem->Load("libWW4JAnalysis.sl");
    }
    else {
      gSystem->Load("libAnlib.so");
      gSystem->Load("../../../../wz/WWStudy/prod/WWSpring.so");
      gSystem->Load("libWW4JAnalysis.so");
   }
  }

  file = new TFile(outputfile,"RECREATE");  	// Output file
  fin  = new TFile(inputfile);            	// Input simulator data

  jsf  = new JSFSteer();			// Create JSF object
  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use. //

  JSFSIMDST    *simdst = new JSFSIMDST();	// Necessary to create SIMDST 
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.
  
  WW4JAnalysis *myanl  = new WW4JAnalysis("WW4JAnalysis","My Analysis");

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

  myanl->SetEvisCut(100.);
  myanl->SetPtCut(20.);
  myanl->SetPlCut(40.);
  myanl->SetCosjetCut(0.98);
  myanl->SetCoswCut(1.);
  myanl->SetMinYcut(0.004);
  myanl->SetM2jCut(12.);
  myanl->SetAcopCut(180.);

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok = 0;
  for (Int_t ev=1; ev <= maxevt; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and WW4JAnalysis.
     if (!(gROOT->IsBatch())) {
        if (nok++%freq == 0) myanl->DrawHist();	// Draw hists, if interactive.
     }
     jsf->Clear();
  }
  if (!(gROOT->IsBatch())) myanl->DrawHist();	// Draw hists, if interactive.

  jsf->Terminate();				// Terminate analysis.

  file->Write();
  return 0;
}
