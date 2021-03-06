Int_t maxevt = 5000;
Int_t freq   = 10;

int anl4J()
{  
  TFile *file;
  TFile *fin;
  JSFSteer *jsf  = new JSFSteer();			// Create JSF object

  Char_t *outputfile="jsf.root";  // A file to output histograms
  Char_t *inputfile="../../../../../wz/WWStudy/prod/wwsim.root"; 
  // Char_t *inputfile="wwsim.root";	// Input simulator file.

      gSystem->Load("libS4Utils.so");
      gSystem->Load("libAnlib.so");
      gSystem->Load("libJSFAnlib.so");
      gSystem->Load("../../../../../wz/WWStudy/prod/WWSpring.so");
      gSystem->Load("libWW4JAnalysis.so");

  file = new TFile(outputfile,"RECREATE");  	// Output file
  fin  = new TFile(inputfile);            	// Input simulator data

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

  myanl->SetNtrackCut(25);
  myanl->SetEtrackCut(0.10);
  myanl->SetEvisCut(150.);
  myanl->SetPtCut(25.);
  myanl->SetPlCut(80.);
  myanl->SetMinYcut(0.004);
  myanl->SetNjetCut(4);
  myanl->SetEjetCut(5.00);
  myanl->SetCosjetCut(0.98);
  myanl->SetCoswCut(1.);
  myanl->SetM2jCut(12.);
  myanl->SetAcopCut(20.);

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
