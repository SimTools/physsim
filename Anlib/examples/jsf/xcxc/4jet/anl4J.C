Int_t freq   = 10;

int anl4J()
{  
  TFile *file;
  TFile *fin;
  JSFSteer *jsf  = new JSFSteer();			// Create JSF object

  Char_t *outputfile="jsf.root";  // A file to output histograms
  Char_t *inputfile="../../../../../susy/XCXCStudy/prod/xcxcsim.root"; 
  // Char_t *inputfile="xcxcsim.root";	// Input simulator file.

      gSystem->Load("libS4Utils.so");
      gSystem->Load("libAnlib.so");
      gSystem->Load("libJSFAnlib.so");
      gSystem->Load("../../../../../susy/XCXCStudy/prod/XCXCSpring.so");
      gSystem->Load("libXCXC4JAnalysis.so");

  file = new TFile(outputfile,"RECREATE");  	// Output file
  fin  = new TFile(inputfile);            	// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  Int_t nevent=jsf->Env()->GetValue("JSFSteer.Nevent",1000000);  
  Int_t minevt=1;
  Int_t maxevt=minevt+nevent;

  // Define modules to use. //

  JSFSIMDST    *simdst = new JSFSIMDST();	// Necessary to create SIMDST 
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.
  
  XCXC4JAnalysis *myanl = new XCXC4JAnalysis("XCXC4JAnalysis","My Analysis");

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

  myanl->SetEvisLoCut(20.);
  myanl->SetEvisHiCut(400.);
  myanl->SetPtCut(0.);
  myanl->SetPlCut(9999.);
  myanl->SetElCut(25.);
  myanl->SetCosjetCut(0.80,0.95);
  myanl->SetCoswCut(0.90);
  myanl->SetMinYcut(0.01);
#if 0
  myanl->SetM2jLoCut(10.);
  myanl->SetM2jHiCut(10.);
#else
  myanl->SetM2jLoCut(30.);
  myanl->SetM2jHiCut(30.);
#endif
  myanl->SetMM1Cut(70.);
  myanl->SetMM2Cut(120.);
  myanl->SetAcopCut(30.);

  jsf->BeginRun(1);      			// Set run number to 1.  
  Int_t nok = 0;
  for (Int_t ev=minevt; ev <= maxevt; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and XCXC4JAnalysis.
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
