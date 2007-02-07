//***************************************************************************
// anlax.C
//
// JSF macro for analyzing e+e- -> AX -> 4 jets process at JLC.
// This program uses the library physsim-99a-1 of K.Fujii.
//
// (Update Record)
//	19 Nov 1999	A.L.C. Sanchez	*Based on the examples provided
//					with the physsim-99a-1 library.
//					Still to be corrected to suit the
//					desired JLC process.
//	22 Nov 1999	A.L.C. Sanchez	*Modified for testing AXAnalysis
//***************************************************************************

Int_t maxnevt = 20000;

Char_t *outputfile = "jsf.root";	// A file to output histograms
Char_t *inputfile  = "../../../../../../xd/RSAXStudy/prod/axsim.root";	// Input simulator file

int anlax()
{
  jsf  = new JSFSteer();			// Create JSF object

  gSystem->Load("libS4Utils.so");
  gSystem->Load("libAnlib.so");
  gSystem->Load("libJSFAnlib.so");
  gSystem->Load("../../../../../../xd/RSAXStudy/prod/RSAXSpring.so");
  gSystem->Load("libAXAnalysis.so");

  file = new TFile(outputfile,"RECREATE");	// Outputfile
  fin  = new TFile(inputfile);			// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use.

  JSFSIMDST     *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  AXAnalysis *myanl  = new AXAnalysis("AXAnalysis","My Analysis");
  
  jsf->Initialize();

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis.

  // Adjust Cut.

  myanl->SetEvisCut(0.);	// minimum Evis
  myanl->SetPtCut(999.);		// maximum Pt
  myanl->SetPlCut(999.);	// maximum Pl
  myanl->SetCosjetCut(0.999);	
  myanl->SetMinYcut(0.004);
  myanl->SetM2jCut(25.);

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok =0;
  for (Int_t evt=1; evt<=maxnevt; evt++) {
    if (!(jsf->GetEvent(evt))) break;		// Read in an event.
    if (!(jsf->Process(evt))) continue;		// Do SIMDST and AXAnalysis.
    if (!(gROOT->IsBatch())) {
    }
    jsf->Clear();
  }
  
  jsf->Terminate();
  return 0;
}
