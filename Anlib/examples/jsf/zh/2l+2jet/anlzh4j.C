//***************************************************************************
// anlzh4j.C
//
// JSF macro for analyzing e+e- -> ZH -> 4 jets process at JLC.
// This program uses the library physsim-99a-1 of K.Fujii.
//
// (Update Record)
//	19 Nov 1999	A.L.C. Sanchez	*Based on the examples provided
//					with the physsim-99a-1 library.
//					Still to be corrected to suit the
//					desired JLC process.
//	22 Nov 1999	A.L.C. Sanchez	*Modified for testing ZH2L2JAnalysis
//***************************************************************************

Int_t maxnevt = 60;
Int_t freq = 10;

Char_t *utils_so = "libS4Utils.so";
Char_t *anlib_so = "libAnlib.so";
Char_t *jsfanlib_so = "libJSFAnlib.so";
Char_t *zhspr_so = "../../../../../higgs/ZHStudy/prod/ZHSpring.so";
Char_t *procanl_so = "libZH2L2JAnalysis.so";

Char_t *outputfile = "jsf.root";	// A file to output histograms
Char_t *inputfile = "../../../../../higgs/ZHStudy/prod/zhsim.root";	// Input simulator file


int anlzh4j()
{
  jsf  = new JSFSteer();			// Create JSF object

  gSystem->Load(utils_so);
  gSystem->Load(anlib_so);
  gSystem->Load(jsfanlib_so);
  gSystem->Load(zhspr_so);
  gSystem->Load(procanl_so);

  file = new TFile(outputfile,"RECREATE");	// Outputfile
  fin  = new TFile(inputfile);			// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use.

  JSFSIMDST     *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  ZH2L2JAnalysis *myanl  = new ZH2L2JAnalysis("ZH2L2JAnalysis","My Analysis");
  
  jsf->Initialize();

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis.

  // Adjust Cut.

  myanl->SetNtrackCut(3);
  myanl->SetEvisCut(0.);	// minimum Evis
  myanl->SetPtCut(0.);		// minimum Pt
  myanl->SetPlCut(999.);	// maximum Pl
  myanl->SetCosjetCut(0.99);	
  myanl->SetMinYcut(0.004);
  myanl->SetM2jCut(15.);
  myanl->SetAcopCut(180.);

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok =0;
  for (Int_t evt=1; evt<=maxnevt; evt++) {
    if (!(jsf->GetEvent(evt))) break;		// Read in an event.
    if (!(jsf->Process(evt))) continue;		// Do SIMDST and ZH2L2JAnalysis.
    if (!(gROOT->IsBatch())) {
      if (nok++%freq == 0) myanl->DrawHist();	// Draw hists, if interactive.
    }
    jsf->Clear();
  }
  if (!(gROOT->IsBatch())) myanl->DrawHist();	// Draw hists, if interactive.
  
  jsf->Terminate();
  return 0;
}
