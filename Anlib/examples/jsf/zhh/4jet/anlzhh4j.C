//***************************************************************************
// anlzhh4j.C
//
// JSF macro for analyzing e+e- -> ZHH -> 6 jets process at JLC.
// This program uses the library physsim-99a-1 of K.Fujii.
//
// (Update Record)
//	19 Nov 1999	A.L.C. Sanchez	*Based on the examples provided
//					with the physsim-99a-1 library.
//					Still to be corrected to suit the
//					desired JLC process.
//	22 Nov 1999	A.L.C. Sanchez	*Modified for testing ZHH4JAnalysis
//***************************************************************************

Int_t maxnevt = 200000;

Char_t *utils_so = "libS4Utils.so";
Char_t *anlib_so = "libAnlib.so";
Char_t *jsfanlib_so = "libJSFAnlib.so";
Char_t *zhhspr_so = "../../../../../higgs/ZHHStudy++/prod/ZHHSpring.so";
Char_t *procanl_so = "libZHH4JAnalysis.so";

Char_t *outputfile = "jsf.root";	// A file to output histograms
Char_t *inputfile = "../../../../../higgs/ZHHStudy++/prod/zhhsim.root";	// Input simulator file


int anlzhh4j()
{
  jsf  = new JSFSteer();			// Create JSF object

  gSystem->Load(utils_so);
  gSystem->Load(anlib_so);
  gSystem->Load(jsfanlib_so);
  gSystem->Load(zhhspr_so);
  gSystem->Load(procanl_so);

  file = new TFile(outputfile,"RECREATE");	// Outputfile
  fin  = new TFile(inputfile);			// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use.

  JSFSIMDST     *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  ZHH4JAnalysis *myanl  = new ZHH4JAnalysis("ZHH4JAnalysis","My Analysis");
  
  jsf->Initialize();

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis.

  // Adjust Cut.

  myanl->SetEvisLoCut(0.);	// minimum Evis
  myanl->SetEvisHiCut(420.);	// minimum Evis
  myanl->SetPtCut(999.);	// minimum Pt
  myanl->SetPlCut(999.);	// maximum Pl
  myanl->SetCosjetCut(0.99);	
  myanl->SetMinYcut(0.004);
  myanl->SetM2jCut( 30.);
  myanl->SetMM1Cut( 50.);
  myanl->SetMM2Cut(250.);

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok =0;
  for (Int_t evt=1; evt<=maxnevt; evt++) {
    if (!(jsf->GetEvent(evt))) break;		// Read in an event.
    if (!(jsf->Process(evt))) continue;		// Do SIMDST and ZHH4JAnalysis.
  }
  
  jsf->Terminate();
  return 0;
}
