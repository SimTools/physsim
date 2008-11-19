//***************************************************************************
// anlzhh6j.C
//
// JSF macro for analyzing e+e- -> ZHH -> 6 jets process at JLC.
// This program uses the library physsim-99a-1 of K.Fujii.
//
// (Update Record)
//	19 Nov 1999	A.L.C. Sanchez	*Based on the examples provided
//					with the physsim-99a-1 library.
//					Still to be corrected to suit the
//					desired JLC process.
//	22 Nov 1999	A.L.C. Sanchez	*Modified for testing ZHH6JAnalysis
//***************************************************************************

Int_t maxnevt = 200000;

Char_t *utils_so = "libS4Utils.so";
Char_t *anlib_so = "libAnlib.so";
Char_t *jsfanlib_so = "libJSFAnlib.so";
Char_t *zhhspr_so = "../../../../../higgs/ZHHStudy++/prod/ZHHSpring.so";
Char_t *procanl_so = "libZHH6JAnalysis.so";

Char_t *outputfile = "jsf.root";	// A file to output histograms
Char_t *inputfile = "../../../../../higgs/ZHHStudy++/prod/zhhsim.root";	// Input simulator file


int anlzhh6j()
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

  ZHH6JAnalysis *myanl  = new ZHH6JAnalysis("ZHH6JAnalysis","My Analysis");
  
  jsf->Initialize();

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis.

  // Adjust Cut.

  myanl->SetEvisLoCut(300.);	// minimum Evis
  myanl->SetEvisHiCut(600.);	// maximum Evis
  myanl->SetPtCut(999.);	// maximum Pt
  myanl->SetPlCut(999.);	// maximum Pl
  myanl->SetCosjetCut(1.00);	
  myanl->SetMinYcut(0.004);
  myanl->SetM2jCut( 80.);
  myanl->SetMM1Cut(-10.);
  myanl->SetMM2Cut(250.);
#if 0
  myanl->SetBtagNsig  (1.0); // loose b-tag used to tag b's
  myanl->SetBtagNoffv (2);   //
  myanl->SetBTtagNsig (3.0); // tight b-tag used to veto any b in W
  myanl->SetBTtagNoffv(99);   // 
#else
  myanl->SetBtagNsig  (2.5); // loose b-tag used to tag b's
  myanl->SetBtagNoffv (0);   //
  myanl->SetBTtagNsig (90.0); // tight b-tag used to veto any b in W
  myanl->SetBTtagNoffv(99);   // 
#endif

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok =0;
  for (Int_t evt=1; evt<=maxnevt; evt++) {
    if (!(jsf->GetEvent(evt))) break;		// Read in an event.
    if (!(jsf->Process(evt))) continue;		// Do SIMDST and ZHH6JAnalysis.
  }
  
  jsf->Terminate();
  return 0;
}
