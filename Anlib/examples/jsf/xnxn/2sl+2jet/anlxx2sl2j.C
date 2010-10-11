//***************************************************************************
// anlxx2sl2j.C
//
// JSF macro for analyzing e+e- -> XNXN -> (stau tau) (stau tau)
//
// (Update Record)
//  2010/10/10  K.Fujii		Original version.
//***************************************************************************

Int_t maxnevt = 20000;

Char_t *utils_so = "libS4Utils.so";
Char_t *anlib_so = "libAnlib.so";
Char_t *jsfanlib_so = "libJSFAnlib.so";
Char_t *xxspr_so = "../../../../../susy/XN1XN1Study/prod/XN1XN1Spring.so";
Char_t *procanl_so = "libXN1XN12SL2JAnalysis.so";

Char_t *outputfile = "jsf.root";
Char_t *inputfile = "../../../../../susy/XN1XN1Study/prod/xn1xn1sim.root";

int anlxx2sl2j()
{
  jsf  = new JSFSteer();			// Create JSF object

  gSystem->Load(utils_so);
  gSystem->Load(anlib_so);
  gSystem->Load(jsfanlib_so);
  gSystem->Load(xxspr_so);
  gSystem->Load(procanl_so);

  file = new TFile(outputfile,"RECREATE");	// Outputfile
  fin  = new TFile(inputfile);			// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use.

  JSFSIMDST     *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  XN1XN12SL2JAnalysis *myanl  = new XN1XN12SL2JAnalysis("XN1XN12SL2JAnalysis","My Analysis");
  Double_t Ecm = 500.;
  myanl->SetEcm(Ecm);
  
  jsf->Initialize();

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis.

  // Adjust Cut.

  myanl->SetNtrackCut(3);
  myanl->SetEvisLoCut(0.);	// minimum Evis
  myanl->SetPtCut(999.);	// minimum Pt
  myanl->SetPlCut(999.);	// maximum Pl
  myanl->SetCosjetCut(0.99);	
  myanl->SetMinYcut(0.004);
  myanl->SetM2jCut(50.);
  myanl->SetAcopCut(0.);

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok =0;
  for (Int_t evt=1; evt<=maxnevt; evt++) {
    if (!(jsf->GetEvent(evt))) break;		// Read in an event.
    if (!(jsf->Process(evt))) continue;		// Do SIMDST and XN1XN12SL2JAnalysis.
    if (!(gROOT->IsBatch())) {
    }
    jsf->Clear();
  }
  
  jsf->Terminate();
  return 0;
}
