//***************************************************************************
// anlnnh2j.C
//
// JSF macro for analyzing e+e- -> NNH -> 2-lepton process at JLC.
//
// (Update Record)
//  2010/04/28 K.Fujii	Original version.
//***************************************************************************

Int_t maxnevt = 20000;
Int_t freq = 10;

Char_t *utils_so = "libS4Utils.so";
Char_t *anlib_so = "libAnlib.so";
Char_t *jsfanlib_so = "libJSFAnlib.so";
Char_t *procanl_so = "libNNH2JAnalysis.so";

Char_t *outputfile = "jsf.root";	// A file to output histograms
#if 1
Char_t *nnhspr_so = "../../../../../higgs/NNHStudy++/prod/NNHSpring.so";
Char_t *inputfile = "../../../../../higgs/NNHStudy++/prod/nnhsim.root";	// Input simulator file
#else
Char_t *nnhspr_so = "../../../../../wz/NNZStudyprod/NNZSpring.so";
Char_t *inputfile = "../../../../../wz/NNZStudyprod/nnzsim.root";	// Input simulator file
#endif


int anlnnh2j()
{
  jsf  = new JSFSteer();			// Create JSF object

  gSystem->Load(utils_so);
  gSystem->Load(anlib_so);
  gSystem->Load(jsfanlib_so);
  gSystem->Load(nnhspr_so);
  gSystem->Load(procanl_so);

  file = new TFile(outputfile,"RECREATE");	// Outputfile
  fin  = new TFile(inputfile);			// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use.

  JSFSIMDST     *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  NNH2JAnalysis *myanl  = new NNH2JAnalysis("NNH2JAnalysis","My Analysis");
  
  jsf->Initialize();

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis.

  // Adjust Cut.

  myanl->SetNtrackCut(5);
  myanl->SetEvisLoCut(0.);	// minimum Evis
  myanl->SetPtCut(999.);	// minimum Pt
  myanl->SetPlCut(999.);	// maximum Pl
  myanl->SetCosjetCut(0.99999);	// |cos(theta_lep)| maximum
  myanl->SetM2jCut(100.);	// mass window

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok =0;
  for (Int_t evt=1; evt<=maxnevt; evt++) {
    if (!(jsf->GetEvent(evt))) break;		// Read in an event.
    if (!(jsf->Process(evt))) continue;		// Do SIMDST and NNH2JAnalysis.
    jsf->Clear();
  }
  
  jsf->Terminate();
  return 0;
}
