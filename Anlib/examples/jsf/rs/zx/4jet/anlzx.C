int anlzx(Int_t maxevt = 10000)
{  
  TFile *file;
  TFile *fin;

  Char_t *ofilen = "jsf.root";  // A file to output histograms
  Char_t *ifilen = "../../../../../../xd/RSZXStudy/prod/zxsim.root"; 
  //Char_t *ifilen = "../../../../../../xd/RSZXStudy/prod/zhsim.root"; 

  gSystem->Load("libJSFGenerator.so");
  gSystem->Load("libBasesSpring.so");
  gSystem->Load("libJSFTools.so");
  gSystem->Load("libJSFQuickSim.so");
  gSystem->Load("libS4Utils.so");
  gSystem->Load("libAnlib.so");
  gSystem->Load("libJSFAnlib.so");
  gSystem->Load("../../../../../../xd/RSZXStudy/prod/RSZXSpring.so");
  gSystem->Load("libZX4JAnalysis.so");

  JSFSteer *jsf  = new JSFSteer();			// Create JSF object

  file = new TFile(ofilen,"RECREATE"); 	// Output file
  fin  = new TFile(ifilen);            	// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use. //

  JSFSIMDST    *simdst  = new JSFSIMDST();	// Necessary to create SIMDST 
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.
  
  ZX4JAnalysis *myanl = new ZX4JAnalysis("ZX4JAnalysis","My Analysis");

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

  myanl->SetEtrackCut( 0.100);
  myanl->SetYcutCut  ( 0.004);
  myanl->SetM2jCut   (12.000);
  cerr << " ------<Selection Cuts>--------------" << endl;
  cerr << " xEtrack  = " << myanl->GetEtrackCut() << endl
       << " xYcut    = " << myanl->GetYcutCut  () << endl 
       << " xM2j     = " << myanl->GetM2jCut   () << endl;
  cerr << " ------------------------------------" << endl;

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok = 0;
  for (Int_t ev=1; ev <= maxevt; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and ZX4JAnalysis.
     if (!(gROOT->IsBatch())) {
     }
     jsf->Clear();
  }
  jsf->Terminate();				// Terminate analysis.

  return 0;
}
