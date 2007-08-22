Int_t maxevt = 10000;

int anlzx()
{  
  TFile *file;
  TFile *fin;

  Char_t *ofilen = "jsf.root";  // A file to output histograms
  Char_t *ifilen = "../../../../../../xd/RSZXStudy/prod/zxsim.root"; 
  //Char_t *ifilen = "../../../../../../xd/RSZXStudy/prod/zhsim.root"; 
  //Char_t *ifilen = "../../../../../../higgs/ZHStudy++/prod/zhsim.root"; 

  gSystem->Load("libJSFGenerator.so");
  gSystem->Load("libBasesSpring.so");
  gSystem->Load("libJSFTools.so");
  gSystem->Load("libJSFQuickSim.so");
  gSystem->Load("libS4Utils.so");
  gSystem->Load("libAnlib.so");
  gSystem->Load("libJSFAnlib.so");
  gSystem->Load("../../../../../../xd/RSZXStudy/prod/RSZXSpring.so");
  //gSystem->Load("../../../../../../higgs/ZHStudy++/prod/ZHSpring.so");
  gSystem->Load("libZX2A2JAnalysis.so");

  JSFSteer *jsf  = new JSFSteer();			// Create JSF object

  file = new TFile(ofilen,"RECREATE"); 	// Output file
  fin  = new TFile(ifilen);            	// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use. //

  JSFSIMDST    *simdst  = new JSFSIMDST();	// Necessary to create SIMDST 
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.
  
  ZX2A2JAnalysis *myanl = new ZX2A2JAnalysis("ZX2A2JAnalysis","My Analysis");

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

  myanl->SetEtrackCut( 0.100);
  myanl->SetEgammaCut(20.000);
  myanl->SetCosCone  ( 0.940);
  myanl->SetEconeCut ( 0.500);
  myanl->SetYcutCut  ( 0.004);
  cerr << " ------<Selection Cuts>--------------" << endl;
  cerr << " xEtrack  = " << myanl->GetEtrackCut() << endl
       << " xEgamma  = " << myanl->GetEgammaCut() << endl
       << " xCosCone = " << myanl->GetCosCone  () << endl
       << " xEcone   = " << myanl->GetEconeCut () << endl
       << " xYcut    = " << myanl->GetYcutCut  () << endl;
  cerr << " ------------------------------------" << endl;

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok = 0;
  for (Int_t ev=1; ev <= maxevt; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and ZX2A2JAnalysis.
     if (!(gROOT->IsBatch())) {
     }
     jsf->Clear();
  }
  jsf->Terminate();				// Terminate analysis.

  return 0;
}
