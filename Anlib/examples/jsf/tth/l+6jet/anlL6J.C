//*************************************************************************
//* ====================
//*  anlL6J.C JSF Macro
//* ====================
//*
//* (Description)
//*   JSF macro for analyze MC ttbar data to select lepton + 6-jet events.
//* (Update Recored)
//*   1999/08/16  K.Ikematsu    Derived from anl6J.C
//*   2001/07/07  K.Ikematsu    Minor update
//*
//* $Id$
//*************************************************************************

Int_t maxevt = 99999999;
//Int_t maxevt = 1000;
Int_t freq   = 10;

int anlL6J()
{
  TFile *file;
  TFile *fin;
  JSFSteer *jsf  = new JSFSteer();			// Create JSF object

  Char_t *outputfile="jsf.root";  // A file to output histograms
  Char_t *inputfile="../../../../../top/TTHStudy/prod/tthsim.root";
  // Char_t *inputfile="tthsim.root";	// Input simulator file.

      gSystem->Load("libS4Utils.so");
      gSystem->Load("libAnlib.so");
      gSystem->Load("libJSFAnlib.so");
      gSystem->Load("../../../../../top/TTHStudy/prod/TTHSpring.so");
      gSystem->Load("libTTHL6JAnalysis.so");

  file = new TFile(outputfile,"RECREATE");  	// Output file
  fin  = new TFile(inputfile);            	// Input simulator data
  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use. //

  JSFSIMDST    *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  TTHL6JAnalysis *myanl  = new TTHL6JAnalysis("TTHL6JAnalysis","My Analysis");

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

  myanl->SetEvisCut(150.);
  myanl->SetPtCut(9999.);
  myanl->SetPlCut(9999.);
  myanl->SetNtrackCut(20);
  myanl->SetConeAngle(15.);
  myanl->SetElpCut(10.);
  myanl->SetEconeCut(5.);
  myanl->SetMinYcut(0.0001);
  myanl->SetNjetCut(6);
  myanl->SetEjetCut(5.);
  myanl->SetCosjetCut(1.);
#if 0
  myanl->SetM2jCut(40.);
  myanl->SetM3jCut(60.);
#else
  myanl->SetM2jCut(900.);
  myanl->SetM3jCut(900.);
#endif
  myanl->SetCosbwCut(1.);
  //myanl->SetThrustCut(0.8);
  myanl->SetThrustCut(1.);
#if 0
  myanl->SetBtagNsig  (2.); // loose b-tag used to tag b's
  myanl->SetBtagNoffv (0);   //
  myanl->SetBTtagNsig (10.0); // tight b-tag used to veto any b in W
  myanl->SetBTtagNoffv(3);   // 
#else
  myanl->SetBtagNsig  (2.0); // loose b-tag used to tag b's
  myanl->SetBtagNoffv (1);   //
  myanl->SetBTtagNsig (3.0); // tight b-tag used to veto any b in W
  myanl->SetBTtagNoffv(2);   // 
#endif

  jsf->BeginRun(1);      			// Set run number to 1.
  Int_t nok = 0;
  for (Int_t ev=1; ev <= maxevt; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and TTHL6JAnalysis.
     if (!(gROOT->IsBatch())) {
     }
     jsf->Clear();
  }

  jsf->Terminate();				// Terminate analysis.

  file->Write();
  return 0;
}
