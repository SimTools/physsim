//*************************************************************************
//* ====================
//*  anl8J.C JSF Macro
//* ====================
//*
//* (Description)
//*   JSF macro for analyze MC ttbar data to select 8-jet events.
//* (Update Recored)
//*   1999/08/16  K.Ikematsu    Derived from anlL6J.C
//*   2001/07/07  K.Ikematsu    Minor update
//*
//* $Id$
//*************************************************************************

Int_t maxevt = 999999;
Int_t freq   = 10;

int anl8Jttz()
{
  TFile *file;
  TFile *fin;
  JSFSteer *jsf  = JSFSteer::Instance();	// Create JSF object

  Char_t *outputfile="jsf.ttz.root";  // A file to output histograms
  Char_t *inputfile="../../../../../top/TTZStudy/prod/ttzsim01.root";
  // Char_t *inputfile="tthsim.root";	// Input simulator file.

      gSystem->Load("libS4Utils.so");
      gSystem->Load("libAnlib.so");
      gSystem->Load("libJSFAnlib.so");
      gSystem->Load("../../../../../top/TTZStudy/prod/TTZSpring.so");
      gSystem->Load("libTTH8JAnalysis.so");
      TTH8JAnalysis::SetBasesName("TTZBases");

  file = new TFile(outputfile,"RECREATE");  	// Output file
  fin  = new TFile(inputfile);            	// Input simulator data
  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use. //

  JSFSIMDST    *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  TTH8JAnalysis *myanl  = new TTH8JAnalysis("TTH8JAnalysis","My Analysis");

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

  myanl->SetNtrackCut(25);
  myanl->SetEvisCut(300.);
  myanl->SetPtCut(100.);
  myanl->SetPlCut(9999.);
  myanl->SetNjetCut(8);
  myanl->SetMinYcut(0.001);
  myanl->SetEjetCut(5);
  myanl->SetCosjetCut(1.);
  myanl->SetCosbwCut(1.);
#if 0
  myanl->SetM2jCut(20.);
  myanl->SetM3jCut(30.);
#else
  myanl->SetM2jCut(30.);
  myanl->SetM3jCut(50.);
#endif
  myanl->SetThrustCut(0.8);

  myanl->SetBtagNsig  (1.0); // loose b-tag used to tag b's
  myanl->SetBtagNoffv (2);   //
  myanl->SetBTtagNsig (3.0); // tight b-tag used to veto any b in W
  myanl->SetBTtagNoffv(2);   // 

  jsf->BeginRun(1);      			// Set run number to 1.
  Int_t nok = 0;
  for (Int_t ev=1; ev <= maxevt; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and TTH8JAnalysis.
     if (!(gROOT->IsBatch())) {
     }
     jsf->Clear();
  }

  jsf->Terminate();				// Terminate analysis.

  file->Write();
  return 0;
}
