//*************************************************************************
//* ===================
//*  anl6J.C JSF Macro
//* ===================
//*
//* (Description)
//*   JSF macro for analyze MC ttbar data to select 6-jet events.
//* (Update Recored)
//*   1999/08/03  K.Fujii       Original version
//*   2001/07/07  K.Ikematsu    Minor update
//*
//* $Id$
//*************************************************************************
//
Int_t maxevt = 50000;
Int_t freq   = 10;

int anl6J()
{
  TFile *file;
  TFile *fin;
  JSFSteer *jsf  = new JSFSteer();			// Create JSF object

  Char_t *outputfile="jsf.root";  // A file to output histograms
  Char_t *inputfile="../../../../../top/TTStudy/prod/ttsim.root";
  // Char_t *inputfile="ttsim.root";	// Input simulator file.

      gSystem->Load("libS4Utils.so");
      gSystem->Load("libAnlib.so");
      gSystem->Load("libJSFAnlib.so");
      gSystem->Load("../../../../../top/TTStudy/prod/TTSpring.so");
      gSystem->Load("libTT6JAnalysis.so");

  file = new TFile(outputfile,"RECREATE");  	// Output file
  fin  = new TFile(inputfile);            	// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use. //

  JSFSIMDST    *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  TT6JAnalysis *myanl  = new TT6JAnalysis("TT6JAnalysis","My Analysis");

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

  myanl->SetEvisCut(250.);
  myanl->SetPtCut(40.);
  //myanl->SetPlCut(20.);
  myanl->SetPlCut(40.);
  myanl->SetMinYcut(0.001);
  myanl->SetCosjetCut(1.);
#if 1
  myanl->SetM2jCut(24.);
  myanl->SetM3jCut(48.);
  //myanl->SetCosbwCut(-0.80);
  myanl->SetCosbwCut(1.);
#else
  myanl->SetM2jCut(12.);
  myanl->SetM3jCut(40.);
  myanl->SetCosbwCut(-0.80);
#endif

  jsf->BeginRun(1);      			// Set run number to 1.
  Int_t nok = 0;
  for (Int_t ev=1; ev <= maxevt; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and TT6JAnalysis.
     if (!(gROOT->IsBatch())) {
        if (nok++%freq == 0) myanl->DrawHist();	// Draw hists, if interactive.
     }
     jsf->Clear();
  }
  if (!(gROOT->IsBatch())) myanl->DrawHist();	// Draw hists, if interactive.

  jsf->Terminate();				// Terminate analysis.

  file->Write();
  return 0;
}
