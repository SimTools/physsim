//* $Id$
//*************************************************************************
//* =====================
//*  anl2L.C JSF Macro
//* =====================
//*
//* (Description)
//*   JSF macro to analyze sfsf data to select acoplanar 2-lepton events.
//* (Update Recored)
//*   2002/11/02  K.Fujii	Original version.
//*
//*************************************************************************

Int_t freq   = 10;

int anl2Lbg()
{
  TFile *file;
  TFile *fin;
  JSFSteer *jsf  = new JSFSteer();			// Create JSF object

  Char_t *outputfile="jsfwwbg.root";  // A file to output histograms
  // Char_t *inputfile="../../../../../wz/WWStudy/prod/wwsim.p+0.90.mutauonly.root";
  Char_t *inputfile="../../../../../wz/WWStudy/prod/wwsim.p+0.00.mutauonly.root";

  gSystem->Load("libS4Utils.so");
  gSystem->Load("libAnlib.so");
  gSystem->Load("libJSFAnlib.so");
  gSystem->Load("../../../../../wz/WWStudy/prod/WWSpring.so");
  gSystem->Load("libSFSF2LAnalysis.so");

  file = new TFile(outputfile,"RECREATE");  	// Output file
  fin  = new TFile(inputfile);            	// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  Int_t nevent=jsf->Env()->GetValue("JSFSteer.Nevent",100000);  
  Int_t minevt=1;
  Int_t maxevt=minevt+nevent;

  // Define modules to use. //

  JSFSIMDST    *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  SFSF2LAnalysis *myanl  = new SFSF2LAnalysis("SFSF2LAnalysis","My Analysis");

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

  myanl->SetAcopCut(30.);
  myanl->SetEvisLoCut(20.);
  myanl->SetEvisHiCut(250.);
  myanl->SetElepLoCut(5.);
  myanl->SetElepHiCut(125.);

  myanl->SetSpringName("WWSpring");

  jsf->BeginRun(1);      			// Set run number to 1.
  Int_t nok = 0;
  for (Int_t ev=1; ev <= nevent; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and SFSF2LAnalysis.
     if (!(gROOT->IsBatch())) {
        if (nok++%freq == 0) myanl->DrawHist();	// Draw hists, if interactive.
     }
     jsf->Clear();
  }
  if (!(gROOT->IsBatch())) myanl->DrawHist();	// Draw hists, if interactive.

  jsf->Terminate();				// Terminate analysis.

  //file->Write();
  return 0;
}
