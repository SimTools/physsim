//#define BG_ANL "ENWZ"
//#define BG_ANL "NNWW"
//#define BG_ANL "EEWW"
//#define BG_ANL "WWZ"
//#define BG_ANL "WW"
Int_t freq   = 10;

int anl4J()
{  
  TFile *file;
  TFile *fin;
  JSFSteer *jsf  = new JSFSteer();			// Create JSF object

#ifndef BG_ANL
  Char_t *outputfile="jsf.whwh.root";  // A file to output histograms
  Char_t *inputfile="../../../../../lht/WHWHStudy/prod/whwhsim.root"; 
#else
#if BG_ANL == "ENWZ"
  Char_t *outputfile="jsf.enwz.root";  // A file to output histograms
  Char_t *inputfile="../../../../../wz/ENWZStudy/prod/enwzsim.root"; 
#endif
#if BG_ANL == "NNWW"
  Char_t *outputfile="jsf.nnww.root";  // A file to output histograms
  Char_t *inputfile="../../../../../wz/NNWWStudy/prod/nnwwsim.root"; 
#endif
#if BG_ANL == "EEWW"
  Char_t *outputfile="jsf.eeww.root";  // A file to output histograms
  Char_t *inputfile="../../../../../wz/EEWWStudy/prod/eewwsim.root"; 
#endif
#if BG_ANL == "WWZ"
  Char_t *outputfile="jsf.wwz.root";  // A file to output histograms
  Char_t *inputfile="../../../../../wz/WWZStudy/prod/wwzsim.root"; 
#endif
#if BG_ANL == "WW"
  Char_t *outputfile="jsf.ww.root";  // A file to output histograms
  Char_t *inputfile="../../../../../wz/WWStudy/prod/wwsim.root"; 
#endif
#endif
      gSystem->Load("libS4Utils.so");
      gSystem->Load("libAnlib.so");
      gSystem->Load("libJSFAnlib.so");
#ifndef BG_ANL
      gSystem->Load("../../../../../lht/WHWHStudy/prod/WHWHSpring.so");
#else
#if BG_ANL == "ENWZ"
      gSystem->Load("../../../../../wz/ENWZStudy/prod/ENWZSpring.so");
#endif
#if BG_ANL == "NNWW"
      gSystem->Load("../../../../../wz/NNWWStudy/prod/NNWWSpring.so");
#endif
#if BG_ANL == "EEWW"
      gSystem->Load("../../../../../wz/EEWWStudy/prod/EEWWSpring.so");
#endif
#if BG_ANL == "WWZ"
      gSystem->Load("../../../../../wz/WWZStudy/prod/WWZSpring.so");
#endif
#if BG_ANL == "WW"
      gSystem->Load("../../../../../wz/WWStudy/prod/WWSpring.so");
#endif
#endif
      gSystem->Load("libWHWH4JAnalysis.so");

  file = new TFile(outputfile,"RECREATE");  	// Output file
  fin  = new TFile(inputfile);            	// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  Int_t nevent=jsf->Env()->GetValue("JSFSteer.Nevent",1000000);  
  Int_t minevt=1;
  Int_t maxevt=minevt+nevent;

  // Define modules to use. //

  JSFSIMDST    *simdst = new JSFSIMDST();	// Necessary to create SIMDST 
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.
  
  WHWH4JAnalysis *myanl = new WHWH4JAnalysis("WHWH4JAnalysis","My Analysis");

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

  myanl->SetEcm(1000.); // Ecm set by hand
  myanl->SetEvisLoCut(20.);
  myanl->SetEvisHiCut(900.);
  myanl->SetPtCut(0.);
  myanl->SetPlCut(9999.);
  myanl->SetElCut(25.);
  myanl->SetCosjetCut(0.95);
  myanl->SetCoswCut(0.95);
  myanl->SetMinYcut(0.01);
  myanl->SetM2jCut(30.);
  myanl->SetMM1Cut(91.);
  myanl->SetMM2Cut(91.);
  myanl->SetAcopCut(30);

  jsf->BeginRun(1);      			// Set run number to 1.  
  Int_t nok = 0;
  for (Int_t ev=minevt; ev <= maxevt; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and WHWH4JAnalysis.
     jsf->Clear();
  }
  jsf->Terminate();				// Terminate analysis.

  //file->Write();
  return 0;
}
