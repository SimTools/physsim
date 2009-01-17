#define BG_ANL "NNH"
//#define BG_ANL "ZH"
//#define BG_ANL "ZZ"
//#define BG_ANL "NNZ"
//#define BG_ANL "TT"
Int_t freq   = 10;

int anl2J()
{  
  TFile *file;
  TFile *fin;
  JSFSteer *jsf  = new JSFSteer();			// Create JSF object

#ifndef BG_ANL
  Char_t *inputfile="../../../../../lht/ZHAHStudy/prod/zhahsim.root"; 
  Char_t *outputfile="jsf.zhah.root";  // A file to output histograms
#else
#if BG_ANL == "NNH"
  Char_t *inputfile="../../../../../higgs/NNHStudy++/prod/nnhsim.root"; 
  Char_t *outputfile="jsf.nnh.root";  // A file to output histograms
#endif
#if BG_ANL == "ZH"
  Char_t *inputfile="../../../../../higgs/ZHStudy++/prod/zhsim.root"; 
  Char_t *outputfile="jsf.zh.root";  // A file to output histograms
#endif
#if BG_ANL == "ZZ"
  Char_t *inputfile="../../../../../wz/ZZStudy/prod/zzsim.root"; 
  Char_t *outputfile="jsf.zz.root";  // A file to output histograms
#endif
#if BG_ANL == "NNZ"
  Char_t *inputfile="../../../../../wz/NNZStudy/prod/nnzsim.root"; 
  Char_t *outputfile="jsf.nnz.root";  // A file to output histograms
#endif
#if BG_ANL == "TT"
  Char_t *inputfile="../../../../../top/TTStudy/prod/ttsim.root"; 
  Char_t *outputfile="jsf.tt.root";  // A file to output histograms
#endif
#endif

      gSystem->Load("libS4Utils.so");
      gSystem->Load("libAnlib.so");
      gSystem->Load("libJSFAnlib.so");
#ifndef BG_ANL
      gSystem->Load("../../../../../lht/ZHAHStudy/prod/ZHAHSpring.so");
#else
#if BG_ANL == "NNH"
      gSystem->Load("../../../../../higgs/NNHStudy++/prod/NNHSpring.so");
#endif
#if BG_ANL == "ZH"
      gSystem->Load("../../../../../higgs/ZHStudy++/prod/ZHSpring.so");
#endif
#if BG_ANL == "ZZ"
      gSystem->Load("../../../../../wz/ZZStudy/prod/ZZSpring.so");
#endif
#if BG_ANL == "NNZ"
      gSystem->Load("../../../../../wz/NNZStudy/prod/NNZSpring.so");
#endif
#if BG_ANL == "TT"
      gSystem->Load("../../../../../top/TTStudy/prod/TTSpring.so");
#endif
#endif
      gSystem->Load("libZHAH2JAnalysis.so");

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
  
  ZHAH2JAnalysis *myanl = new ZHAH2JAnalysis("ZHAH2JAnalysis","My Analysis");

  Double_t ecm = 500.;
  myanl->SetEcm(ecm);

  jsf->Initialize();             		// JSF Module initialization.

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis. //

  // Adjust Cut //

#if 0
  myanl->SetNtrackCut(25.);
  myanl->SetEvisLoCut(20.);
  myanl->SetEvisHiCut(400.);
  myanl->SetPtCut(80.);
  myanl->SetPlCut(9999.);
  myanl->SetElCut(20.); 
  myanl->SetCosjetCut(0.999);
  myanl->SetCoshCut(0.8);
  myanl->SetMinYcut(0.01);
#else
  myanl->SetNtrackCut(5.);
  myanl->SetEvisLoCut(0.);
  myanl->SetEvisHiCut(500.);
  myanl->SetPtCut(0.);
  myanl->SetPlCut(99999.);
  myanl->SetElCut(99999.); 
  myanl->SetCosjetCut(1.0);
  myanl->SetCoshCut(1.0);
  myanl->SetMinYcut(0.01);
#endif
  myanl->SetM2jLoCut(100.);
  myanl->SetM2jHiCut(140.);
#if 0
  myanl->SetMM1Cut(70.);
  myanl->SetMM2Cut(110.);
#else
  myanl->SetMM1Cut(91.);
  myanl->SetMM2Cut(91.);
#endif
  myanl->SetBtagNsig(2.);
  myanl->SetBtagNoffv(1);

  jsf->BeginRun(1);      			// Set run number to 1.  

  Int_t nok = 0;
  for (Int_t ev=minevt; ev <= maxevt; ev++) {
     if (!(jsf->GetEvent(ev))) break;		// Read in an event.
     if (!(jsf->Process(ev))) continue;		// Do SIMDST and ZHAH2JAnalysis.
     jsf->Clear();
  }
  jsf->Terminate();				// Terminate analysis.

  file->Write();
  return 0;
}
