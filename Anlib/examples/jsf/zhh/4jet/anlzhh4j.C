//#define BG_PROC "ZH"
//***************************************************************************
// anlzhh4j.C
//
// JSF macro for analyzing e+e- -> ZHH -> 6 jets process at JLC.
// This program uses the library physsim-99a-1 of K.Fujii.
//
// (Update Record)
//	19 Nov 1999	A.L.C. Sanchez	*Based on the examples provided
//					with the physsim-99a-1 library.
//					Still to be corrected to suit the
//					desired JLC process.
//	22 Nov 1999	A.L.C. Sanchez	*Modified for testing ZHH4JAnalysis
//***************************************************************************

Int_t maxnevt = 20000000;

Char_t *utils_so = "libS4Utils.so";
Char_t *anlib_so = "libAnlib.so";
Char_t *jsfanlib_so = "libJSFAnlib.so";
Char_t *procanl_so = "libZHH4JAnalysis.so";
#ifndef BG_PROC
Char_t *zhhspr_so = "../../../../../higgs/ZHHStudy++/prod/ZHHSpring.so";
Char_t *outputfile = "jsf.zhh.root";	// A file to output histograms
const Int_t kNfiles = 2;
Char_t *inputfile[] = {"../../../../../higgs/ZHHStudy++/prod/zhhsim.root", // Input simulator file
                       "../../../../../higgs/ZHHStudy++/prod/zhhsim001.root",
                       "../../../../../higgs/ZHHStudy++/prod/zhhsim002.root",
                       "../../../../../higgs/ZHHStudy++/prod/zhhsim003.root",
                       "../../../../../higgs/ZHHStudy++/prod/zhhsim004.root",
                       "../../../../../higgs/ZHHStudy++/prod/zhhsim005.root",
                       "../../../../../higgs/ZHHStudy++/prod/zhhsim006.root",
                       "../../../../../higgs/ZHHStudy++/prod/zhhsim007.root",
                       "../../../../../higgs/ZHHStudy++/prod/zhhsim008.root",
                       "../../../../../higgs/ZHHStudy++/prod/zhhsim009.root"};
#else
#if BG_PROC == "ZH"
Char_t *zhhspr_so = "../../../../../higgs/ZHStudy++/prod/ZHSpring.so";
Char_t *outputfile = "jsf.zh.root";	// A file to output histograms
const Int_t kNfiles = 1;
Char_t *inputfile[] = {"../../../../../higgs/ZHStudy++/prod/zhsim.root", // Input simulator file
                       "../../../../../higgs/ZHStudy++/prod/zhsim001.root",
                       "../../../../../higgs/ZHStudy++/prod/zhsim002.root",
                       "../../../../../higgs/ZHStudy++/prod/zhsim003.root",
                       "../../../../../higgs/ZHStudy++/prod/zhsim004.root",
                       "../../../../../higgs/ZHStudy++/prod/zhsim005.root",
                       "../../../../../higgs/ZHStudy++/prod/zhsim006.root",
                       "../../../../../higgs/ZHStudy++/prod/zhsim007.root",
                       "../../../../../higgs/ZHStudy++/prod/zhsim008.root",
                       "../../../../../higgs/ZHStudy++/prod/zhsim009.root"};
#elif BG_PROC == "TT"
Char_t *zhhspr_so = "../../../../../top/TTStudy/prod/TTSpring.so";
Char_t *outputfile = "jsf.tt.root";	// A file to output histograms
const Int_t kNfiles = 1;
Char_t *inputfile[] = {"../../../top/TTStudy/prod/data/ttsim.root",	// Input simulator file
                       "../../../top/TTStudy/prod/data/ttsim002.root",
                       "../../../top/TTStudy/prod/data/ttsim003.root",
                       "../../../top/TTStudy/prod/data/ttsim004.root",
                       "../../../top/TTStudy/prod/data/ttsim005.root",
                       "../../../top/TTStudy/prod/data/ttsim006.root",
                       "../../../top/TTStudy/prod/data/ttsim007.root",
                       "../../../top/TTStudy/prod/data/ttsim008.root",
                       "../../../top/TTStudy/prod/data/ttsim009.root",
                       "../../../top/TTStudy/prod/data/ttsim010.root"};
#endif
#endif

int anlzhh4j()
{
  jsf  = new JSFSteer();			// Create JSF object

  gSystem->Load(utils_so);
  gSystem->Load(anlib_so);
  gSystem->Load(jsfanlib_so);
  gSystem->Load(zhhspr_so);
  gSystem->Load(procanl_so);

  Int_t ifile = 0;
  file = new TFile(outputfile,"RECREATE");	// Outputfile
  fin  = new TFile(inputfile[ifile]);			// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use.

  JSFSIMDST     *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.

  ZHH4JAnalysis *myanl  = new ZHH4JAnalysis("ZHH4JAnalysis","My Analysis");
  Double_t ecm = 500.; // should be correctly set by hand
  myanl->SetEcm(ecm);
  
  jsf->Initialize();

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis.

  // Adjust Cut.

  myanl->SetEvisLoCut(0.);	// minimum Evis
  myanl->SetEvisHiCut(420.);	// minimum Evis
  myanl->SetPtCut(999.);	// minimum Pt
  myanl->SetPlCut(999.);	// maximum Pl
  myanl->SetCosjetCut(0.99);	
  myanl->SetMinYcut(0.01);
  myanl->SetM2jCut( 60.);
  myanl->SetMM1Cut( 50.);
  myanl->SetMM2Cut(250.);
  myanl->SetBtagNsig(2.0);
  myanl->SetBtagNoffv(1);

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok = 0;
  Int_t nrd = 0;
  while (nok < maxnevt) {
    if (!(jsf->GetEvent(++nrd))) {		// Read in an event.
      JSFSteer::EJSFReturnCode iret = jsf->GetReturnCode();
      if (iret & jsf->kJSFEOF) {
        cerr << " EOF of " << inputfile[ifile] << endl;
        cerr << "  reached after reading " << nrd-1 << " events" << endl;
	ifile++;
        if (ifile >= kNfiles) break;
        fin  = new TFile(inputfile[ifile]);			// Input simulator data
	jsf->SetInput(*fin);
	nrd = 1;
      } else {
        cerr << " Error reading event " << nrd << " of " << inputfile[ifile] << endl;
	cerr << " Terminate analysis" << endl;
	break;
      }
    } else {
        ++nok;
    }
    if (!(jsf->Process(nrd))) continue;	// Do SIMDST and ZHH4JAnalysis.
  }
  cerr << " Analyzed  " << nok << " events" << endl;
  
  jsf->Terminate();
  return 0;
}
