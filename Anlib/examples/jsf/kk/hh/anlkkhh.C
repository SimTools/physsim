//***************************************************************************
// anlkkhh.C
//
// JSF macro for analyzing e+e- -> hh -> 4 jets process at JLC.
//
//***************************************************************************

Int_t maxnevt = 2000000;
//Int_t maxnevt = 200;

Int_t freq = 1;

Char_t *utils_so    = "libS4Utils.so";
Char_t *anlib_so    = "libAnlib.so";
Char_t *jsfanlib_so = "libJSFAnlib.so";

Char_t *kkhhspr_so  = "../../../../../xd/KKHHStudy/prod/libKKhhSpring.so";
Char_t *procanl_so  = "libKKHHAnalysis.so";

Char_t *zzspr_so    = "../../../../../wz/ZZStudy/prod/ZZSpring.so";
Char_t *zhspr_so    = "../../../../../higgs/ZHStudy/prod/ZHSpring.so";

Char_t *wwspr_so    = "../../../../../wz/WWStudy/prod/WWSpring.so";



//mh120
Double_t crossSection=5.77169;
Char_t *inputfile = "signal.mh120.root";	// Input simulator file
Char_t *outputfile = "analysis-kkhh-signal-mh120.root";	// A file to output histograms

/*
//mh160
Char_t *outputfile = "analysis-kkhh-signal-mh160.root";	// A file to output histograms
Double_t crossSection=5.05525;
Char_t *inputfile = "signal.mh160.root";	// Input simulator file
*/
/*
//mh200
Char_t *outputfile = "analysis-kkhh-signal-mh200.root";	// A file to output histograms
Double_t crossSection=4.21926;
Char_t *inputfile = "signal.mh200.root";	// Input simulator file
*/

/*
//mh120 ms2500
Char_t *outputfile = "analysis-kkhh-signal-2500.root";	// A file to output histograms
Double_t crossSection=0.968328;
Char_t *inputfile = "/data4/jlc/nicolas/datasim/kkhh.jsf.sim.ms2500.root";	// Input simulator file
*/

/*
//mh120 ms3000
Char_t *outputfile = "analysis-kkhh-signal-3000.root";	// A file to output histograms
Double_t crossSection=0.225202;
Char_t *inputfile = "/data4/jlc/nicolas/datasim/kkhh.jsf.sim.ms3000.root";	// Input simulator file
*/


/*
Char_t *outputfile = "analysis-bkgd-zz.root";	// A file to output histograms
Char_t *inputfile = "zzsim.root";	// Input simulator file
Double_t crossSection=206.666;
*/
/*
Char_t *outputfile = "analysis-bkdg-zh-mh120.root";	// A file to output histograms
Char_t *inputfile = "zhsim.mh120.root";	// Input simulator file
Double_t crossSection=18.3948;
*/
/*
Char_t *outputfile = "analysis-bkdg-zh-mh160.root";	// A file to output histograms
Char_t *inputfile = "zhsim.mh160.root";	// Input simulator file
Double_t crossSection=16.2704;
*/
/*
Char_t *outputfile = "analysis-bkdg-zh-mh200.root";	// A file to output histograms
Char_t *inputfile = "zhsim.mh200.root";	// Input simulator file
Double_t crossSection=0.100826;
*/
/*
Char_t *outputfile = "analysis-bkdg-ww.root";	// A file to output histograms
Char_t *inputfile = "wwsim.root";	// Input simulator file
Double_t crossSection=3833.33;
*/

/*
Char_t *outputfile = "analysis-bkdg-4b.root";	// A file to output histograms
Char_t *inputfile = "jsf.4b.root";	// Input simulator file
Double_t crossSection= 0.3779416E+01;
*/



int anlkkhh()
{
  jsf  = new JSFSteer();			// Create JSF object

  gSystem->Load(utils_so);
  gSystem->Load(anlib_so);
  gSystem->Load(jsfanlib_so);

  gSystem->Load(kkhhspr_so);
#if 0
  gSystem->Load(zzspr_so);
  gSystem->Load(zhspr_so);
  gSystem->Load(wwspr_so);
#endif
  gSystem->Load(procanl_so);

  file = new TFile(outputfile,"RECREATE");	// Outputfile
  fin  = new TFile(inputfile);			// Input simulator data

  jsf->SetInput(*fin);
  jsf->SetOutput(*file);

  // Define modules to use.
  
  JSFSIMDST     *simdst = new JSFSIMDST();	// Necessary to create SIMDST
  simdst->SetFile(file);			// since we analyze SIMDST
  simdst->NoReadWrite();			// instead of QuickSim data.
  
  KKHHAnalysis *myanl  = new KKHHAnalysis("KKHHAnalysis","Analysis");
  
  jsf->Initialize();

  JSFQuickSim *sim = (JSFQuickSim*)jsf->FindModule("JSFQuickSim");
  simdst->SetQuickSimParam(sim->Param());

  // Start analysis.

  // Adjust Cut.

  myanl->SetEvisCut(0.);	// minimum Evis
  myanl->SetPtCut(9999.);		// minimum Pt
  myanl->SetPlCut(9999.);	// maximum Pl
  //  myanl->SetCosjetCut(0.9);	
  myanl->SetCosjetCut(1.05);	
  myanl->SetMinYcut(0.004);
  myanl->SetM2jCut(16.);
  //myanl->SetM2jCut(32.);
  myanl->SetAcopCut(180.);

  myanl->SetMassTotDist(10);

  myanl->SetCrossSection(crossSection);

  jsf->BeginRun(1);				// Set run number to 1.
  Int_t nok =0;
  for (Int_t evt=1; evt<=maxnevt; evt++) {
    //    cout << "New event... " << endl;
    if (!(jsf->GetEvent(evt))) break;		// Read in an event.
    if (!(jsf->Process(evt))) continue;		// Do SIMDST and ZH4JAnalysis.
    if (!(gROOT->IsBatch())) {
      if (nok++%freq == 0) {
	myanl->DrawHist();	// Draw hists, if interactive.
	cout << "Event " << nok << endl;
      }
    }
    jsf->Clear();
  }

  if (!(gROOT->IsBatch())) myanl->DrawHist();	// Draw hists, if interactive.
  
  jsf->Terminate();
  file->Write();
  return 0;
}
