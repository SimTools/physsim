//*************************************************************************
//* ================
//*  UserAnalysis.C
//* ================
//*
//* (Description)
//*    A very primitive sample script to study e+e- -> ttbar.
//* (Usage)
//*    	 $ jsf gui.C
//*    gui.C then invokes this script from within before event loop
//*    starts.
//* (Update Recored)
//*    1999/08/05  K.Fujii	Original version for 6-jet analysis.
//*
//*************************************************************************
//
#include <iostream.h>
//*------------------------*//
//* User Analysis          *//
//*------------------------*//

#if 1
Bool_t gDEBUG = kFALSE;
#else
Bool_t gDEBUG = kTRUE;
#endif

// Canvas //

TCanvas *cHist;
static const Int_t    kZoneX   = 4;	// No. X Zones in the Canvas
static const Int_t    kZoneY   = 4;	// No. Y Zones in the Canvas

// Hists //

TH1F *hStat;
TH1F *hNtracks;
TH1F *hEvis;
TH1F *hPt;
TH1F *hNjets;
TH1F *hEjet;
TH1F *hCosjet;
TH1F *hNsols;
TH1F *hChi2;
TH2F *hEw1Ew2;
TH2F *hCosbw1Cosbw2;
TH2F *hMw1Mw2;
TH2F *hMt1Mt2;
TH2F *hEvisPl;

// Event Information //
Double_t fEcm = 350.0;	// CM energy (GeV)
Int_t    fNtracks;	// track multiplicity
Double_t fEvis;		// visible energy
Double_t fPt;		// Pt
Double_t fPl;		// Pl
Int_t    fNjets;	// jet multiplicity
Double_t fYcut;		// y_cut to force the event to 4 jets

// Event selection condition //

static const Int_t maxcut = 50;
Char_t cutName[maxcut][100];	// Cut names

Int_t    xNtracks  =     25; 	// No. of tracks
Double_t xEtrack   =   0.10;	// track energy
Double_t xEvis     = 250.00;	// Minimum visible energy
Double_t xPt       =  40.00;	// Pt maximum
Double_t xPl       =  20.00;	// Pl maximum
Double_t xYcut     =  0.001;	// y_cut to force the event to 4 jets
Int_t    xNjets    =      6;	// No. of jets
Double_t xEjet     =   5.00;	// E_jet minimum
Double_t xCosjet   =   1.00;	// |cos(theta_j)| maximum
Double_t xCosbw    =  -0.80;	// cos(theta_bw) maximum
Double_t xM2j      =  12.00;	// |m_jj-m_W| maximum
Double_t xM3j      =  40.00;	// |m_3j-m_t| maximum

// Some constants //

static const Double_t kMassW   = 80.00; // W mass
static const Double_t kMassZ   = 91.19; // Z mass
static const Double_t kMasst   = 170.0; // top mass
static const Double_t kSigmaMw =   4.0; // W mass resolution
static const Double_t kSigmaMz =   4.0; // W mass resolution
static const Double_t kSigmaMt =  15.0; // top mass resolution

// Some global variables //

static Int_t Ngoods = 0;	// Number of good events

//_____________________________________________________________________
void UserInitialize()
{
  //  This function is called at the begining of the job or when
  //  "reset hist" action is selected in the gui menu.
  //  This is used to define/reset histograms.

  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat         = new TH1F("hStat","Cut Statistics",  20,   0.0,  20.0);
  hNtracks      = new TH1F("hNtracks","No. tracks" ,  50,   0.0, 200.0);
  hEvis         = new TH1F("hEvis","Visible energy",  80,   0.0, 400.0);
  hPt           = new TH1F("hPt","Missing Pt"      ,  60,   0.0, 120.0);
  hNjets        = new TH1F("hNjets","No. jets"     ,  20,   0.0,  20.0);
  hEjet         = new TH1F("hEjet","Jet energy"    ,  50,   0.0, 100.0);
  hCosjet       = new TH1F("hCosjet","cos(theta_j)",  50,  -1.0,  +1.0);
  hNsols        = new TH1F("hNsols","No. solutions",  20,   0.0,  20.0);
  hChi2         = new TH1F("hChi2","Chi2"          ,  50,   0.0,  50.0);
  hEw1Ew2       = new TH2F("hEw1Ew2","(E_w1,E_w2)" ,  
  				    50,  0.0, 200.0,  50,   0.0, 200.0);
  hCosbw1Cosbw2 = new TH2F("hCosbw1Cosbw2","(cos_bw1,cow_bw2)",
				    50, -1.0,  +1.0,  50,  -1.0,  +1.0);
  hMw1Mw2       = new TH2F("hMw1Mw2","(m_w1,m_w2)" , 
  				    60, 50.0, 110.0,  60,  50.0, 110.0);
  hEvisPl       = new TH2F("hEvisPl","(Evis,Pl)"   , 
  				    60,  0.0, 600.0,  50,-100.0,+100.0);
  hMt1Mt2       = new TH2F("hMt1Mt2","(m_t1,m_t2)" , 
  				    50,120.0, 220.0,  50, 120.0, 220.0);

  last->cd();
  return;
}
//_________________________________________________________
void DrawHist()
{
  //  This function is called to draw histograms during the interactive 
  //  session.  Thus you can see the accumulation of the histogram
  //  interactively.  

  TDirectory *last = gDirectory;
  if( !cHist ) {
    cHist = new TCanvas("cHist","Canvas 1",100, 100, kZoneX*200, kZoneY*200);
    cHist->Divide(kZoneX,kZoneY);
  } else {
    cHist->cd();
  }
  cHist->cd(1);		hStat->Draw();
  cHist->cd(2);		hNtracks->Draw();
  cHist->cd(3);		hEvis->Draw();
  cHist->cd(4);		hPt->Draw();
  cHist->cd(5);		hNjets->Draw();
  cHist->cd(6);		hEjet->Draw();
  cHist->cd(7);		hCosjet->Draw();
  cHist->cd(8);		hNsols->Draw();
  cHist->cd(9);		hChi2->Draw();
  cHist->cd(10);	hEw1Ew2->Draw();
  cHist->cd(11);	hCosbw1Cosbw2->Draw();
  cHist->cd(12);	hMw1Mw2->Draw();
  cHist->cd(13);	hEvisPl->Draw();
  cHist->cd(14);	hMt1Mt2->Draw();

  cHist->Update();
  last->cd();
}

//_________________________________________________________
void UserSetOptions()
{
  // This function is called only once, soon after jsf is started.
  // This function can be used to define parameters which is not 
  // defined in jsf.conf file.
}

//_________________________________________________________
void UserAnalysis()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  Char_t msg[60];

  // ---------------------
  // Analysis starts here.
  // ---------------------
  
  Float_t selid = -0.5;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) strcpy(&cutName[(Int_t)selid][0],"No cut");

  // Get event buffer and make combined tracks accessible.

  JSFSIMDST     *sds     = (JSFSIMDST*)jsf->FindModule("JSFSIMDST");
  JSFSIMDSTBuf  *evt     = (JSFSIMDSTBuf*)sds->EventBuf();
  Int_t          ntrks   = evt->GetNLTKCLTracks(); 	// No. of tracks 
  TClonesArray  *trks    = evt->GetLTKCLTracks(); 	// combined tracks
   
  // Select good tracks and store them in "TObjArray tracks".
  
  ANL4DVector qsum;
  TObjArray tracks(1000);
  fNtracks = 0;
  for ( Int_t i = 0; i < ntrks; i++ ) {
    JSFLTKCLTrack *t = (JSFLTKCLTrack*)trks->UncheckedAt(i);
    if ( t->GetE() > xEtrack ) {
      ANLTrack *qt = new ANLTrack(t);
      tracks.Add(qt); 		// track 4-momentum
      qsum += *qt;		// total 4-momentum
      fNtracks++;
    }				// *qt stays.
  }

  // Cut on No. of tracks.
  
  hNtracks->Fill(fNtracks);
  if ( fNtracks < xNtracks ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"N_tracks > %i",xNtracks);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  fEvis = qsum.E();		// E_vis
  fPt   = qsum.GetPt();		// P_t
  fPl   = qsum.Pz();		// P_l

  // Cut on Evis.

  hEvis->Fill(fEvis);
  if ( fEvis < xEvis ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"E_vis > %g",xEvis);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
 
  // Cut on Pt.

  hPt->Fill(fPt);
  if ( fPt > xPt ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Pt <= %g",xPt);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
 
  // Cut on Pl.

  if ( TMath::Abs(fPl) > xPl ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|Pl| <= %g",xPl);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Find jets.

  fYcut = xYcut;
  ANLJadeEJetFinder jclust(fYcut);
  jclust.Initialize(tracks);
  jclust.FindJets();
  fYcut  = jclust.GetYcut();
  fNjets = jclust.GetNjets();

  // Cut on No. of jets.
    
  hNjets->Fill(fNjets);
  if ( fNjets < xNjets ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Njets >= %i for Ycut = %g",xNjets,xYcut);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  // Now force the event to be xNjets.
  
  jclust.ForceNJets(xNjets);
  fNjets = jclust.GetNjets();
  fYcut  = jclust.GetYcut();

  // Make sure that No. of jets is xNjets.
    
  if ( fNjets != xNjets ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Njets = %i",xNjets);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  // Loop over jets and decide Ejet_min and |cos(theta_j)|_max.
  
  TObjArray &jets = jclust.GetJets();	// jets is an array of ANLJet's
  TIter nextjet(&jets);			// and nextjet is an iterator for it
  ANLJet *jetp;
  Double_t ejetmin = 999999.;
  Double_t cosjmax = 0.;
  Int_t nj = 0;
  while ((jetp = (ANLJet *)nextjet())) {
    ANLJet &jet = *jetp;
    if (gDEBUG) {
       cerr << "Jet " << ++nj << " : ";
       jet.DebugPrint();
    }
    Double_t ejet = jet().E();
    if (ejet < ejetmin) ejetmin = ejet;
    hEjet->Fill(ejet);			// Ejet
    Double_t cosj = jet.CosTheta();
    if (TMath::Abs(cosj) > TMath::Abs(cosjmax)) cosjmax = cosj;
    hCosjet->Fill(cosj);		// cos(theta_jet)
  }

  // Cut on Ejet_min.
  
  if ( ejetmin < xEjet ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Ejet > %g",xEjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on |cos(theta_j)|_max.
    
  if ( TMath::Abs(cosjmax) > xCosjet ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_j)| <= %g",xCosjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Find W and top candidates in given mass windows.
  // Avoid using indices since there might be empty slots.
    
  TObjArray solutions(10);
  ANLPairCombiner w1candidates(jets,jets);

  if (gDEBUG/* || kTRUE */) {
    cerr << "------------------------------------------" << endl;
    cerr << "- w1candidates:" << endl;
    w1candidates.DebugPrint();
  }
  Int_t nsol = 0;

  ANLPair *w1p, *w2p, *bbp;
  while ((w1p = (ANLPair *)w1candidates())) {
    ANLPair &w1 = *w1p;
    Double_t w1mass = w1().GetMass();
    if (TMath::Abs(w1mass - kMassW) > xM2j) continue;	// w1 candidate found
    w1.LockChildren();					// now lock w1 daughters
    ANLPairCombiner w2candidates(w1candidates);		// w2 after w1
    if (gDEBUG) {
      cerr << "-- w2candidates:" << endl;
      w2candidates.DebugPrint();
    }
    while ((w2p = (ANLPair *)w2candidates())) {
      ANLPair &w2 = *w2p;
      if (w2.IsLocked()) continue;			// skip if locked
      Double_t w2mass = w2().GetMass();
      if (TMath::Abs(w2mass - kMassW) > xM2j) continue;	// w2 candidate found
      w2.LockChildren();				// now lock w2 daughters
      ANLPairCombiner bbcandidates(w2candidates);	// bbar from the rest
      bbcandidates.Reset();				// bb can be before w1/2
      if (gDEBUG) {
        cerr << "---- bbcandidates:" << endl;
        bbcandidates.DebugPrint();
      }
      while ((bbp = (ANLPair *)bbcandidates())) {
        ANLPair &bb = *bbp;
        if (bb.IsLocked()) continue;			// skip if locked
        for (Int_t i = 0; i < 2; i++) {
          ANL4DVector *b1p = (ANL4DVector *)bb[i];	// b for w1
          ANL4DVector *b2p = (ANL4DVector *)bb[1-i];	// b for w2
          ANLPair *bw1p = new ANLPair(b1p,w1p);
          ANLPair *bw2p = new ANLPair(b2p,w2p);
          ANLPair &bw1  = *bw1p;
          ANLPair &bw2  = *bw2p;
          Double_t t1mass = bw1().GetMass();
          Double_t t2mass = bw2().GetMass();
          if (TMath::Abs(t1mass - kMasst) > xM3j ||
              TMath::Abs(t2mass - kMasst) > xM3j ) {
            delete bw1p;
            delete bw2p;
            continue;
          }						// tt candidate found
          Double_t chi2 = TMath::Power((w1mass - kMassW)/kSigmaMw,2.)
                        + TMath::Power((w2mass - kMassW)/kSigmaMw,2.)
                        + TMath::Power((t1mass - kMasst)/kSigmaMt,2.)
                        + TMath::Power((t2mass - kMasst)/kSigmaMt,2.);
          ANLPair *solp = new ANLPair(bw1p,bw2p,chi2);
          solutions.Add(solp);	// save this solution
          if (gDEBUG) { 
            Double_t cosbw1 = b1p->CosTheta(w1());
            Double_t cosbw2 = b2p->CosTheta(w2());
            cerr << " Solution ---" << ++nsol << " Chi2 = " << chi2 << endl;
            cerr << " (cosbw1,cosbw2) = (" << cosbw1 << ", " << cosbw2 << ")";
            cerr << " M_bb = " << bb().GetMass() << endl;
            cerr << " (M_w1,M_w2) = (" << w1mass << ", " << w2mass << ")";
            cerr << " (M_t1,M_t2) = (" << t1mass << ", " << t2mass << ")" << endl;
            cerr << " b1:"; b1p->DebugPrint();
            cerr << " b2:"; b2p->DebugPrint();
            cerr << " w1:"; w1().DebugPrint();
            cerr << " w2:"; w2().DebugPrint();
            // cerr << " b1p = " << ((void *)b1p) << " b2p = " << ((void *)b2p);
            // cerr << " w1p = " << ((void *)w1p) << " w2p = " << ((void *)w2p);
            // cerr << endl;
        }
#if 0
          hEw1Ew2->Fill(w1().E(),w2().E(),1.0);
          hMw1Mw2->Fill(w1mass,w2mass,1.0);
          hMt1Mt2->Fill(t1mass,t2mass,1.0);
#endif
        }
      }
      w2.UnlockChildren();
    }
    w1.UnlockChildren();
  }
  
  // Cut on No. of solutions.

  if ( !solutions.GetEntries() ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|m_jj - m_W| <= %g && |m_3j - m_t| <= %g",xM2j,xM3j);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on cos(theta_bW).
  
  nsol = 0;
    
  TIter nextsol(&solutions);			// iterator for solutions
  ANLPair *solp;
  while ((solp = (ANLPair *)nextsol())) {
    ANLPair &sol = *solp;
    ANLPair &bbww1 = *((ANLPair *)sol[0]);
    ANLPair &bbww2 = *((ANLPair *)sol[1]);
    ANL4DVector &bb1 = *(ANL4DVector *)bbww1[0]; // Dunno why, but had to
    ANL4DVector &ww1 = *(ANL4DVector *)bbww1[1]; // use new refs instead
    ANL4DVector &bb2  = *(ANL4DVector *)bbww2[0];// of b1, b2, w1, and w2
    ANL4DVector &ww2 = *(ANL4DVector *)bbww2[1]; // here in RCINT????
    Double_t cosbw1 = bb1.CosTheta(ww1);		 // Looks like the refs
    Double_t cosbw2 = bb2.CosTheta(ww2);		 // are in global scope??
    hCosbw1Cosbw2->Fill(cosbw1,cosbw2,1.0);
    if (cosbw1 > xCosbw || cosbw2 > xCosbw) {
      solutions.Remove(solp);			// discard this solution
      delete solp;
    }
    if (gDEBUG) {
            ANL4DVector *b1p = (ANL4DVector *)bbww1[0];
            ANLPair         *w1p = (ANLPair *)bbww1[1];
            ANL4DVector *b2p = (ANL4DVector *)bbww2[0];
            ANLPair         *w2p = (ANLPair *)bbww2[1];
            cerr << " Solution ---" << ++nsol << endl;
            cerr << " b1p = " << ((void *)b1p) << " b2p = " << ((void *)b2p);
            cerr << " w1p = " << ((void *)w1p) << " w2p = " << ((void *)w2p);
            cerr << endl;
    }
  }
  if ( !solutions.GetEntries() ) { tracks.Delete(); last->cd(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_bw)| > %g",xCosbw);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  if (gDEBUG/* || kTRUE */) {
    cerr << "------------------------------------------" << endl;
    cerr << "- w1candidates after :" << endl;
    w1candidates.DebugPrint();
  }
  
  // ----------------------
  // End of event selection
  // ----------------------

  if ( Ngoods == 0 ) {
    selid++;
    sprintf(msg,"END");
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  Ngoods++;
  
  cerr << "------------------------------------------" << endl
       << "Event " << jsf->GetEventNumber() 
       << ": Number of solutions = " << solutions.GetEntries() << endl
       << "------------------------------------------" << endl;

  // Sort the solutions in the ascending order of chi2 vlues.
  
  solutions.Sort();

  // Hists and plots for selected events.

  if (gDEBUG) {
    Int_t nj = 0;
    nextjet.Reset();
    while ((jetp = (ANLJet *)nextjet())) {
       cerr << "Jet " << ++nj << " : ";
       jetp->DebugPrint();
    }
  }

  hNsols->Fill(solutions.GetEntries());
  hEvisPl->Fill(fEvis,fPl,1.);
  
  nextsol.Reset();
  Int_t nsols = 0;
  while ((solp = (ANLPair *)nextsol())) {
    if ( nsols++ ) break;				// choose the best
    ANLPair &sol = *solp;
    ANLPair &bbbwww1 = *((ANLPair *)sol[0]);
    ANLPair &bbbwww2 = *((ANLPair *)sol[1]);
    ANL4DVector &www1 = *(ANL4DVector *)bbbwww1[1];
    ANL4DVector &www2 = *(ANL4DVector *)bbbwww2[1];
    Double_t chi2   = sol.GetQuality();
    Double_t w1mass = www1.GetMass();
    Double_t w2mass = www2.GetMass();
    Double_t t1mass = bbbwww1().GetMass();
    Double_t t2mass = bbbwww2().GetMass();
          hChi2->Fill(chi2);
#if 1
          hEw1Ew2->Fill(www1.E(),www2.E(),1.0);
          hMw1Mw2->Fill(w1mass,w2mass,1.0);
          hMt1Mt2->Fill(t1mass,t2mass,1.0);
#endif
  }

  // Clean up

  solutions.Delete();
  tracks.Delete();

  last->cd();
  return;
}

//_________________________________________________________
void UserTerminate()
{
  // This function is called at the end of job.
  cout << endl;
  cout << "  =============" << endl;
  cout << "   Cut Summary " << endl;
  cout << "  =============" << endl;
  cout << endl;
  cout << "  -----------------------------------------------------------" << endl;
  cout << "   ID   No.Events    Cut Description" << endl;
  cout << "  -----------------------------------------------------------" << endl;
  Int_t i;
  for ( i = 0; strncmp(&cutName[i][0],"END",4) && i < maxcut ; i++ ) {
    printf("  %3d  %10d  : %s\n",i,(int)hStat->GetBinContent(i+1),&cutName[i][0]);
  } 
  cout << "  -----------------------------------------------------------" << endl;
  return;
}

