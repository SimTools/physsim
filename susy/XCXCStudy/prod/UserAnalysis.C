//*************************************************************************
//* ================
//*  UserAnalysis.C
//* ================
//*
//* (Description)
//*    A very primitive sample script to study e+e- -> X+X-.
//* (Usage)
//*    	 $ jsf gui.C
//*    gui.C then invokes this script from within before event loop
//*    starts.
//* (Update Recored)
//*    1999/05/26  K.Fujii	Original version.
//*
//*************************************************************************
//
#include <iostream.h>
//*------------------------*//
//* User Analysis          *//
//*------------------------*//

Bool_t gDEBUG = kFALSE;

// Canvas //

TCanvas *cHist;

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
TH2F *hCosw1Cosw2;
TH2F *hMw1Mw2;
TH2F *hEvisPl;
TH1F *hAcop;

// Event Information //

Double_t fEcm = 500.0;	// CM energy (GeV)
Int_t    fNtracks;	// track multiplicity
Double_t fEvis;		// visible energy
Double_t fPt;		// Pt
Double_t fPl;		// Pl
Int_t    fNjets;	// jet multiplicity
Double_t fYcut;		// y_cut to force the event to 4 jets
Double_t fMjj[2];	// 1st and 2nd W candidate masses
Double_t fEjj[2];	// 1st and 2nd W candidate energies
Double_t fCosjj[2];	// cos(theta) of 1st and 2nd W candidates
Double_t fAcop;		// Acoplanarity for w1 and w2

// Event selection condition //

static const Int_t maxcut = 50;
Char_t cutName[maxcut][100];	// Cut names

Int_t    xNtracks  =     25; 	// No. of tracks
Double_t xEtrack   =   0.10;	// track energy
Double_t xEvis     =  80.00;	// Minimum visible energy
Double_t xPt       =  10.00;	// Pt minimum
Double_t xPl       = 999.00;	// Pl maximum
Double_t xYcut     =  0.004;	// y_cut to force the event to 4 jets
Int_t    xNjets    =      4;	// No. of jets
Double_t xEjet     =   5.00;	// E_jet minimum
Double_t xCosjet   =   0.98;	// |cos(theta_j)| maximum
Double_t xCosw     =   0.95;	// |cos(theta_w)| maximum
Double_t xM2j      =  12.00;	// |m_jj-m_W| maximum
Double_t xAcop     =  30.00;	// Acoplanarity minimum

// Some constants //

static const Double_t kMassW   = 80.00; // W mass
static const Double_t kMassZ   = 91.19; // Z mass
static const Double_t kSigmaMw =   4.0; // W mass resolution
static const Double_t kSigmaMz =   4.0; // W mass resolution
static const Int_t    kZoneX   =     4;	// No. X Zones in the Canvas
static const Int_t    kZoneY   =     4;	// No. Y Zones in the Canvas

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

  hStat       = new TH1F("hStat","Cut Statistics",  20,   0.0,  20.0);
  hNtracks    = new TH1F("hNtracks","No. tracks" ,  50,   0.0, 100.0);
  hEvis       = new TH1F("hEvis","Visible energy",  50,   0.0, 500.0);
  hPt         = new TH1F("hPt","Missing Pt"      ,  50,   0.0, 250.0);
  hNjets      = new TH1F("hNjets","No. jets"     ,  20,   0.0,  20.0);
  hEjet       = new TH1F("hEjet","Jet energy"    ,  50,   0.0, 100.0);
  hCosjet     = new TH1F("hCosjet","cos(theta_j)",  50,  -1.0,  +1.0);
  hNsols      = new TH1F("hNsols","No. solutions",  20,   0.0,  20.0);
  hChi2       = new TH1F("hChi2","Chi2"          ,  50,   0.0,  50.0);
  hEw1Ew2     = new TH2F("hEw1Ew2","(E_w1,E_w2)" ,  
  				  50,  0.0, 200.0,  50,   0.0, 200.0);
  hCosw1Cosw2 = new TH2F("hCosw1Cosw2","(cos_w1,cow_w2)",
				  50, -1.0,  +1.0,  50,  -1.0,  +1.0);
  hMw1Mw2     = new TH2F("hMw1Mw2","(m_w1,m_w2)" , 
  				  60, 50.0, 110.0,  60,  50.0, 110.0);
  hEvisPl     = new TH2F("hEvisPl","(Evis,Pl)"   , 
  				  60,  0.0, 600.0,  50,-100.0,+100.0);
  hAcop       = new TH1F("hAcop","Acoplanarity"  ,  90,   0.0, 180.0);

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
  } 
  else {
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
  cHist->cd(11);	hCosw1Cosw2->Draw();
  cHist->cd(12);	hMw1Mw2->Draw();
  cHist->cd(13);	hEvisPl->Draw();
  cHist->cd(14);	hAcop->Draw();

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
  // Debug a particular event.
  
  Int_t eventToDebug = -1;
  if (jsf->GetEventNumber() == eventToDebug) {
  	gDEBUG = kTRUE;
        cerr << "------------------------------------------" << endl;
        cerr << "Event " << jsf->GetEventNumber();
        cerr << endl;
  } else {
  	gDEBUG = kFALSE;
  }

  Char_t msg[60];

  // Analysis starts here.
  
  Float_t selid = -0.5;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) strcpy(&cutName[(Int_t)selid][0],"No cut");

  // Get event buffer and make combined tracks accessible.

  JSFSIMDST     *sds     = (JSFSIMDST*)jsf->FindModule("JSFSIMDST");
  JSFSIMDSTBuf  *evt     = (JSFSIMDSTBuf*)sds->EventBuf();
  Int_t          ntrks   = evt->GetNLTKCLTracks(); 	// No. of tracks 
  TClonesArray  *trks    = evt->GetLTKCLTracks(); 	// combined tracks
  
  ANL4DVector qsum;
  TObjArray tracks(500);
   
  // Select good tracks

  fNtracks = 0;
  for ( Int_t i = 0; i < ntrks; i++ ) {
    JSFLTKCLTrack *t = (JSFLTKCLTrack*)trks->UncheckedAt(i);
    if ( t->GetE() > xEtrack ) {
      ANL4DVector *qt = new ANL4DVector(t->GetPV());
      tracks.Add(qt); 		// track 4-momentum
      qsum += *qt;		// total 4-momentum
      fNtracks++;
    }				// *qt stays.
  }

  // Cut on No. of tracks.
  
  hNtracks->Fill(fNtracks);
  if ( fNtracks < xNtracks ) { tracks.Delete(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"N_tracks > %g",xNtracks);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  fEvis = qsum.E();		// E_vis
  fPt   = qsum.GetPt();		// P_t
  fPl   = qsum.Pz();		// P_l

  // Cut on Evis.

  hEvis->Fill(fEvis);
  if ( fEvis < xEvis ) { tracks.Delete(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"E_vis > %g",xEvis);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
 
  // Cut on Pt.

  hPt->Fill(fPt);
  if ( fPt < xPt ) { tracks.Delete(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Pt > %g",xPt);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
 
  // Cut on Pl.

  if ( TMath::Abs(fPl) > xPl ) { tracks.Delete(); return; }
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
  if ( fNjets < xNjets ) { tracks.Delete(); return; }
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
    
  if ( fNjets != xNjets ) { tracks.Delete(); return; }
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
  while ((jetp = (ANLJet *)nextjet())) {
    ANLJet &jet = *jetp;
    if (gDEBUG && kFALSE) jet.DebugPrint();
    Double_t ejet = jet().E();
    if (ejet < ejetmin) ejetmin = ejet;
    hEjet->Fill(ejet);			// Ejet
    Double_t cosj = jet.CosTheta();
    if (TMath::Abs(cosj) > TMath::Abs(cosjmax)) cosjmax = cosj;
    hCosjet->Fill(cosj);		// cos(theta_jet)
  }

  // Cut on Ejet_min.
  
  if ( ejetmin < xEjet ) { tracks.Delete(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Ejet > %g",xEjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on |cos(theta_j)|_max.
    
  if ( TMath::Abs(cosjmax) > xCosjet ) { tracks.Delete(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_j)| <= %g",xCosjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Find W candidates in given mass windows.
  // Avoid using indices since there might be empty slots.
  
  TObjArray solutions(10);
  ANLPairCombiner w1candidates(jets,jets);
  if (gDEBUG) {
    cerr << "------------------------------------------" << endl;
    cerr << "- w1candidates:" << endl;
    w1candidates.DebugPrint();
  }
  ANLPair *w1p, *w2p;
  while ((w1p = (ANLPair *)w1candidates())) {
    ANLPair &w1 = *w1p;
    Double_t w1mass = w1().GetMass();
    if (TMath::Abs(w1mass - kMassW) > xM2j) continue;	// w1 candidate found
    w1.LockChildren();					// now lock w1 daughters
    ANLPairCombiner w2candidates(w1candidates);		// w2 after w1
    while ((w2p = (ANLPair *)w2candidates())) {
      ANLPair &w2 = *w2p;
      if (w2.IsLocked()) continue;			// skip if locked
      Double_t w2mass = w2().GetMass();
      if (TMath::Abs(w2mass - kMassW) > xM2j) continue;	// w2 candidate found
      Double_t chi2 = TMath::Power((w1mass - kMassW)/kSigmaMw,2.)
                    + TMath::Power((w2mass - kMassW)/kSigmaMw,2.);
      solutions.Add(new ANLPair(w1p,w2p,chi2));
      // hMw1Mw2->Fill(w1mass,w2mass,1.0);
    }
    w1.UnlockChildren();
  }
  
  // Cut on No. of solutions.

  if ( !solutions.GetEntries() ) { tracks.Delete(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|m_jj - m_W| <= %g",xM2j);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  if (gDEBUG) {
    cerr << "------------------------------------------" << endl;
    cerr << "- w1candidates after 1:" << endl;
    w1candidates.DebugPrint();
  }

  // Cut on cos(theta_W).

  TIter nextsol(&solutions);
  ANLPair *sol;
  while ((sol = (ANLPair *)nextsol())) {
    ANL4DVector  &ww1 = *(ANL4DVector *)(*sol)[0];
    ANL4DVector  &ww2 = *(ANL4DVector *)(*sol)[1];
    Double_t ew1 = ww1.E();
    Double_t ew2 = ww2.E();
    hEw1Ew2->Fill(ew1,ew2,1.0);
    Double_t cosw1 = ww1.CosTheta();
    Double_t cosw2 = ww2.CosTheta();
    hCosw1Cosw2->Fill(cosw1,cosw2,1.0);
    if (TMath::Abs(cosw1) > xCosw || TMath::Abs(cosw2) > xCosw) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { tracks.Delete(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_w)| <= %g",xCosw);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  if (gDEBUG) {
    cerr << "------------------------------------------" << endl;
    cerr << "- w1candidates after 2:" << endl;
    w1candidates.DebugPrint();
  }
  
  // Cut on Acop.

  nextsol.Reset();
  while ((sol = (ANLPair *)nextsol())) {
    ANL4DVector  &www1 = *(ANL4DVector *)(*sol)[0];
    ANL4DVector  &www2 = *(ANL4DVector *)(*sol)[1];
    Double_t acop = www1.Acop(www2);
    hAcop->Fill(acop);
    if (acop < xAcop) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { tracks.Delete(); return; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Acop > %g",xAcop);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  if (gDEBUG) {
    cerr << "------------------------------------------" << endl;
    cerr << "- w1candidates after 3:" << endl;
    w1candidates.DebugPrint();
  }

  // ------------------------
  //  End of event selection
  // ------------------------

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

  if (gDEBUG && kFALSE) {
    Int_t nj = 0;
    nextjet.Reset();
    while ((jetp = (ANLJet *)nextjet())) {
       cerr << "------" << endl
            << "Jet " << ++nj << endl
            << "------" << endl;
       jetp->DebugPrint();
    }
  }

  hNsols->Fill(solutions.GetEntries());
  hEvisPl->Fill(fEvis,fPl,1.);
  
  nextsol.Reset();
  Int_t nsols = 0;
  while ((sol = (ANLPair *)nextsol())) {
    if ( nsols++ ) break;				// choose the best
    ANL4DVector  &wwww1 = *(ANL4DVector *)(*sol)[0];
    ANL4DVector  &wwww2 = *(ANL4DVector *)(*sol)[1];
    Double_t chi2   = sol->GetQuality();
    Double_t w1mass = wwww1.GetMass();
    Double_t w2mass = wwww2.GetMass();
          hChi2->Fill(chi2);
          hMw1Mw2->Fill(w1mass,w2mass,1.0);
  }
  
  // Clean up
  
  solutions.Delete();
  tracks.Delete();

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
