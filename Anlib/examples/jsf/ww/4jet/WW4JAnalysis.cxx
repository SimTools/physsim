//*************************************************************************
//* ======================
//*  WW4JAnalysis Classes
//* ======================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC W boson pair data. 
//* (Requires)
//*     library Anlib
//*     library WWStudy
//* (Provides)
//*     class WW4JAnalysis
//*     class WW4JAnalysisBuf
//* (Usage)
//*   Take a look at anl4J.C.
//* (Update Recored)
//*   1999/08/01  K.Fujii       Original version.
//*   2000/05/03  K.Ikematsu    Minor update.
//*
//*************************************************************************
//
#include "WW4JAnalysis.h"
#include "ANLTrack.h"

static const Double_t kMassW   = 80.00; // W mass
static const Double_t kMassZ   = 91.19; // Z mass
static const Double_t kSigmaMw =   4.0; // W mass resolution
static const Double_t kSigmaMz =   4.0; // W mass resolution
static const Int_t    kZoneX   = 4;     // No. X Zones in the Canvas
static const Int_t    kZoneY   = 4;     // No. Y Zones in the Canvas

Int_t WW4JAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

//_____________________________________________________________________
//  ---------------------
//  WW4JAnalysisBuf Class
//  ---------------------
//
//
ClassImp(WW4JAnalysisBuf)

//_________________________________________________________
WW4JAnalysisBuf::WW4JAnalysisBuf(const Char_t *name, const Char_t *title, 
   WW4JAnalysis *module) : JSFEventBuf(name, title, (JSFModule*)module) {}

//_________________________________________________________
WW4JAnalysisBuf::WW4JAnalysisBuf(WW4JAnalysis *module, const Char_t *name,
   const Char_t *title) : JSFEventBuf(name, title, (JSFModule*)module) {}

//_____________________________________________________________________
//  ------------------
//  WW4JAnalysis Class
//  ------------------
//
//

ClassImp(WW4JAnalysis)

WW4JAnalysis::WW4JAnalysis(const Char_t *name, const Char_t *title)
	       : JSFModule(name, title) 
{
  fEventBuf = new WW4JAnalysisBuf(this);
  SetBufferSize(2000);  // buffer size for event data.
  cout << "WW4JAnalysisBuf is created...fEventBuf is " 
       << (Int_t)fEventBuf << endl;
}

//_____________________________________________________________________
WW4JAnalysis::~WW4JAnalysis()
{
  cout << "WW4JAnalysisBuf will be deleted...fEventBuf is " 
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}

//_____________________________________________________________________
void WW4JAnalysis::CleanUp(TObjArray *objs)
{
#if 0
  TIter next(objs);
  TObject *obj;
  while ((obj = next())) {
  	objs->Remove(obj);
  	delete obj;
  }
#else
  objs->SetOwner();
#endif
}

//_____________________________________________________________________
Bool_t WW4JAnalysis::Initialize()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat       = new TH1F("hStat","Cut Statistics",  20,   0.0,  20.0);
  hNtracks    = new TH1F("hNtracks","No. tracks" ,  50,   0.0, 100.0);
  hEvis       = new TH1F("hEvis","Visible energy",  50,   0.0, 500.0);
  hPt         = new TH1F("hPt","Missing Pt"      ,  50,   0.0, 250.0);
  hNjets      = new TH1F("hNjets","No. jets"     ,  20,   0.0,  20.0);
  hYcut       = new TH1F("hYcut","Ymax"           ,100,   0.0,   1.0);
  hEjet       = new TH1F("hEjet","Jet energy"    ,  50,   0.0, 200.0);
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

  xNtracks  =     25;   // No. of tracks
  xEtrack   =   0.10;   // track energy
  xEvis     =  80.00;   // Minimum visible energy
  xPt       =  10.00;   // Pt maximum
  xPl       = 999.00;   // Pl maximum
  xYcut     =  0.004;   // y_cut to force the event to 4 jets
  xNjets    =      4;   // No. of jets
  xEjet     =   5.00;	// E_jet minimum
  xCosjet   =   0.99;	// |cos(theta_j)| maximum
  xCosw     =   0.97;	// |cos(theta_j)| maximum
  xM2j      =  18.00;	// |m_jj-m_W| maximum
  xAcop     =  30.00;	// Acoplanarity maximum

  last->cd();
  return 0;
}

//_________________________________________________________
void WW4JAnalysis::DrawHist()
{
  TDirectory *last = gDirectory;
  if (!cHist) {
    cHist = new TCanvas("cHist","Canvas 1",10, 10, kZoneX*200, kZoneY*200);
    cHist->Divide(kZoneX,kZoneY);
  } else {
    cHist->cd();
  }

  Int_t Ihist = 0;
  cHist->cd(++Ihist);    hStat->Draw();
  cHist->cd(++Ihist);    hNtracks->Draw();
  cHist->cd(++Ihist);    hEvis->Draw();
  cHist->cd(++Ihist);    hPt->Draw();
  cHist->cd(++Ihist);    hNjets->Draw();
  cHist->cd(++Ihist);    hYcut->Draw();
  cHist->cd(++Ihist);    hEjet->Draw();
  cHist->cd(++Ihist);    hCosjet->Draw();
  cHist->cd(++Ihist);    hNsols->Draw();
  cHist->cd(++Ihist);    hChi2->Draw();
  cHist->cd(++Ihist);    hEw1Ew2->Draw();
  cHist->cd(++Ihist);    hCosw1Cosw2->Draw();
  cHist->cd(++Ihist);    hMw1Mw2->Draw();
  cHist->cd(++Ihist);    hEvisPl->Draw();
  cHist->cd(++Ihist);    hAcop->Draw();

  cHist->Update();
  
  last->cd();
}

//_________________________________________________________
Bool_t WW4JAnalysis::Process(Int_t ev)
{
  // Local copies of WW4JAnalysisBuf data members.

  Int_t     	fNtracks;	// track multiplicity
  Double_t  	fEvis;		// visible energy
  Double_t  	fPt;		// Pt
  Double_t  	fPl;		// Pl
  Double_t  	fYcut;		// y_cut to force the event to 4 jets
  Int_t        	fNjets;		// jet multiplicity
  Double_t      fAcop;          // Acoplanarity

  // Remember the previous directory.

  TDirectory *last = gDirectory;
  gFile->cd("/");

  Char_t msg[60];

  // ---------------------
  // Analysis starts here.
  // ---------------------

  if (gDEBUG) {
    cerr << "------------------------------------------" << endl
         << "Event " << gJSF->GetEventNumber() 	       << endl
         << "------------------------------------------" << endl;
  }

  Float_t selid = -0.5;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) strcpy(&cutName[(Int_t)selid][0],"No cut");

  // Get event buffer and make combined tracks accessible.

  JSFSIMDST     *sds     = (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
  JSFSIMDSTBuf  *evt     = (JSFSIMDSTBuf*)sds->EventBuf();
  WW4JAnalysisBuf *ua    = (WW4JAnalysisBuf *)fEventBuf;
  WW4JAnalysisBuf &a     = *ua;

  Int_t          ntrks   = evt->GetNLTKCLTracks(); 	// No. of tracks 
  TObjArray     *trks    = evt->GetLTKCLTracks(); 	// combined tracks

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
  if (gDEBUG) cerr << "Ntracks = " << fNtracks << endl;

  // Cut on No. of tracks.

  hNtracks->Fill(fNtracks);
  if ( fNtracks < xNtracks ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"N_tracks > %i",xNtracks);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  fEvis = qsum(0);              // E_vis
  fPt   = qsum.GetPt();         // P_t
  fPl   = qsum(3);              // P_l

  if (gDEBUG) cerr << "Evis = " << fEvis << " Pt = " 
       << fPt << " Pl = " << fPl << endl;

  // Cut on Evis.

  hEvis->Fill(fEvis);
  if ( fEvis < xEvis ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"E_vis > %g",xEvis);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Pt.

  hPt->Fill(fPt);
  if ( fPt > xPt ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Pt <= %g",xPt);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Pl.

  if ( TMath::Abs(fPl) > xPl ) { CleanUp(&tracks); return kFALSE; }
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
  fYcut  = jclust.GetYmax();
  fNjets = jclust.GetNjets();

  if (gDEBUG) cerr << "Ycut = " << fYcut << " Njets = " << fNjets << endl;

  // Cut on No. of jets.
    
  hNjets->Fill(fNjets);
  if ( fNjets < xNjets ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Njets >= %i for Ycut = %g",xNjets,xYcut);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Now force the event to be xNjets.

  jclust.ForceNJets(xNjets);
  fNjets = jclust.GetNjets();
  fYcut  = jclust.GetYmax();
  hYcut->Fill(fYcut);

  if (gDEBUG) cerr << "Ycut = " << fYcut << " Njets = " << fNjets << endl;

  // Make sure that No. of jets is xNjets.

  if ( fNjets != xNjets ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Njets = %i",xNjets);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Loop over jets and decide Ejet_min and |cos(theta_j)|_max.

  TObjArray &jets = jclust.GetJets();
  TIter nextjet(&jets);
  ANLJet *jetp;
  Double_t ejetmin = 999999.;
  Double_t cosjmax = 0.;
  while ((jetp = (ANLJet *)nextjet())) {
    ANLJet &jet = *jetp;
    if (gDEBUG) jet.DebugPrint();
    Double_t ejet = jet()(0);
    if (ejet < ejetmin) ejetmin = ejet;
    hEjet->Fill(ejet);
    Double_t cosj = jet.CosTheta();
    if (TMath::Abs(cosj) > TMath::Abs(cosjmax)) cosjmax = cosj;
    hCosjet->Fill(cosj);
  }

  // Cut on Ejet_min.

  if ( ejetmin < xEjet ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Ejet > %g",xEjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on |cos(theta_j)|_max.

  if ( TMath::Abs(cosjmax) > xCosjet ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_j)| <= %g",xCosjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Find W candidates in a given mass window.

  TObjArray solutions(10);
  ANLPairCombiner w1candidates(jets,jets);
  ANLPair *w1p, *w2p;
  while ((w1p = (ANLPair *)w1candidates())) {
    ANLPair &w1 = *w1p;
    Double_t w1mass = w1().GetMass();
    if (TMath::Abs(w1mass - kMassW) > xM2j) continue;
    w1.LockChildren();
    ANLPairCombiner w2candidates(w1candidates);
    while ((w2p = (ANLPair *)w2candidates())) {
      ANLPair &w2 = *w2p;
      if (w2.IsLocked()) continue;
      Double_t w2mass = w2().GetMass();
      if (TMath::Abs(w2mass - kMassW) > xM2j) continue;
      if (gDEBUG) { 
        cerr << " M_w1 = " << w1mass << " M_w2 = " << w2mass << endl;
        cerr << " w1p  = " << (void *)w1p 
             << " w2p  = " << (void *)w2p << endl;
        cerr << " w1[0] = " << (void *)w1[0] 
             << " w1[1] = " << (void *)w1[1]
             << " w2[0] = " << (void *)w2[0] 
             << " w2[1] = " << (void *)w2[1] << endl;
      }
      Double_t chi2 = TMath::Power((w1mass - kMassW)/kSigmaMw,2.)
                    + TMath::Power((w2mass - kMassW)/kSigmaMw,2.);
      solutions.Add(new ANLPair(w1p,w2p,chi2));
    }
    w1.UnlockChildren();
  }

  // Cut on No. of solutions.

  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|m_jj - m_W| <= %g",xM2j);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on cos(theta_W).

  TIter nextsol(&solutions);
  ANLPair *sol;
  while ((sol = (ANLPair *)nextsol())) {
    ANL4DVector  &w1 = *(ANL4DVector *)(*sol)[0];
    ANL4DVector  &w2 = *(ANL4DVector *)(*sol)[1];
    Double_t ew1 = w1(0);
    Double_t ew2 = w2(0);
    hEw1Ew2->Fill(ew1,ew2,1.0);
    Double_t cosw1 = w1.CosTheta();
    Double_t cosw2 = w2.CosTheta();
    hCosw1Cosw2->Fill(cosw1,cosw2,1.0);
    if (TMath::Abs(cosw1) > xCosw || TMath::Abs(cosw2) > xCosw) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_w)| <= %g",xCosw);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Acop.

  fAcop = 0.;
  nextsol.Reset();
  while ((sol = (ANLPair *)nextsol())) {
    ANL4DVector  &w1 = *(ANL4DVector *)(*sol)[0];
    ANL4DVector  &w2 = *(ANL4DVector *)(*sol)[1];
    fAcop = w1.Acop(w2);
    hAcop->Fill(fAcop);
    if (fAcop > xAcop) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Acop <= %g",xAcop);
    strcpy(&cutName[(Int_t)selid][0],msg);
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
       << "Event " << gJSF->GetEventNumber()
       << ": Number of solutions = " << solutions.GetEntries() << endl
       << "------------------------------------------" << endl;

  // Sort the solutions in the ascending order of chi2 vlues.

  solutions.Sort();

  // Now store this in WW4JAnalysisBuf.

  a.fNtracks	= fNtracks;
  a.fEvis	= fEvis;
  a.fPt		= fPt;
  a.fPl		= fPl;
  a.fYcut	= fYcut;
  a.fNjets	= fNjets;
  a.fAcop	= fAcop;

  // Hists and plots for selected events.

  if (gDEBUG) {
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
    ANLPair &w1 = *(ANLPair *)(*sol)[0];
    ANLPair &w2 = *(ANLPair *)(*sol)[1];
    Double_t chi2   = sol->GetQuality();
    Double_t w1mass = w1.GetMass();
    Double_t w2mass = w2.GetMass();
          hChi2->Fill(chi2);
          hMw1Mw2->Fill(w1mass,w2mass,1.0);
  }

  // Clean up

  nextsol.Reset();
#if 0
// Bug:
// we should not delete w's contained in sol since they
// belong to ANLPairCombiner w1candidates which will delete
// them in its dtor.
//
  while ((sol = (ANLPair *)nextsol())) sol->Delete();
#endif
  CleanUp(&solutions);
  CleanUp(&tracks);

  last->cd();
  return kTRUE;
}

//_________________________________________________________
Bool_t WW4JAnalysis::Terminate()
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
  for ( i = 0; strncmp(&cutName[i][0],"END",4) && i < MAXCUT ; i++ ) {
    printf("  %3d  %10d  : %s\n",i,(int)hStat->GetBinContent(i+1),&cutName[i][0]);
  } 
  cout << "  -----------------------------------------------------------" << endl;
  return 0;
}
