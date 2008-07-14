//*************************************************************************
//* ========================
//*  XCXC4JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC chargino pair data. 
//* (Requires)
//* 	library Anlib
//* 	library XCXCStudy
//* (Provides)
//* 	class XCXC4JAnalysis
//* 	class XCXC4JAnalysisBuf
//* (Usage)
//*   Take a look at Anl.C.
//* (Update Recored)
//*   1999/08/01  K.Fujii	Original version.
//*
//*************************************************************************
//
#include "XCXC4JAnalysis.h"
#include "ANLTrack.h"

static const Double_t kMassW   = 80.00; // W mass
static const Double_t kMassZ   = 91.19; // Z mass
static const Double_t kSigmaMw =   4.0; // W mass resolution
static const Double_t kSigmaMz =   4.0; // W mass resolution
static const Int_t    kZoneX   = 4;	// No. X Zones in the Canvas
static const Int_t    kZoneY   = 4;	// No. Y Zones in the Canvas

Int_t XCXC4JAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

typedef enum { kElectron = 11, kMuon = 13 } EPID;

//_____________________________________________________________________
//  -----------------------
//  XCXC4JAnalysisBuf Class
//  -----------------------
//
//
ClassImp(XCXC4JAnalysisBuf)

//_________________________________________________________
XCXC4JAnalysisBuf::XCXC4JAnalysisBuf(const Char_t *name, const Char_t *title, 
   XCXC4JAnalysis *module) : JSFEventBuf(name, title, (JSFModule*)module) {}

//_________________________________________________________
XCXC4JAnalysisBuf::XCXC4JAnalysisBuf(XCXC4JAnalysis *module, const Char_t *name,
   const Char_t *title) : JSFEventBuf(name, title, (JSFModule*)module) {}


//_____________________________________________________________________
//  --------------------
//  XCXC4JAnalysis Class
//  --------------------
//
//

ClassImp(XCXC4JAnalysis)

XCXC4JAnalysis::XCXC4JAnalysis(const Char_t *name, const Char_t *title)
	       : JSFModule(name, title) 
{
  fEventBuf = new XCXC4JAnalysisBuf(this);
  SetBufferSize(2000);  // buffer size for event data.
  cout << "XCXC4JAnalysisBuf is created...fEventBuf is " 
       << (Int_t)fEventBuf << endl;
}

//_____________________________________________________________________
XCXC4JAnalysis::~XCXC4JAnalysis()
{
  cout << "XCXC4JAnalysisBuf will be deleted...fEventBuf is " 
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}

//_____________________________________________________________________
void XCXC4JAnalysis::CleanUp(TObjArray *objs)
{
  TIter next(objs);
  TObject *obj;
  while ((obj = next())) {
  	objs->Remove(obj);
  	delete obj;
  }
}

//_____________________________________________________________________
Bool_t XCXC4JAnalysis::Initialize()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat       = new TH1F("hStat","Cut Statistics",  20,   0.0,  20.0);
  hNtracks    = new TH1F("hNtracks","No. tracks" ,  50,   0.0, 100.0);
  hEvis       = new TH1F("hEvis","Visible energy", 400,   0.0, 400.0);
  hPt         = new TH1F("hPt","Missing Pt"      ,  50,   0.0, 250.0);
  hNjets      = new TH1F("hNjets","No. jets"     ,  20,   0.0,  20.0);
  hEjet       = new TH1F("hEjet","Jet energy"    ,  50,   0.0, 100.0);
  hElep       = new TH1F("hElep","Lepton energy" ,  50,   0.0, 100.0);
  hCosjet     = new TH1F("hCosjet","cos(theta_j)",  50,  -1.0,  +1.0);
  hNsols      = new TH1F("hNsols","No. solutions",  20,   0.0,  20.0);
  hChi2       = new TH1F("hChi2","Chi2"          ,  50,   0.0,  50.0);
  hEw1Ew2     = new TH2F("hEw1Ew2","(E_w1,E_w2)" ,  
  				 200,  0.0, 200.0, 200,   0.0, 200.0);
  hCosw1Cosw2 = new TH2F("hCosw1Cosw2","(cos_w1,cow_w2)",
				  50, -1.0,  +1.0,  50,  -1.0,  +1.0);
  hMw1Mw2     = new TH2F("hMw1Mw2","(m_w1,m_w2)" , 
  				  60, 50.0, 110.0,  60,  50.0, 110.0);
  hEvisPl     = new TH2F("hEvisPl","(Evis,Pl)"   , 
  				  60,  0.0, 600.0,  50,-100.0,+100.0);
  hMM         = new TH1F("hMM","mm_ww"           ,  80, 100.0, 500.0);
  hAcop       = new TH1F("hAcop","Acoplanarity"  ,  90,   0.0, 180.0);
  hEw         = new TH1F("hEw","E_w"             , 400,   0.0, 200.0);

  xNtracks  =      25;   // No. of tracks
  xEtrack   =    0.10;   // track energy
  xEvisLo   =   20.00;   // Minimum visible energy
  xEvisHi   =  400.00;   // Maximum visible energy
  xPt       =    0.00;   // Pt minimum
  xPl       = 9999.00;   // Pl maximum
  xEl       =   25.00;   // Maximum lepton energy
  xYcut     =   0.005;   // y_cut to force the event to 4 jets
  xNjets    =       4;   // No. of jets
  xEjet     =    5.00;	// E_jet minimum
  xCosjet1  =    0.80;	// |cos(theta_j)| maximum
  xCosjet2  =    0.95;	// |cos(theta_j)| maximum
  xCosw     =    0.90;	// |cos(theta_j)| maximum
  xM2jLo    =   10.00;	// |m_jj-m_W| maximum
  xM2jHi    =   20.00;	// |m_jj-m_W| maximum
  xMM1      =   70.00;   // missing mass cut against WW
  xMM2      =  120.00;   // missing mass cut against WW
  xAcop     =   30.00;	// Acoplanarity minimum

  last->cd();
  return 0;
}

//_________________________________________________________
void XCXC4JAnalysis::DrawHist()
{
  TDirectory *last = gDirectory;
  if (!cHist) {
    cHist = new TCanvas("cHist","Canvas 1",10, 10, kZoneX*200, kZoneY*200);
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
  cHist->cd(11);	hCosw1Cosw2->Draw();
  cHist->cd(12);	hMw1Mw2->Draw();
  cHist->cd(13);	hEvisPl->Draw();
  cHist->cd(14);	hAcop->Draw();

  cHist->Update();
  
  last->cd();
}

//_________________________________________________________
Bool_t XCXC4JAnalysis::Process(Int_t ev)
{
  // Local copies of XCXC4JAnalysisBuf data members.

  Int_t     	fNtracks;	// track multiplicity
  Double_t  	fEvis;		// visible energy
  Double_t  	fPt;		// Pt
  Double_t  	fPl;		// Pl
  Double_t  	fElmax = 0.;	// Elmax
  Double_t  	fYcut;		// y_cut to force the event to 4 jets
  Int_t        	fNjets;		// jet multiplicity
  Double_t      fEcm;		// Ecm
  Double_t      fMM;		// missing mass
 
  // Remember the previous directory.
  
  TDirectory *last = gDirectory;
  gFile->cd("/");

  Char_t msg[60];

  // Analysis starts here.
  
  Float_t selid = -0.5;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) strcpy(&cutName[(Int_t)selid][0],"No cut");

  // Get event buffer and make combined tracks accessible.

  JSFSIMDST     *sds     = (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
  JSFSIMDSTBuf  *evt     = (JSFSIMDSTBuf*)sds->EventBuf();
  XCXC4JAnalysisBuf *ua  = (XCXC4JAnalysisBuf *)fEventBuf;
  XCXC4JAnalysisBuf &a   = *ua;

  Int_t          ntrks   = evt->GetNLTKCLTracks(); 	// No. of tracks 
  TObjArray     *trks    = evt->GetLTKCLTracks(); 	// combined tracks

  fEcm = evt->GetEcm();

  // Select good tracks

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
      if (t->GetType() == kMuon || t->GetType() == kElectron) {
        if (t->GetE() > fElmax) fElmax = t->GetE();
      }
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

  fEvis = qsum(0);		// E_vis
  fPt   = qsum.GetPt();	// P_t
  fPl   = qsum(3);		// P_l
  
  if (gDEBUG) cerr << "Evis = " << fEvis << " Pt = " 
       << fPt << " Pl = " << fPl << endl;

  // Cut on Evis.

  hEvis->Fill(fEvis);
  if (fEvis < xEvisLo || fEvis > xEvisHi) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"%g GeV < E_vis < %g GeV", xEvisLo, xEvisHi);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
 
  // Cut on Pt.

  hPt->Fill(fPt);
  if ( fPt < xPt ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Pt > %g GeV",xPt);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
 
  // Cut on Pl.

  if ( TMath::Abs(fPl) > xPl ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|Pl| <= %g GeV",xPl);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Elepton.

  hElep->Fill(fElmax);
  if (fElmax > xEl) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Elep <= %g GeV",xEl);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Find jets.

  fYcut = xYcut;
  ANLJadeEJetFinder jclust(fYcut);
  jclust.Initialize(tracks);
  jclust.FindJets();
  fYcut  = jclust.GetYcut();
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
  fYcut  = jclust.GetYcut();

  if (gDEBUG) cerr << "Ycut = " << fYcut << " Njets = " << fNjets << endl;

  // Make sure that No. of jets is xNjets.
    
  if ( fNjets != xNjets ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Njets = %i",xNjets);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  TObjArray &jets = jclust.GetJets();
  TIter nextjet(&jets);
  ANLJet *jetp;
  Double_t ejetmin = 999999.;
  Double_t cosjmax = 0.;
  Int_t ncos = 0;
  while ((jetp = (ANLJet *)nextjet())) {
    ANLJet &jet  = *jetp;
    if (gDEBUG) jet.DebugPrint();
    Double_t ejet = jet()(0);
    if (ejet < ejetmin) ejetmin = ejet;
    hEjet->Fill(ejet);
    Double_t cosj = jet.CosTheta();
    if (TMath::Abs(cosj) > TMath::Abs(cosjmax)) cosjmax = cosj;
    if (TMath::Abs(cosj) < xCosjet1) ncos++;
    hCosjet->Fill(cosj);
  }

  // Cut on Ejet_min.
  
  if ( ejetmin < xEjet ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Ejet > %g GeV",xEjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on |cos(theta_j)|.
    
  if (ncos < 2) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"# jets with |cos(theta_j)| < %g >= 2",xCosjet1);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on |cos(theta_j)|_max.
    
  if (TMath::Abs(cosjmax) > xCosjet2) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_j)| <= %g",xCosjet2);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Find W candidates in a given mass window.
  
  TObjArray solutions(10);
  ANLPairCombiner w1candidates(jets,jets);
  ANLPair *w1p, *w2p;
  while ((w1p = (ANLPair *)w1candidates())) {
    ANLPair &w1 = *w1p;
    Double_t w1mass = w1().GetMass();
    if (w1mass < kMassW-xM2jLo || w1mass > kMassW+xM2jHi) continue;
    w1.LockChildren();
    ANLPairCombiner w2candidates(w1candidates);
    while ((w2p = (ANLPair *)w2candidates())) {
      ANLPair &w2 = *w2p;
      if (w2.IsLocked()) continue;
      Double_t w2mass = w2().GetMass();
      if (w2mass < kMassW-xM2jLo || w2mass > kMassW+xM2jHi) continue;
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
      // hMw1Mw2->Fill(w1mass,w2mass,1.0);
    }
    w1.UnlockChildren();
  }
  
  // Cut on No. of solutions.

  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"m_W - %g GeV < m_jj < m_W + %g GeV", xM2jLo, xM2jHi);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  // Cut on cos(theta_W).

  TIter nextsol(&solutions);
  ANLPair *sol;
  while ((sol = (ANLPair *)nextsol())) {
    ANL4DVector  &w1 = *(ANL4DVector *)(*sol)[0];
    ANL4DVector  &w2 = *(ANL4DVector *)(*sol)[1];
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

  // Cut on missing mass.

  ANL4DVector qcm(fEcm);
  ANL4DVector qmm = qcm -qsum;
  fMM = qmm.GetMass();
  hMM->Fill(fMM,1.);

  if (fMM > xMM1 && fMM < xMM2) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"mm_WW <= %g GeV or mm_WW >= %g GeV",xMM1, xMM2);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  // Cut on Acop.

  nextsol.Reset();
  while ((sol = (ANLPair *)nextsol())) {
    ANL4DVector  &w1 = *(ANL4DVector *)(*sol)[0];
    ANL4DVector  &w2 = *(ANL4DVector *)(*sol)[1];
    Double_t acop = w1.Acop(w2);
    hAcop->Fill(acop);
    if (acop < xAcop) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Acop > %g deg.",xAcop);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  // End of event selection

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

  // Now store this in XCXC4JAnalysisBuf.
  
  a.fNtracks	= fNtracks;
  a.fEvis	= fEvis;
  a.fPt		= fPt;
  a.fPl		= fPl;
  a.fYcut	= fYcut;
  a.fNjets	= fNjets;

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
#if 0
    Double_t ew1 = w1()(0);
    Double_t ew2 = w2()(0);
#else
    Double_t ew1 = TMath::Sqrt(kMassW*kMassW + w1().GetMag2());
    Double_t ew2 = TMath::Sqrt(kMassW*kMassW + w2().GetMag2());
#endif
          hChi2->Fill(chi2);
          hMw1Mw2->Fill(w1mass,w2mass,1.0);
          hEw1Ew2->Fill(ew1,ew2,1.0);
          hEw->Fill(ew1,1.0);
          hEw->Fill(ew2,1.0);
  }
  
  // Clean up
  
  CleanUp(&solutions);
  CleanUp(&tracks);
  
  last->cd();
  return kTRUE;
}

//_________________________________________________________
Bool_t XCXC4JAnalysis::Terminate()
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
