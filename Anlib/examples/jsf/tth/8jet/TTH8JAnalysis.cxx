//*************************************************************************
//* =======================
//*  TTH8JAnalysis Classes
//* =======================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC ttbar data to select 8-jet events.
//* (Requires)
//* 	library Anlib
//* 	library TTHStudy
//* (Provides)
//* 	class TTH8JAnalysis
//* 	class TTH8JAnalysisBuf
//* (Usage)
//*   Take a look at anl8J.C.
//* (Update Recored)
//*   2002/08/15  K.Fujii	Original version.
//*************************************************************************
//
#include "TTH8JAnalysis.h"
#include "TTHSpring.h"

static const Double_t kMassW   = 80.00; 	// W mass
static const Double_t kMassZ   = 91.19; 	// Z mass
static const Double_t kMasst   = 175.0; 	// top mass
static const Double_t kMassH   = 120.0;		// Higgs mass
static const Double_t kSigmaMw =   8.0; 	// W mass resolution
static const Double_t kSigmaMz =   8.0; 	// Z mass resolution
static const Double_t kSigmaMt =  15.0; 	// top mass resolution
static const Double_t kSigmaMh =   8.0; 	// H mass resolution
static const Int_t    kZoneX   = 4;		// No. X Zones in the Canvas
static const Int_t    kZoneY   = 4;		// No. Y Zones in the Canvas

Int_t TTH8JAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

//_____________________________________________________________________
//  ----------------------
//  TTH8JAnalysisBuf Class
//  ----------------------
//
//
ClassImp(TTH8JAnalysisBuf)

//_________________________________________________________
TTH8JAnalysisBuf::TTH8JAnalysisBuf(const Char_t *name, const Char_t *title,
   TTH8JAnalysis *module) : JSFEventBuf(name, title, (JSFModule*)module) {}

//_________________________________________________________
TTH8JAnalysisBuf::TTH8JAnalysisBuf(TTH8JAnalysis *module, const Char_t *name,
   const Char_t *title) : JSFEventBuf(name, title, (JSFModule*)module) {}


//_____________________________________________________________________
//  ------------------
//  TTH8JAnalysis Class
//  ------------------
//
//

ClassImp(TTH8JAnalysis)

TTH8JAnalysis::TTH8JAnalysis(const Char_t *name, const Char_t *title)
  : JSFModule(name, title), cHist(0)
{
  fEventBuf = new TTH8JAnalysisBuf(this);
  SetBufferSize(2000);  // buffer size for event data.
  cout << "TTH8JAnalysisBuf is created...fEventBuf is "
       << (Int_t)fEventBuf << endl;
}

//_____________________________________________________________________
TTH8JAnalysis::~TTH8JAnalysis()
{
  cout << "TTH8JAnalysisBuf will be deleted...fEventBuf is "
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}

//_____________________________________________________________________
void TTH8JAnalysis::CleanUp(TObjArray *objs)
{
  objs->SetOwner();
}

//_____________________________________________________________________
Bool_t TTH8JAnalysis::Initialize()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat         = new TH1F("hStat","Cut Statistics",  20,   0.0,  20.0);
  hNtracks      = new TH1F("hNtracks","No. tracks" ,  50,   0.0, 400.0);
  hEvis         = new TH1F("hEvis","Visible energy",  80,   0.0, 800.0);
  hPt           = new TH1F("hPt","Missing Pt"      ,  60,   0.0, 240.0);
  hNjets        = new TH1F("hNjets","No. jets"     ,  20,   0.0,  20.0);
  hEjet         = new TH1F("hEjet","Jet energy"    ,  50,   0.0, 200.0);
  hCosjet       = new TH1F("hCosjet","cos(theta_j)",  50,  -1.0,  +1.0);
  hNsols        = new TH1F("hNsols","No. solutions",  20,   0.0,  20.0);
  hChi2         = new TH1F("hChi2","Chi2"          ,  50,   0.0,  50.0);
  hEw1Ew2       = new TH2F("hEw1Ew2","(E_w1,E_w2)" ,
  				    50,  0.0, 400.0,  50,   0.0, 400.0);
  hCosbw1Cosbw2 = new TH2F("hCosbw1Cosbw2","(cos_bw1,cow_bw2)",
				    50, -1.0,  +1.0,  50,  -1.0,  +1.0);
  hMw1Mw2       = new TH2F("hMw1Mw2","(m_w1,m_w2)" ,
  				    60, 50.0, 110.0,  60,  50.0, 110.0);
  hEvisPl       = new TH2F("hEvisPl","(Evis,Pl)"   ,
  				    80,  0.0, 800.0,  50,-200.0,+200.0);
  hMt1Mt2       = new TH2F("hMt1Mt2","(m_t1,m_t2)" ,
  				    50,120.0, 220.0,  50, 120.0, 220.0);
  hMw2Mh        = new TH2F("hMw2MH","(m_w2,m_H)" ,
  				    60, 50.0, 110.0,  50,  70.0, 170.0);
  hThrust       = new TH1F("hThrust","Thrust"      ,  50,   0.0,   1.0);
  hYcut         = new TH1F("hYcut","Ymax"          ,1000,   0.0,  0.02);

  xNtracks  =     25;   // No. of tracks
  xEtrack   =   0.10;   // track energy
  xEvis     = 500.00;   // Minimum visible energy
  xPt       =  50.00;   // Pt maximum
  xPl       = 999.00;   // Pl maximum
  xYcut     =  0.004;   // y_cut to force the event to 4 jets
  xNjets    =      8;   // No. of jets
  xEjet     =   5.00;	// E_jet minimum
  xCosjet   =   0.99;	// |cos(theta_j)| maximum
  xCosbw    =    1.0;	// cos(theta_bw) maximum
  xM2j      =  18.00;	// |m_jj-m_W| maximum
  xM3j      =  24.00;	// |m_3j-m_t| maximum
  xThrust   =   0.80;   // Thrust maximum

  for (Int_t i = 0; i < MAXCUT; i++) {
    strcpy(&cutName[i][0],"     ");
  }

   //--
   //  Read in Generator info.
   //--
   gJSF->GetInput()->cd("/conf/init");
   TTHBases *bsp = (TTHBases *)gROOT->FindObject("TTHBases");
   cerr << "------------------------------------" << endl
        << " Ecm = " << bsp->GetRoots() << " GeV" << endl
        << "------------------------------------" << endl;
    TTH8JAnalysisBuf *bufp = (TTH8JAnalysisBuf *)EventBuf();
    bufp->SetEcm(bsp->GetRoots()); 

  last->cd();
  return 0;
}

//_________________________________________________________
void TTH8JAnalysis::DrawHist()
{
  TDirectory *last = gDirectory;
  if (!cHist) {
    cHist = new TCanvas("cHist","Canvas 1",10, 10, kZoneX*200, kZoneY*200);
    cHist->Divide(kZoneX,kZoneY);
  } else {
    cHist->cd();
  }

  Int_t Ihist = 0;
  cHist->cd(++Ihist);		hStat->Draw();
  cHist->cd(++Ihist);		hNtracks->Draw();
  cHist->cd(++Ihist);		hEvis->Draw();
  cHist->cd(++Ihist);		hPt->Draw();
  cHist->cd(++Ihist);		hNjets->Draw();
  cHist->cd(++Ihist);		hEjet->Draw();
  cHist->cd(++Ihist);		hCosjet->Draw();
  cHist->cd(++Ihist);		hNsols->Draw();
  cHist->cd(++Ihist);		hChi2->Draw();
  cHist->cd(++Ihist);		hEw1Ew2->Draw();
  cHist->cd(++Ihist);		hCosbw1Cosbw2->Draw();
  cHist->cd(++Ihist);		hMw1Mw2->Draw();
  cHist->cd(++Ihist);		hEvisPl->Draw();
  cHist->cd(++Ihist);		hMt1Mt2->Draw();
  cHist->cd(++Ihist);		hThrust->Draw();
  cHist->cd(++Ihist);		hYcut->Draw();

  cHist->Update();

  last->cd();
}

//_________________________________________________________
Bool_t TTH8JAnalysis::Process(Int_t ev)
{
  // Local copies of TTH8JAnalysisBuf data members.

  Int_t     	fNtracks;	// track multiplicity
  Double_t  	fEvis;		// visible energy
  Double_t  	fPt;		// Pt
  Double_t  	fPl;		// Pl
  Double_t  	fYcut;		// y_cut to force the event to 4 jets
  Int_t        	fNjets;		// jet multiplicity
  Double_t  	fThrust;	// thrust

  // Remember the previous directory.

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

  JSFSIMDST     *sds     = (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
  JSFSIMDSTBuf  *evt     = (JSFSIMDSTBuf*)sds->EventBuf();
  TTH8JAnalysisBuf *ua    = (TTH8JAnalysisBuf *)fEventBuf;
  TTH8JAnalysisBuf &a     = *ua;

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

  fEvis = qsum(0);		// E_vis
  fPt   = qsum.GetPt();		// P_t
  fPl   = qsum(3);		// P_l

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

  if (gDEBUG) cerr << "Ycut = " << fYcut << " Njets = " << fNjets << endl;

  // Make sure that No. of jets is xNjets.

  hYcut->Fill(fYcut);
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

  // Find W and top candidates in given mass windows.

  TObjArray solutions(10);
  ANLPairCombiner w1candidates(jets,jets);
  if (gDEBUG) {
    cerr << "------------------------------------------" << endl;
    cerr << "- w1candidates:" << endl;
    w1candidates.DebugPrint();
  }
  ANLPair *w1p, *w2p, *bbp, *hp;
  
  ANLVTXTagger btag(2.,2);
  
  while ((w1p = (ANLPair *)w1candidates())) {
    ANLPair &w1 = *w1p;
    Double_t w1mass = w1().GetMass();
    if (TMath::Abs(w1mass - kMassW) > xM2j) continue;
    w1.LockChildren();
    ANLPairCombiner w2candidates(w1candidates);
    if (gDEBUG) {
      cerr << "-- w2candidates:" << endl;
      w2candidates.DebugPrint();
    }
    while ((w2p = (ANLPair *)w2candidates())) {
      ANLPair &w2 = *w2p;
      if (w2.IsLocked()) continue;
      Double_t w2mass = w2().GetMass();
      if (TMath::Abs(w2mass - kMassW) > xM2j) continue;
      w2.LockChildren();
      ANLPairCombiner bbcandidates(w2candidates);
      bbcandidates.Reset();
      if (gDEBUG) {
        cerr << "---- bbcandidates:" << endl;
        bbcandidates.DebugPrint();
      }
      while ((bbp = (ANLPair *)bbcandidates())) {
        ANLPair &bb = *bbp;
        if (bb.IsLocked()) continue;
        if (!btag(*(ANLJet *)bb[0]) || !btag(*(ANLJet *)bb[1])) continue;
        bb.LockChildren();
        for (Int_t i = 0; i < 2; i++) {
          ANL4DVector *b1p = (ANL4DVector *)bb[i];
          ANL4DVector *b2p = (ANL4DVector *)bb[1-i];
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
          }
          if (gDEBUG) {
            cerr << " M_w1 = " << w1mass << " M_w2 = " << w2mass << endl
                 << " M_t1 = " << t1mass << " M_t2 = " << t2mass << endl;
            cerr << " w1p  = " << (void *)w1p
                 << " w2p  = " << (void *)w2p
                 << " bbp  = " << (void *)bbp << endl;
            cerr << " w1[0] = " << (void *)w1[0]
                 << " w1[1] = " << (void *)w1[1]
                 << " w2[0] = " << (void *)w2[0]
                 << " w2[1] = " << (void *)w2[1]
                 << " bb[0] = " << (void *)bb[0]
                 << " bb[1] = " << (void *)bb[1] << endl;
          }
          ANLPairCombiner hcandidates(bbcandidates);
	  hcandidates.Reset();
	  Bool_t ok = kFALSE;
	  while ((hp = (ANLPair *)hcandidates())) { 
	    ANLPair &h = *hp;
	    Double_t hmass = h().GetMass();
	    if (h.IsLocked() ||
	       !btag(*(ANLJet *)h[0]) || !btag(*(ANLJet *)h[1]) ||
	        TMath::Abs(hmass - kMassH) > xM2j) continue;
            Double_t chi2 = TMath::Power((w1mass - kMassW)/kSigmaMw,2.)
                          + TMath::Power((w2mass - kMassW)/kSigmaMw,2.)
                          + TMath::Power((t1mass - kMasst)/kSigmaMt,2.)
                          + TMath::Power((t2mass - kMasst)/kSigmaMt,2.)
                          + TMath::Power((hmass  - kMassH)/kSigmaMh,2.);
	    ANLPair *ttp = new ANLPair(bw1p,bw2p);
	    ANLPair *hpp = new ANLPair(h);
            solutions.Add(new ANLPair(ttp,hpp,chi2));
            ok = kTRUE;
            // hEw1Ew2->Fill(w1()(0),w2()(0),1.0);
            // hMw1Mw2->Fill(w1mass,w2mass,1.0);
            // hMt1Mt2->Fill(t1mass,t2mass,1.0);
          }
          if (!ok) {
	    delete bw1p;
	    delete bw2p;
	  }
        }
        bb.UnlockChildren();
      }
      w2.UnlockChildren();
    }
    w1.UnlockChildren();
  }

  // Cut on No. of solutions.

  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|m_jj - m_W| <= %g && |m_3j - m_t| <= %g",xM2j,xM3j);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on cos(theta_bW).

  TIter nextsol(&solutions);
  ANLPair *sol;
  while ((sol = (ANLPair *)nextsol())) {
    ANLPair &tt  = *(ANLPair *)(*sol)[0];
    ANLPair &bw1 = *(ANLPair *)tt[0];
    ANLPair &bw2 = *(ANLPair *)tt[1];
    ANL4DVector &b1 = *(ANL4DVector *)bw1[0];
    ANL4DVector &w1 = *(ANL4DVector *)bw1[1];
    ANL4DVector &b2 = *(ANL4DVector *)bw2[0];
    ANL4DVector &w2 = *(ANL4DVector *)bw2[1];
    Double_t cosbw1 = b1.CosTheta(w1);
    Double_t cosbw2 = b2.CosTheta(w2);
    hCosbw1Cosbw2->Fill(cosbw1,cosbw2,1.0);
    if (cosbw1 > xCosbw || cosbw2 > xCosbw) {
      solutions.Remove(sol);
      sol->Delete();
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_bw)| > %g",xCosbw);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Thrust.

  ANLEventShape eshape;
  eshape.Initialize(tracks);
  fThrust = eshape.GetThrust();
  hThrust->Fill(fThrust);
  if ( fThrust > xThrust ) {
  	nextsol.Reset();
  	while ((sol = (ANLPair *)nextsol())) sol->Delete();
  	CleanUp(&solutions);
  	CleanUp(&tracks);
  	return kFALSE;
  }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Thrust <= %g",xThrust);
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

  // Now store this in TTH8JAnalysisBuf.

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
    ANLPair &tt  = *(ANLPair *)(*sol)[0];
    ANLPair &h   = *(ANLPair *)(*sol)[1];
    ANLPair &bw1 = *(ANLPair *)tt[0];
    ANLPair &bw2 = *(ANLPair *)tt[1];
    ANL4DVector &w1 = *(ANL4DVector *)bw1[1];
    ANL4DVector &w2 = *(ANL4DVector *)bw2[1];
    Double_t chi2   = sol->GetQuality();
    Double_t w1mass = w1.GetMass();
    Double_t w2mass = w2.GetMass();
    Double_t t1mass = bw1().GetMass();
    Double_t t2mass = bw2().GetMass();
    Double_t hmass  = h.GetMass();
          hChi2->Fill(chi2);
          hEw1Ew2->Fill(w1(0),w2(0),1.0);
          hMw1Mw2->Fill(w1mass,w2mass,1.0);
          hMt1Mt2->Fill(t1mass,t2mass,1.0);
          hMw2Mh->Fill(w2mass,hmass,1.0);
  }

  // Clean up

  nextsol.Reset();
  while ((sol = (ANLPair *)nextsol())) sol->Delete();
  CleanUp(&solutions);
  CleanUp(&tracks);

  last->cd();
  return kTRUE;
}

//_________________________________________________________
Bool_t TTH8JAnalysis::Terminate()
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


