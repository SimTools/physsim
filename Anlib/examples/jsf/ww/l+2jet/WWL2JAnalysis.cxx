//*************************************************************************
//* =======================
//*  WWL2JAnalysis Classes
//* =======================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC W boson pair data. 
//* (Requires)
//*     library Anlib
//*     library WWStudy
//* (Provides)
//*     class WWL2JAnalysis
//*     class WWL2JAnalysisBuf
//* (Usage)
//*   Take a look at anlL2J.C.  
//* (Update Recored)
//*   2000/04/28  K.Ikematsu    Derived from WW4JAnalysis.h.
//*
//*************************************************************************
//
#include "WWL2JAnalysis.h"

static const Double_t kMassW   = 80.00; // W mass
static const Double_t kMassZ   = 91.19; // Z mass
static const Double_t kSigmaMw =   4.0; // W mass resolution
static const Double_t kSigmaMz =   4.0; // W mass resolution
static const Int_t    kZoneX   = 4;     // No. X Zones in the Canvas
static const Int_t    kZoneY   = 4;     // No. Y Zones in the Canvas

Int_t WWL2JAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

//_____________________________________________________________________
//  ----------------------
//  WWL2JAnalysisBuf Class
//  ----------------------
//
//
ClassImp(WWL2JAnalysisBuf)

//_________________________________________________________
WWL2JAnalysisBuf::WWL2JAnalysisBuf(const Char_t *name, const Char_t *title, 
   WWL2JAnalysis *module) : JSFEventBuf(name, title, (JSFModule*)module) {}

//_________________________________________________________
WWL2JAnalysisBuf::WWL2JAnalysisBuf(WWL2JAnalysis *module, const Char_t *name,
   const Char_t *title) : JSFEventBuf(name, title, (JSFModule*)module) {}

//_____________________________________________________________________
//  -------------------
//  WWL2JAnalysis Class
//  -------------------
//
//

ClassImp(WWL2JAnalysis)

WWL2JAnalysis::WWL2JAnalysis(const Char_t *name, const Char_t *title)
	       : JSFModule(name, title) 
{
  fEventBuf = new WWL2JAnalysisBuf(this);
  SetBufferSize(2000);  // buffer size for event data.
  cout << "WWL2JAnalysisBuf is created...fEventBuf is " 
       << (Int_t)fEventBuf << endl;
}

//_____________________________________________________________________
WWL2JAnalysis::~WWL2JAnalysis()
{
  cout << "WWL2JAnalysisBuf will be deleted...fEventBuf is " 
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}

//_____________________________________________________________________
void WWL2JAnalysis::CleanUp(TObjArray *objs)
{
  TIter next(objs);
  TObject *obj;
  while ((obj = next())) {
  	objs->Remove(obj);
  	delete obj;
  }
}

//_____________________________________________________________________
Bool_t WWL2JAnalysis::Initialize()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat       = new TH1F("hStat","Cut Statistics",  20,   0.0,  20.0);
  hNtracks    = new TH1F("hNtracks","No. tracks" ,  50,   0.0, 100.0);
  hEvis       = new TH1F("hEvis","Visible energy",  50,   0.0, 500.0);
  hPt         = new TH1F("hPt","Missing Pt"      ,  50,   0.0, 250.0);
  hNlptracks  = new TH1F("hNlptracks","No.lptrks",  10,   0.0,  10.0);
  hYcut       = new TH1F("hYcut","Ymax"           ,100,   0.0,   1.0);
  hNjets      = new TH1F("hNjets","No. jets"     ,  20,   0.0,  20.0);
  hEjet       = new TH1F("hEjet","Jet energy"    ,  50,   0.0, 200.0);
  hCosjet     = new TH1F("hCosjet","cos(theta_j)",  50,  -1.0,  +1.0);
  hEw1Ew2     = new TH2F("hEw1Ew2","(E_w1,E_w2)" ,  
  				  50,  0.0, 400.0,  50,   0.0, 400.0);
  hCosw1Cosw2 = new TH2F("hCosw1Cosw2","(cos_w1,cow_w2)",
				  50, -1.0,  +1.0,  50,  -1.0,  +1.0);
  hMw1Mw2     = new TH2F("hMw1Mw2","(m_w1,m_w2)" , 
  				  70, 50.0, 130.0,  60,  50.0, 110.0);
  hEvisPl     = new TH2F("hEvisPl","(Evis,Pl)"   , 
  				  60,  0.0, 600.0,  50,-100.0,+100.0);
  hAcop       = new TH1F("hAcop","Acoplanarity"  ,  90,   0.0, 180.0);

  xNtracks  =      0;   // No. of tracks
  xEtrack   =    0.0;   // track energy
  xEvis     =    0.0;   // Maximum visible energy
  xPt       =    0.0;   // Pt maximum
  xPl       =    0.0;   // Pl maximum
  xElepton  =    0.0;   // Elepton mimimum
  xCosCone  =    0.0;   // cos(theta_cone)
  xEcone    =    0.0;   // Econe maximum
  xYcut     =    0.0;   // y_cut to force the event to 2 jets
  xNjets    =      0;   // No. of jets
  xEjet     =    0.0;   // E_jet minimum
  xCosjet   =    0.0;   // |cos(theta_j)| maximum
  xCosw     =    0.0;	// |cos(theta_w)| maximum
  xM2j      =    0.0;	// |m_jj-m_W| maximum
  xAcop     =    0.0;	// Acoplanarity maximum

  for (Int_t i = 0; i < MAXCUT; i++) {
    strcpy(&cutName[i][0],"     ");
  }

  last->cd();
  return 0;
}

//_________________________________________________________
void WWL2JAnalysis::DrawHist()
{
  TDirectory *last = gDirectory;
  if (!cHist) {
    cHist = new TCanvas("cHist","Canvas 1",10, 10, kZoneX*200, kZoneY*200);
    cHist->Divide(kZoneX,kZoneY);
  } else {
    cHist->cd();
  }

  Int_t Ihist = 0;
  cHist->cd(++Ihist);	hStat->Draw();
  cHist->cd(++Ihist);	hNtracks->Draw();
  cHist->cd(++Ihist);	hEvis->Draw();
  cHist->cd(++Ihist);	hPt->Draw();
  cHist->cd(++Ihist);	hNlptracks->Draw();
  cHist->cd(++Ihist);	hNjets->Draw();
  cHist->cd(++Ihist);	hYcut->Draw();
  cHist->cd(++Ihist);	hEjet->Draw();
  cHist->cd(++Ihist);	hCosjet->Draw();
  cHist->cd(++Ihist);	hEw1Ew2->Draw();
  cHist->cd(++Ihist);	hCosw1Cosw2->Draw();
  cHist->cd(++Ihist);	hMw1Mw2->Draw();
  cHist->cd(++Ihist);	hEvisPl->Draw();
  cHist->cd(++Ihist);	hAcop->Draw();

  cHist->Update();

  last->cd();
}

//_________________________________________________________
Bool_t WWL2JAnalysis::Process(Int_t ev)
{
  // Local copies of WWL2JAnalysisBuf data members.

  Int_t     	fNtracks;	// track multiplicity
  Int_t     	fNlptracks;     // Isolated lepton multiplicity
  Double_t  	fEvis;		// visible energy
  Double_t  	fPt;		// Pt
  Double_t  	fPl;		// Pl
  Double_t  	fYcut;		// y_cut to force the event to 2 jets
  Int_t        	fNjets;		// jet multiplicity
  Double_t      fAcop;          // Acoplanarity

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
  WWL2JAnalysisBuf *ua    = (WWL2JAnalysisBuf *)fEventBuf;
  WWL2JAnalysisBuf &a     = *ua;

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
  if ( fEvis > xEvis ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"E_vis < %g",xEvis);
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

  // Find Isolated Lepton.

  TObjArray lptracks(20);
  fNlptracks = 0;
  TIter nexttrk (&tracks);
  ANLTrack *trkp;
  Double_t lpcharge = 0.;
  while ((trkp = (ANLTrack *)nexttrk())) {
    ANLTrack &trk = *trkp;
    if ( ! trk.IsLepton() ) continue;
    Double_t elepton = trk.E();
    if ( elepton < xElepton ) continue;
    Double_t econe = trk.GetConeEnergy(xCosCone, &tracks);
    if ( econe <= xEcone ) {
      lptracks.Add(trkp);
      fNlptracks++;
      trk.Lock();
      lpcharge = trk.GetCharge();
    }
  }

  // Require only one Isolated Lepton.

  hNlptracks->Fill(fNlptracks);
  if ( fNlptracks != 1 ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Nlptracks = 1");
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

  // Find leptonic decayed W and hadronic decayed W.

  ANL4DVector qcm(a.GetEcm());
  ANL4DVector qnu = qcm - qsum;

  ANLJet lpjet(lptracks);
  ANLJet nujet;
  nujet.Add(&qnu);
  ANLPair w1(&lpjet, &nujet);
  ANLPair w2(jets[0],jets[1]);

  // Cut on cos(theta_W).

  Double_t ew1 = w1.E();
  Double_t ew2 = w2.E();
  hEw1Ew2->Fill(ew1,ew2,1.0);
  Double_t cosw1 = w1.CosTheta();
  Double_t cosw2 = w2.CosTheta();
  hCosw1Cosw2->Fill(cosw1,cosw2,1.0);  
  if (TMath::Abs(cosw1) > xCosw || TMath::Abs(cosw2) > xCosw) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_w)| <= %g",xCosw);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Acop.

  fAcop = w1.Acop(w2);
  hAcop->Fill(fAcop);
  if (fAcop > xAcop) { CleanUp(&tracks); return kFALSE; }
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
       << "Event " << gJSF->GetEventNumber() << endl
       << "------------------------------------------" << endl;

  // Now store this in WWL2JAnalysisBuf.

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

  hEvisPl->Fill(fEvis,fPl,1.);
  Double_t w1mass = w1.GetMass();
  Double_t w2mass = w2.GetMass();
  hMw1Mw2->Fill(w1mass,w2mass,1.0);

  // Clean up

  CleanUp(&tracks);

  last->cd();
  return kTRUE;
}

//_________________________________________________________
Bool_t WWL2JAnalysis::Terminate()
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
