//*************************************************************************
//* ========================
//*  TT2L2JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC ttbar data to select 2-lepton + 2-jet events.
//* (Requires)
//* 	library Anlib
//* 	library TTStudy
//* (Provides)
//* 	class TT2L2JAnalysis
//* 	class TT2L2JAnalysisBuf
//* (Usage)
//*   Take a look at anl2L2J.C.
//* (Update Recored)
//*   1999/08/19  K.Ikematsu    Derived from TTL4JAnalysis.cxx.
//*   2001/07/07  K.Ikematsu    Modified for MacOS X.
//*
//* $Id$
//*************************************************************************
//
#include "TT2L2JAnalysis.h"

static const Double_t kMassW   = 80.00; 	// W mass
static const Double_t kMassZ   = 91.19; 	// Z mass
static const Double_t kMasst   = 170.0; 	// top mass
static const Double_t kSigmaMw =   4.0; 	// W mass resolution
static const Double_t kSigmaMz =   4.0; 	// Z mass resolution
static const Double_t kSigmaMt =  15.0; 	// top mass resolution
static const Int_t    kZoneX   =     4;		// No. X Zones in the Canvas
static const Int_t    kZoneY   =     4;		// No. Y Zones in the Canvas

Int_t TT2L2JAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

//_____________________________________________________________________
//  -----------------------
//  TT2L2JAnalysisBuf Class
//  -----------------------
//
//
ClassImp(TT2L2JAnalysisBuf)

//_________________________________________________________
TT2L2JAnalysisBuf::TT2L2JAnalysisBuf(const Char_t *name, const Char_t *title,
   TT2L2JAnalysis *module) : JSFEventBuf(name, title, (JSFModule*)module) {}

//_________________________________________________________
TT2L2JAnalysisBuf::TT2L2JAnalysisBuf(TT2L2JAnalysis *module, const Char_t *name,
   const Char_t *title) : JSFEventBuf(name, title, (JSFModule*)module) {}


//_____________________________________________________________________
//  --------------------
//  TT2L2JAnalysis Class
//  --------------------
//
//

ClassImp(TT2L2JAnalysis)

TT2L2JAnalysis::TT2L2JAnalysis(const Char_t *name, const Char_t *title)
  : JSFModule(name, title), cHist(0)
{
  fEventBuf = new TT2L2JAnalysisBuf(this);
  SetBufferSize(2000);  // buffer size for event data.
  cout << "TT2L2JAnalysisBuf is created...fEventBuf is "
       << (Int_t)fEventBuf << endl;
}

//_____________________________________________________________________
TT2L2JAnalysis::~TT2L2JAnalysis()
{
  cout << "TT2L2JAnalysisBuf will be deleted...fEventBuf is "
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}

//_____________________________________________________________________
void TT2L2JAnalysis::CleanUp(TObjArray *objs)
{
  objs->SetOwner();
}

//_____________________________________________________________________
Bool_t TT2L2JAnalysis::Initialize()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat         = new TH1F("hStat","Cut Statistics" ,  20,   0.0,  20.0);
  hNtracks      = new TH1F("hNtracks","No. tracks"  ,  50,   0.0, 200.0);
  hEvis         = new TH1F("hEvis","Visible energy" ,  80,   0.0, 400.0);
  hPt           = new TH1F("hPt","Missing Pt"       ,  60,   0.0, 120.0);
  hNlptracks    = new TH1F("hNlptracks","No.lp trks",  10,   0.0,  10.0);
  hNlpmtracks   = new TH1F("hNlpmtracks","# l- trks",   5,   0.0,   5.0);
  hNlpptracks   = new TH1F("hNlpptracks","# l+ trks",   5,   0.0,   5.0);
  hElpmlpp      = new TH2F("hElpmlpp","(E_l-,E_l+)" ,
   				     50,  0.0, 100.0,  50,   0.0, 100.0);
  hCoslpmlpp    = new TH2F("hCoslpmlpp","(Cos_l-,Cos_l+)",
                                     50, -1.0,  +1.0,  50,  -1.0,  +1.0);
  hNjets        = new TH1F("hNjets","No. jets"      ,  20,   0.0,  20.0);
  hEjet         = new TH1F("hEjet","Jet energy"     ,  50,   0.0, 100.0);
  hCosjet       = new TH1F("hCosjet","cos(theta_j)" ,  50,  -1.0,  +1.0);
  hEvisPl       = new TH2F("hEvisPl","(Evis,Pl)"    ,
  				     60,  0.0, 600.0,  50,-100.0,+100.0);
  hThrust       = new TH1F("hThrust","Thrust"       ,  50,   0.0,   1.0);
  hYcut         = new TH1F("hYcut","Ymax"           , 100,   0.0,   1.0);

  xNtracks  =     10;   // No. of tracks
  xEtrack   =   0.10;   // track energy
  xEvis     = 300.00;   // Maximum visible energy
  xPt       = 999.00;   // Pt maximum
  xPl       = 999.00;   // Pl maximum
  xElepton  =  20.00;   // Elepton mimimum
  xCosCone  =   0.94;   // cos(theta_cone)
  xEcone    =  10.00;   // Ecome maximum
  xYcut     =  0.004;   // y_cut to force the event to 2 jets
  xNjets    =      2;   // No. of jets
  xEjet     =   5.00;	// E_jet minimum
  xCosjet   =   0.99;	// |cos(theta_j)| maximum
  xThrust   =   1.00;   // Thrust maximum

  for (Int_t i = 0; i < MAXCUT; i++) {
    strcpy(&cutName[i][0],"     ");
  }

  last->cd();
  return 0;
}

//_________________________________________________________
void TT2L2JAnalysis::DrawHist()
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
  cHist->cd(++Ihist);   hNlptracks->Draw();
  cHist->cd(++Ihist);   hNlpmtracks->Draw();
  cHist->cd(++Ihist);   hNlpptracks->Draw();
  cHist->cd(++Ihist);	hElpmlpp->Draw();
  cHist->cd(++Ihist);	hCoslpmlpp->Draw();
  cHist->cd(++Ihist);	hNjets->Draw();
  cHist->cd(++Ihist);	hEjet->Draw();
  cHist->cd(++Ihist);	hCosjet->Draw();
  cHist->cd(++Ihist);	hEvisPl->Draw();
  cHist->cd(++Ihist);	hThrust->Draw();
  cHist->cd(++Ihist);	hYcut->Draw();

  cHist->Update();

  last->cd();
}

//_________________________________________________________
Bool_t TT2L2JAnalysis::Process(Int_t ev)
{
  // Local copies of TT2L2JAnalysisBuf data members.

  Int_t     	fNtracks;	// track multiplicity
  Double_t  	fEvis;		// visible energy
  Double_t  	fPt;		// Pt
  Double_t  	fPl;		// Pl
  Int_t         fNlptracks;     // Isolated lepton multiplicity
  Int_t         fNlpmtracks;    // No. of positive charged leptons
  Int_t         fNlpptracks;    // No. of naegative charged leptons
  Double_t  	fYcut;		// y_cut to force the event to 2 jets
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
  TT2L2JAnalysisBuf *ua    = (TT2L2JAnalysisBuf *)fEventBuf;
  TT2L2JAnalysisBuf &a     = *ua;

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
  fNlpptracks = 0;
  fNlpmtracks = 0;
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
      if ( trk.GetCharge() < 0. ) {
	//	cerr << "charge is " << trk.GetCharge() << endl;
	fNlpmtracks++;
      } else {
	//	cerr << "charge is " << trk.GetCharge() << endl;
	fNlpptracks++;
      }
      trk.Lock();
      lpcharge = trk.GetCharge();
    }
  }

  // Require two Isolated Leptons.

  hNlptracks->Fill(fNlptracks);
  if ( fNlptracks != 2 ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Nlptracks = 2");
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Require Charge oppositeness of two Isolated Leptons.

  hNlpmtracks->Fill(fNlpmtracks);
  hNlpptracks->Fill(fNlpptracks);
  if ( fNlpmtracks != 1 && fNlpptracks != 1) {
    CleanUp(&tracks); return kFALSE;
  }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Charge oppositeness of two Leptons");
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

  // Cut on Thrust.

  ANLEventShape eshape;
  eshape.Initialize(tracks);
  fThrust = eshape.GetThrust();
  hThrust->Fill(fThrust);
  if ( fThrust > xThrust ) { CleanUp(&tracks); return kFALSE; }
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
       << "Event " << gJSF->GetEventNumber() << endl
       << "------------------------------------------" << endl;

  // Now store this in TT2L2JAnalysisBuf.

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

  ANLTrack qm, qp;

  if ( ((ANLTrack *)lptracks.UncheckedAt(0))->GetCharge() < 0. ) {
    qm = *((ANLTrack *)lptracks.UncheckedAt(0));
    qp = *((ANLTrack *)lptracks.UncheckedAt(1));
  } else {
    qp = *((ANLTrack *)lptracks.UncheckedAt(0));
    qm = *((ANLTrack *)lptracks.UncheckedAt(1));
  }
  Double_t elpm = qm.E();
  Double_t coslpm = qm.CosTheta();
  Double_t elpp = qp.E();
  Double_t coslpp = qp.CosTheta();
  hElpmlpp->Fill(elpm,elpp,1.);
  hCoslpmlpp->Fill(coslpm,coslpp,1.);

  //  Double_t eelm = qm.IsElectron();
  //  Double_t eelp = qp.IsElectron();
  //  Double_t emum = qm.IsMuon();
  //  Double_t emup = qp.IsMuon();

  hEvisPl->Fill(fEvis,fPl,1.);

  CleanUp(&tracks);

  last->cd();
  return kTRUE;
}

//_________________________________________________________
Bool_t TT2L2JAnalysis::Terminate()
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
