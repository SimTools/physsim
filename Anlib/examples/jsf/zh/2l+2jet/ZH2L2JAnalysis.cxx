//***************************************************************************
//*  ====================
//*  ZH2L2JAnalysis Classes
//*  ====================
//*
//*  (Description)
//*	A user analysis class for JLC analyses.
//*	This reads and analyzes MC e+e- -> ZH data.
//*  (Requires)
//*	library Anlib (in physsim-99a-1 if K.Fujii)
//*	library ZHStudy (also in physsim)
//*  (Provides)
//*	class ZH2L2JAnalysis
//*	class ZH2L2JAnalysisBuf
//*  (Usage)
//*	...
//*  (Update Record)
//*	18 Nov 1999	A.L.C.Sanchez	Patterned after examples (WW4JAnalysis)
//*					with physsim-99a-1 (still
//*					largely for correction!)
//*	22 Nov 1999	A.L.C.Sanchez	Modified to suit ZH 4jet analysis
//***************************************************************************

#include "ZH2L2JAnalysis.h"

static const Double_t kMassH   = 120.0; // H mass
static const Double_t kMassZ   = 91.19;	// Z mass
static const Double_t kSigmaMh =   4.0;	// H mass resolution
static const Double_t kSigmaMz =   4.0;	// Z mass resolution
static const Int_t    kZoneX   =     4;	// No. of X Zones in the Canvas
static const Int_t    kZoneY   =     4;	// No. of Y Zones in the Canvas

Int_t ZH2L2JAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ---------------------
//  ZH2L2JAnalysisBuf Class
//  ---------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ClassImp(ZH2L2JAnalysisBuf)

// Constructors
ZH2L2JAnalysisBuf::ZH2L2JAnalysisBuf(const Char_t *name, const Char_t *title,
   ZH2L2JAnalysis *module) : JSFEventBuf(name, title, (JSFModule*)module) {}

ZH2L2JAnalysisBuf::ZH2L2JAnalysisBuf(ZH2L2JAnalysis *module, const Char_t *name,
   const Char_t *title) : JSFEventBuf(name, title, (JSFModule*)module) {}

ZH2L2JAnalysisBuf::~ZH2L2JAnalysisBuf() {}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ------------------
//  ZH2L2JAnalysis Class
//  ------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ClassImp(ZH2L2JAnalysis)

// Constructor
ZH2L2JAnalysis::ZH2L2JAnalysis(const Char_t *name, const Char_t *title)
               : JSFModule(name, title)
{
  fEventBuf = new ZH2L2JAnalysisBuf(this);
  SetBufferSize(2000);	// buffer size for event data
  cout << "ZH2L2JAnalysis is created... fEventBuf is "
       << (Int_t)fEventBuf << endl;
}

// Destructor
ZH2L2JAnalysis::~ZH2L2JAnalysis()
{
  cout << "ZH2L2JAnalysisBuf will be deleted... fEventBuf is "
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}


//####### Functions ###########

// *****CleanUp()***** //
void ZH2L2JAnalysis::CleanUp(TObjArray *objs)
{
  TIter next(objs);
  TObject *obj;
  while ( (obj = next()) ) {
    objs->Remove(obj);
    delete obj;
  }
}
// *****End of CleanUp()***** //


// ******Initialize()****** //
Bool_t ZH2L2JAnalysis::Initialize()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat       = new TH1F("hStat","Cut Statistics"   ,  20,   0.0,   20.0);
  hNtracks    = new TH1F("hNtracks","No. of Tracks" ,  50,   0.0,  100.0);
  hEvis       = new TH1F("hEvis","Visible Energy"   ,  50,   0.0,  500.0);
  hPt         = new TH1F("hPt","Missing Pt"         ,  50,   0.0,  250.0);
  hNlptracks    = new TH1F("hNlptracks","No.lp trks", 10,   0.0,  10.0);
  hNjets      = new TH1F("hNjets","No. of Jets"     ,  20,   0.0,   20.0);
  hEjet       = new TH1F("hEjet","Jet Energy"       ,  50,   0.0,  200.0);
  hCosjet     = new TH1F("hCosjet","cos(theta_j)"   ,  50,  -1.0,   +1.0);
  hNsols      = new TH1F("hNsols","No. of solutions",  20,   0.0,   20.0);
  hChi2       = new TH1F("hChi2","Chi2"             ,  50,   0.0,   50.0);
  hEvisPl     = new TH2F("hEvisPl","(Evis,Pl)"      ,
                                     60,  0.0, 600.0,  50,-100.0, +100.0);
  hEzEh       = new TH2F("hEzEh","E_Z vs. E_H"      ,
				     50,  0.0, 500.0,  50,   0.0,  500.0);
  hCoszCosh   = new TH2F("hCoszCosh","Cos(theta_Z) vs. Cos(theta_H)",
                                     50, -1.0,  +1.0,  50,  -1.0,   +1.0);
  hAcop       = new TH1F("hAcop","Acoplanarity"     ,  90,   0.0,  180.0);
  hMassh      = new TH1F("hMassh","M_H Distribution", 600,  90.0,  150.0);
  hMassz      = new TH1F("hMassz","M_Z Distribution", 240,  60.0,  120.0);

  xNtracks  =     25;	// No. of Tracks
  xEtrack   =   0.10;	// track energy
  xEvis     = 100.00;	// Minimum visible energy
  xPt       =  10.00;	// Pt minimum
  xPl       = 999.00;	// Pl maximum
  xElepton  =  10.00;   // Elepton mimimum
  xCosCone  =   0.94;   // cos(theta_cone)
  xEcone    =  10.00;   // Ecome maximum
  xYcut     =  0.004;	// y_cut to force the event to 4 jets
  xNjets    =      2;	// No. of Jets
  xEjet     =   5.00;	// E_jet minimum
  xCosjet   =   0.99;	// |cos(theta_j)| maximum
  xCoszh    =   0.99;	// |cos(theta_Z)| and |cos(theta_H)| maximum
  xM2j      =  18.00;	// |m_jj-m_S| maximum , S = Z,H
  xAcop     =  30.00;	// Acoplanarity maximum

  last->cd();
  return 0;
}
// *****End of Initialize()***** //


// *****DrawHist()***** //
void ZH2L2JAnalysis::DrawHist()
{
  TDirectory *last = gDirectory;
  if (!cHist) {
    cHist = new TCanvas("cHist","Canvas 1", 10, 10, kZoneX*200, kZoneY*200);
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
  cHist->cd(13);	hEvis->Draw();
  cHist->cd(14);	hAcop->Draw();

  cHist->Update();

  last->cd();
}
// *****End of DrawHist()***** //


// *****Process()***** //
Bool_t ZH2L2JAnalysis::Process(Int_t ev)
{
  // Local copies of ZH2L2JAnalysisBuf data members.

  Int_t		fNtracks;	// track multiplicity
  Int_t         fNlptracks;     // Isolated lepton multiplicity
  Double_t	fEvis;		// visible energy
  Double_t	fPt;		// Pt
  Double_t	fPl;		// Pl
  Double_t	fYcut;		// y_cut to force the event to 4 jets
  Int_t		fNjets;		// jet multiplicity

  // Remember the previous directory.
  
  TDirectory *last = gDirectory;
  gFile->cd("/");

  Char_t msg[60];

  // Analysis starts here.

  Float_t selid = -0.5;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) strcpy(&cutName[(Int_t)selid][0],"No Cut");

  // Get event buffer and make combined tracks accessible.

  JSFSIMDST     *sds	= (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
  JSFSIMDSTBuf	*evt	= (JSFSIMDSTBuf*)sds->EventBuf();
  ZH2L2JAnalysisBuf *ua	= (ZH2L2JAnalysisBuf*)fEventBuf;
  ZH2L2JAnalysisBuf &a	= *ua;

  Int_t		ntrks	= evt->GetNLTKCLTracks();	// No. of tracks
  TObjArray	*trks	= evt->GetLTKCLTracks();	// Combined tracks

  // Select good tracks
  ANL4DVector qsum;
  TObjArray tracks(1000);
  fNtracks = 0;
  for ( Int_t i = 0; i < ntrks; i++) {
    JSFLTKCLTrack *t = (JSFLTKCLTrack*)trks->UncheckedAt(i);
    if ( t->GetE() > xEtrack ) {
      ANLTrack *qt = new ANLTrack(t);
      tracks.Add(qt);		// track 4-momentum
      qsum += *qt;		// total 4-mometum
      fNtracks++;		// *qt stays.
    }
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
  
  fEvis = qsum(0);	// E_vis
  fPt	= qsum.GetPt(); // P_t
  fPl	= qsum(3);	// P_l

  if (gDEBUG) cerr << "Evis = " << fEvis << " Pt = "
	<< fPt << " Pl = " << fPl << endl;

  // Cut on Evis.

  hEvis->Fill(fEvis);
  if ( fEvis < xEvis ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Evis > %g",xEvis);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Pt.
  
  hPt->Fill(fPt);
  if ( fPt < xPt ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg," Pt <= %g",xPt);
    strcpy(&cutName[(Int_t)selid][0], msg);
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
      lpcharge += trk.GetCharge();
    }
  }

  // Require only two Isolated Leptons.

  hNlptracks->Fill(fNlptracks);
  if ( fNlptracks != 2 || lpcharge ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Nlptracks = 2");
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
  if ( fNjets < xNjets ) { CleanUp(&tracks); return kFALSE;}
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Njets >= %i for Ycut = %g",xNjets,xYcut);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Now force the event to be xNjets.

  jclust.ForceNJets(xNjets);
  fNjets = jclust.GetNjets();
  fYcut  = jclust.GetYcut();

  if(gDEBUG) cerr << "Ycut = " << fYcut << " Njets = " << fNjets << endl;

  // Make sure the No. of Jets is xNjets.
  
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
  while ( (jetp = (ANLJet*)nextjet()) ) {
    ANLJet &jet = *jetp;
    if (gDEBUG) jet.DebugPrint();
    Double_t ejet = jet()(0);
    if ( ejet < ejetmin ) ejetmin = ejet;
    hEjet->Fill(ejet);
    Double_t cosj = jet.CosTheta();
    if ( TMath::Abs(cosj) > TMath::Abs(cosjmax) ) cosjmax = cosj;
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
    sprintf(msg,"|cos(theta_j)| <= %g", xCosjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Unlock lepton tracks.

  TIter nextlptrk(&lptracks);
  while ((trkp = (ANLTrack *)nextlptrk())) trkp->Unlock();

  // Find ZH candidates in given mass window.
  
  TObjArray solutions(10);
  ANLPairCombiner zcandidates(jets,jets);
  ANLPair *hp, *zp;
  while ( (zp = (ANLPair*)zcandidates()) ) {
    ANLPair &z = *zp;
    Double_t zmass = z().GetMass();
    if (TMath::Abs(zmass - kMassZ) > xM2j) continue;
    ANLPairCombiner hcandidates(lptracks,lptracks);
    while ( (hp = (ANLPair*)hcandidates()) ) {
      ANLPair &h = *hp;
      if (h.IsLocked()) continue;
      Double_t hmass = h().GetMass();
      if (TMath::Abs(hmass - kMassH) > xM2j) continue;
      if (gDEBUG) {
        cerr << " M_Z = " << zmass << " M_H = " << hmass << endl;
        cerr << " zp  = " << (void*)zp
             << " hp  = " << (void*)hp << endl;
        cerr << " z[0] = " << (void*)z[0]
             << " z[1] = " << (void*)z[1]
             << " h[0] = " << (void*)h[0]
             << " h[1] = " << (void*)h[1] << endl;
      }
      Double_t chi2 = TMath::Power((zmass - kMassZ)/kSigmaMz,2.)
                    + TMath::Power((hmass - kMassH)/kSigmaMh,2.);
      solutions.Add(new ANLPair(zp,hp,chi2));
    }
    z.UnlockChildren();
  }

  // Cut on number of solutions.
 
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|m_jj - m_Z| <= %g && |m_jj - m_H| <= %g",xM2j,xM2j);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on cos(theta_ZH).

  TIter nextsol(&solutions);
  ANLPair *sol;
  while ( (sol = (ANLPair*)nextsol()) ) {
    ANL4DVector &z = *(ANL4DVector*)(*sol)[0];
    ANL4DVector &h = *(ANL4DVector*)(*sol)[1];
    Double_t ez = z(0);
    Double_t eh = h(0);
    hEzEh->Fill(ez,eh,1.0);
    Double_t cosz = z.CosTheta();
    Double_t cosh = h.CosTheta();
    hCoszCosh->Fill(cosz,cosh,1.);
    if (TMath::Abs(cosz) > xCoszh || TMath::Abs(cosh) > xCoszh) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_z)| and |cos(theta_h)| <= %g", xCoszh);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Acop.

  nextsol.Reset();
  while ( (sol = (ANLPair*)nextsol()) ) {
    ANL4DVector &z = *(ANL4DVector*)(*sol)[0];
    ANL4DVector &h = *(ANL4DVector*)(*sol)[1];
    Double_t acop = z.Acop(h);
    hAcop->Fill(acop);
    if (acop > xAcop) {
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
  
  // Sort the solutions in the ascending order of chi2 values.

  solutions.Sort();

  // Now store this in ZH2L2JAnalysisBuf.

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
    while ( (jetp = (ANLJet*)nextjet()) ) {
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
  while ( (sol = (ANLPair*)nextsol()) ) {
    if ( nsols++ ) break;
    ANLPair &z = *(ANLPair*)(*sol)[0];
    ANLPair &h = *(ANLPair*)(*sol)[1];
    Double_t chi2   = sol->GetQuality();
    Double_t zmass = z.GetMass();
    hMassz->Fill(zmass);
    Double_t hmass = h.GetMass();
    hMassh->Fill(hmass);
    hChi2->Fill(chi2);
  }

  // Clean up

  CleanUp(&solutions);
  CleanUp(&tracks);

  last->cd();
  return kTRUE;
}
// *****End of Process()***** //


// *****Terminate()***** //
Bool_t ZH2L2JAnalysis::Terminate()
{
  // This function is called at the end of job.
  cout << endl;
  cout << "  =============" << endl;
  cout << "   Cut Summary " << endl;
  cout << "  =============" << endl;
  cout << endl;
  cout << "  ---------------------------------------------------------" <<endl;   
  cout << "  ID   No.Events     Cut Description" << endl;
  cout << "  ---------------------------------------------------------" <<endl;
  Int_t i;
  for ( i = 0; strncmp(&cutName[i][0],"END",4) && i < MAXCUT ; i++ ) {
    printf("  %3d  %10d  : %s\n",i,(int)hStat->GetBinContent(i+1),&cutName[i][0]);
  } 
  cout << "  ---------------------------------------------------------";
  return 0;
}
// *****End of Terminate()***** //
