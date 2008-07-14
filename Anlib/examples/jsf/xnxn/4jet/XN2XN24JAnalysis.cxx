//*************************************************************************
//* ========================
//*  XN2XN24JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC chargino pair data. 
//* (Requires)
//* 	library Anlib
//* 	library XN2XN2Study
//* (Provides)
//* 	class XN2XN24JAnalysis
//* 	class XN2XN24JAnalysisBuf
//* (Usage)
//*   Take a look at Anl.C.
//* (Update Recored)
//*   1999/08/01  K.Fujii	Original version.
//*
//*************************************************************************
//
#include "XN2XN24JAnalysis.h"
#include "ANLTrack.h"

static const Double_t kMassW   = 80.00; // W mass
static const Double_t kMassZ   = 91.19; // Z mass
static const Double_t kSigmaMw =   4.0; // W mass resolution
static const Double_t kSigmaMz =   4.0; // W mass resolution
static const Int_t    kZoneX   = 4;	// No. X Zones in the Canvas
static const Int_t    kZoneY   = 4;	// No. Y Zones in the Canvas

Int_t XN2XN24JAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

typedef enum { kElectron = 11, kMuon = 13 } EPID;

//_____________________________________________________________________
//  -----------------------
//  XN2XN24JAnalysisBuf Class
//  -----------------------
//
//
ClassImp(XN2XN24JAnalysisBuf)

//_________________________________________________________
XN2XN24JAnalysisBuf::XN2XN24JAnalysisBuf(const Char_t *name, const Char_t *title, 
   XN2XN24JAnalysis *module) : JSFEventBuf(name, title, (JSFModule*)module) {}

//_________________________________________________________
XN2XN24JAnalysisBuf::XN2XN24JAnalysisBuf(XN2XN24JAnalysis *module, const Char_t *name,
   const Char_t *title) : JSFEventBuf(name, title, (JSFModule*)module) {}


//_____________________________________________________________________
//  --------------------
//  XN2XN24JAnalysis Class
//  --------------------
//
//

ClassImp(XN2XN24JAnalysis)

XN2XN24JAnalysis::XN2XN24JAnalysis(const Char_t *name, const Char_t *title)
	       : JSFModule(name, title) 
{
  fEventBuf = new XN2XN24JAnalysisBuf(this);
  SetBufferSize(2000);  // buffer size for event data.
  cout << "XN2XN24JAnalysisBuf is created...fEventBuf is " 
       << (Int_t)fEventBuf << endl;
}

//_____________________________________________________________________
XN2XN24JAnalysis::~XN2XN24JAnalysis()
{
  cout << "XN2XN24JAnalysisBuf will be deleted...fEventBuf is " 
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}

//_____________________________________________________________________
void XN2XN24JAnalysis::CleanUp(TObjArray *objs)
{
  TIter next(objs);
  TObject *obj;
  while ((obj = next())) {
  	objs->Remove(obj);
  	delete obj;
  }
}

//_____________________________________________________________________
Bool_t XN2XN24JAnalysis::Initialize()
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
  hEz1Ez2     = new TH2F("hEz1Ez2","(E_z1,E_z2)" ,  
  				 200,  0.0, 200.0, 200,   0.0, 200.0);
  hCosz1Cosz2 = new TH2F("hCosz1Cosz2","(cos_z1,cos_z2)",
				  50, -1.0,  +1.0,  50,  -1.0,  +1.0);
  hMz1Mz2     = new TH2F("hMz1Mz2","(m_z1,m_z2)" , 
  				  60, 50.0, 110.0,  60,  50.0, 110.0);
  hEvisPl     = new TH2F("hEvisPl","(Evis,Pl)"   , 
  				  60,  0.0, 600.0,  50,-100.0,+100.0);
  hMM         = new TH1F("hMM","mm_zz"           ,  80, 100.0, 500.0);
  hAcop       = new TH1F("hAcop","Acoplanarity"  ,  90,   0.0, 180.0);
  hEz         = new TH1F("hEz","E_z"             , 400,   0.0, 200.0);

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
  xCosz     =    0.90;	// |cos(theta_j)| maximum
  xM2jLo    =   10.00;	// |m_jj-m_Z| maximum
  xM2jHi    =   20.00;	// |m_jj-m_Z| maximum
  xMM1      =   70.00;   // missing mass cut against ZZ
  xMM2      =  120.00;   // missing mass cut against ZZ
  xAcop     =   30.00;	// Acoplanarity minimum

  last->cd();
  return 0;
}

//_________________________________________________________
void XN2XN24JAnalysis::DrawHist()
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
  cHist->cd(10);	hEz1Ez2->Draw();
  cHist->cd(11);	hCosz1Cosz2->Draw();
  cHist->cd(12);	hMz1Mz2->Draw();
  cHist->cd(13);	hEvisPl->Draw();
  cHist->cd(14);	hAcop->Draw();

  cHist->Update();
  
  last->cd();
}

//_________________________________________________________
Bool_t XN2XN24JAnalysis::Process(Int_t ev)
{
  // Local copies of XN2XN24JAnalysisBuf data members.

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
  XN2XN24JAnalysisBuf *ua  = (XN2XN24JAnalysisBuf *)fEventBuf;
  XN2XN24JAnalysisBuf &a   = *ua;

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

  // Find Z candidates in a given mass window.
  
  TObjArray solutions(10);
  ANLPairCombiner z1candidates(jets,jets);
  ANLPair *z1p, *z2p;
  while ((z1p = (ANLPair *)z1candidates())) {
    ANLPair &z1 = *z1p;
    Double_t z1mass = z1().GetMass();
    if (z1mass < kMassZ-xM2jLo || z1mass > kMassZ+xM2jHi) continue;
    z1.LockChildren();
    ANLPairCombiner z2candidates(z1candidates);
    while ((z2p = (ANLPair *)z2candidates())) {
      ANLPair &z2 = *z2p;
      if (z2.IsLocked()) continue;
      Double_t z2mass = z2().GetMass();
      if (z2mass < kMassZ-xM2jLo || z2mass > kMassZ+xM2jHi) continue;
      if (gDEBUG) { 
        cerr << " M_z1 = " << z1mass << " M_z2 = " << z2mass << endl;
        cerr << " z1p  = " << (void *)z1p 
             << " z2p  = " << (void *)z2p << endl;
        cerr << " z1[0] = " << (void *)z1[0] 
             << " z1[1] = " << (void *)z1[1]
             << " z2[0] = " << (void *)z2[0] 
             << " z2[1] = " << (void *)z2[1] << endl;
      }
      Double_t chi2 = TMath::Power((z1mass - kMassZ)/kSigmaMz,2.)
                    + TMath::Power((z2mass - kMassZ)/kSigmaMz,2.);
      solutions.Add(new ANLPair(z1p,z2p,chi2));
      // hMz1Mz2->Fill(z1mass,z2mass,1.0);
    }
    z1.UnlockChildren();
  }
  
  // Cut on No. of solutions.

  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"m_Z - %g GeV < m_jj < m_Z + %g GeV", xM2jLo, xM2jHi);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  // Cut on cos(theta_Z).

  TIter nextsol(&solutions);
  ANLPair *sol;
  while ((sol = (ANLPair *)nextsol())) {
    ANL4DVector  &z1 = *(ANL4DVector *)(*sol)[0];
    ANL4DVector  &z2 = *(ANL4DVector *)(*sol)[1];
    Double_t cosz1 = z1.CosTheta();
    Double_t cosz2 = z2.CosTheta();
    hCosz1Cosz2->Fill(cosz1,cosz2,1.0);
    if (TMath::Abs(cosz1) > xCosz || TMath::Abs(cosz2) > xCosz) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_z)| <= %g",xCosz);
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
    sprintf(msg,"mm_ZZ <= %g GeV or mm_ZZ >= %g GeV",xMM1, xMM2);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  // Cut on Acop.

  nextsol.Reset();
  while ((sol = (ANLPair *)nextsol())) {
    ANL4DVector  &z1 = *(ANL4DVector *)(*sol)[0];
    ANL4DVector  &z2 = *(ANL4DVector *)(*sol)[1];
    Double_t acop = z1.Acop(z2);
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

  // Now store this in XN2XN24JAnalysisBuf.
  
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
    ANLPair &z1 = *(ANLPair *)(*sol)[0];
    ANLPair &z2 = *(ANLPair *)(*sol)[1];
    Double_t chi2   = sol->GetQuality();
    Double_t z1mass = z1.GetMass();
    Double_t z2mass = z2.GetMass();
#if 0
    Double_t ez1 = z1()(0);
    Double_t ez2 = z2()(0);
#else
    Double_t ez1 = TMath::Sqrt(kMassZ*kMassZ + z1().GetMag2());
    Double_t ez2 = TMath::Sqrt(kMassZ*kMassZ + z2().GetMag2());
#endif
          hChi2->Fill(chi2);
          hMz1Mz2->Fill(z1mass,z2mass,1.0);
          hEz1Ez2->Fill(ez1,ez2,1.0);
          hEz->Fill(ez1,1.0);
          hEz->Fill(ez2,1.0);
  }
  
  // Clean up
  
  CleanUp(&solutions);
  CleanUp(&tracks);
  
  last->cd();
  return kTRUE;
}

//_________________________________________________________
Bool_t XN2XN24JAnalysis::Terminate()
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
