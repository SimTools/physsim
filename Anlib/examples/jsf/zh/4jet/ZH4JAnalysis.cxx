//***************************************************************************
//*  ====================
//*  ZH4JAnalysis Classes
//*  ====================
//*
//*  (Description)
//*	A user analysis class for JLC analyses.
//*	This reads and analyzes MC e+e- -> ZH data.
//*  (Requires)
//*	library Anlib (in physsim-99a-1 if K.Fujii)
//*	library ZHStudy (also in physsim)
//*  (Provides)
//*	class ZH4JAnalysis
//*	class ZH4JAnalysisBuf
//*  (Usage)
//*	...
//*  (Update Record)
//*	18 Nov 1999	A.L.C.Sanchez	Patterned after examples (WW4JAnalysis)
//*					with physsim-99a-1 (still
//*					largely for correction!)
//*	22 Nov 1999	A.L.C.Sanchez	Modified to suit ZH 4jet analysis
//***************************************************************************

#include "ZH4JAnalysis.h"
#include "ANLVTXTagger.h"
#include "TNtupleD.h"
#include <sstream>

static const Double_t kMassH   = 120.0; // H mass
static const Double_t kMassZ   = 91.19;	// Z mass
static const Double_t kSigmaMh =   4.0;	// H mass resolution
static const Double_t kSigmaMz =   4.0;	// Z mass resolution
static const Int_t    kZoneX   =     4;	// No. of X Zones in the Canvas
static const Int_t    kZoneY   =     4;	// No. of Y Zones in the Canvas

Int_t ZH4JAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ---------------------
//  ZH4JAnalysisBuf Class
//  ---------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ClassImp(ZH4JAnalysisBuf)

// Constructors
ZH4JAnalysisBuf::ZH4JAnalysisBuf(const Char_t *name, const Char_t *title,
   ZH4JAnalysis *module) : JSFEventBuf(name, title, (JSFModule*)module) {}

ZH4JAnalysisBuf::ZH4JAnalysisBuf(ZH4JAnalysis *module, const Char_t *name,
   const Char_t *title) : JSFEventBuf(name, title, (JSFModule*)module) {}

ZH4JAnalysisBuf::~ZH4JAnalysisBuf() {}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ------------------
//  ZH4JAnalysis Class
//  ------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ClassImp(ZH4JAnalysis)

// Constructor
ZH4JAnalysis::ZH4JAnalysis(const Char_t *name, const Char_t *title)
               : JSFModule(name, title)
{
  fEventBuf = new ZH4JAnalysisBuf(this);
  SetBufferSize(2000);	// buffer size for event data
  cout << "ZH4JAnalysis is created... fEventBuf is "
       << (Int_t)fEventBuf << endl;
}

// Destructor
ZH4JAnalysis::~ZH4JAnalysis()
{
  cout << "ZH4JAnalysisBuf will be deleted... fEventBuf is "
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}


//####### Functions ###########

// *****CleanUp()***** //
void ZH4JAnalysis::CleanUp(TObjArray *objs)
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
Bool_t ZH4JAnalysis::Initialize()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat       = new TH1F("hStat","Cut Statistics"   ,  20,   0.0,   20.0);
  hNtracks    = new TH1F("hNtracks","No. of Tracks" ,  50,   0.0,  100.0);
  hEvis       = new TH1F("hEvis","Visible Energy"   ,  50,   0.0,  500.0);
  hPt         = new TH1F("hPt","Missing Pt"         ,  50,   0.0,  250.0);
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
  hMassh      = new TH1F("hMassh","M_H Distribution",  50,   0.0,  300.0);
  hMassz      = new TH1F("hMassz","M_Z Distribution",  50,   0.0,  300.0);

  xNtracks  =     25;	// No. of Tracks
  xEtrack   =   0.05;	// track energy
  xEvis     = 100.00;	// Minimum visible energy
  xPt       = 999.00;	// Pt maximum
  xPl       = 999.00;	// Pl maximum
  xYcut     =  0.004;	// y_cut to force the event to 4 jets
  xNjets    =      4;	// No. of Jets
  xEjet     =   5.00;	// E_jet minimum
  xCosjet   =   0.99;	// |cos(theta_j)| maximum
  xCoszh    =   0.99;	// |cos(theta_Z)| and |cos(theta_H)| maximum
  xM2j      =  36.00;	// |m_jj-m_S| maximum , S = Z,H
  xAcop     =  30.00;	// Acoplanarity maximum

  last->cd();
  return 0;
}
// *****End of Initialize()***** //


// *****DrawHist()***** //
void ZH4JAnalysis::DrawHist()
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
Bool_t ZH4JAnalysis::Process(Int_t ev)
{
  // Local copies of ZH4JAnalysisBuf data members.

  Int_t		fNtracks;	// track multiplicity
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
  ZH4JAnalysisBuf *ua	= (ZH4JAnalysisBuf*)fEventBuf;
  ZH4JAnalysisBuf &a	= *ua;

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
  if ( fPt > xPt ) { CleanUp(&tracks); return kFALSE; }
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

  // Find ZH candidates in given mass window.

  ANLVTXTagger btag(3,3);
  
  TObjArray solutions(10);
  ANLPairCombiner zcandidates(jets,jets);
  ANLPair *hp, *zp;
  while ( (zp = (ANLPair*)zcandidates()) ) {
    ANLPair &z = *zp;
    Double_t zmass = z().GetMass();
    if (TMath::Abs(zmass - kMassZ) > xM2j) continue;
    z.LockChildren();
    ANLPairCombiner hcandidates(zcandidates);
    hcandidates.Reset();
    while ( (hp = (ANLPair*)hcandidates()) ) {
      ANLPair &h = *hp;
      if (h.IsLocked()) continue;
      if (!btag(*(ANLJet *)h[0]) || !btag(*(ANLJet *)h[1])) continue;
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
    sprintf(msg,"|m_jj - m_W| <= %g",xM2j);
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

  // Now store this in ZH4JAnalysisBuf.

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

  static TNtupleD *tupp = 0;
  if (!tupp) {
    stringstream tupstr;
    tupstr << "ntracks:evis:pt:pl:ycut:chi2"                            << ":"
           << "npb1:eb1:pb1x:pb1y:pb1z:npb2:eb2:pb2x:pb2y:pb2z"         << ":"
           << "npj1:ej1:pj1x:pj1y:pj1z:npj2:ej2:pj2x:pj2y:pj2z"         << ":"
           << "zmass:hmass:acop"                                        << ":"
           << "cosb1h:phib1h:cosb2h:phib2h:cosj1h:phij1h:cosj2h:phij2h" << ends;

    tupp = new TNtupleD("hEvt","ZHAnalysis",tupstr.str().data());
  }
  Double_t data[100];
 
  nextsol.Reset();
  Int_t nsols = 0;
  while ( (sol = (ANLPair*)nextsol()) ) {
    if ( nsols++ ) break;
    ANLPair &z  = *(ANLPair*)(*sol)[0];
    ANLJet  &j1 = *static_cast<ANLJet  *>(z[0]);
    ANLJet  &j2 = *static_cast<ANLJet  *>(z[1]);
    ANLPair &h  = *(ANLPair*)(*sol)[1];
    ANLJet  &b1 = *static_cast<ANLJet  *>(h[0]);
    ANLJet  &b2 = *static_cast<ANLJet  *>(h[1]);

    Double_t chi2   = sol->GetQuality();
    Double_t zmass = z.GetMass();
    hMassz->Fill(zmass);
    Double_t hmass = h.GetMass();
    hMassh->Fill(hmass);
    hChi2->Fill(chi2);
    //--
    // Calculate helicity angles.
    //--
    TVector3    ez    = TVector3(0., 0., 1.);
    TVector3    ehz   = h.Vect().Unit();
    TVector3    ehx   = ehz.Cross(ez).Unit();
    TVector3    ehy   = ehz.Cross(ehx);

    TVector3    bsth  = TVector3(0., 0., h.Vect().Mag()/h.E());
    ANL4DVector k1h   = ANL4DVector(b1.E(), b1.Vect()*ehx,
                                            b1.Vect()*ehy,
                                            b1.Vect()*ehz);
    k1h.Boost(-bsth);
    Double_t    csb1h = k1h.CosTheta();
    Double_t    fib1h = k1h.Phi();

    ANL4DVector k2h   = ANL4DVector(b2.E(), b2.Vect()*ehx,
                                            b2.Vect()*ehy,
                                            b2.Vect()*ehz);
    k2h.Boost(-bsth);
    Double_t    csb2h = k2h.CosTheta();
    Double_t    fib2h = k2h.Phi();

    TVector3    ezz   = z.Vect().Unit();
    TVector3    ezx   = ezz.Cross(ez).Unit();
    TVector3    ezy   = ezz.Cross(ezx);
    TVector3    bstz  = TVector3(0., 0., z.Vect().Mag()/z.E());
    ANL4DVector p1h   = ANL4DVector(j1.E(), j1.Vect()*ezx,
                                            j1.Vect()*ezy,
                                            j1.Vect()*ezz);
    p1h.Boost(-bstz);
    Double_t    csj1h = p1h.CosTheta();
    Double_t    fij1h = p1h.Phi();

    ANL4DVector p2h   = ANL4DVector(j2.E(), j2.Vect()*ezx,
                                            j2.Vect()*ezy,
                                            j2.Vect()*ezz);
    p2h.Boost(-bstz);
    Double_t    csj2h = p2h.CosTheta();
    Double_t    fij2h = p2h.Phi();

    //--
    // Fill up Ntuple.
    //--
    data[ 0] = fNtracks;
    data[ 1] = fEvis;
    data[ 2] = fPt;
    data[ 3] = fPl;
    data[ 4] = fYcut;
    data[ 5] = sol->GetQuality();
    data[ 6] = b1.GetNparticles();
    data[ 7] = b1()(0);
    data[ 8] = b1()(1);
    data[ 9] = b1()(2);
    data[10] = b1()(3);
    data[11] = b2.GetNparticles();
    data[12] = b2()(0);
    data[13] = b2()(1);
    data[14] = b2()(2);
    data[15] = b2()(3);
    data[16] = j1.GetNparticles();
    data[17] = j1()(0);
    data[18] = j1()(1);
    data[19] = j1()(2);
    data[20] = j1()(3);
    data[21] = j2.GetNparticles();
    data[22] = j2()(0);
    data[23] = j2()(1);
    data[24] = j2()(2);
    data[25] = j2()(3);
    data[26] = z().GetMass();
    data[27] = h().GetMass();
    data[28] = z.Acop(h);
    data[29] = csb1h;
    data[30] = fib1h;
    data[31] = csb2h;
    data[32] = fib2h;
    data[33] = csj1h;
    data[34] = fij1h;
    data[35] = csj2h;
    data[36] = fij2h;
    tupp->Fill(data);
  }

  // Clean up

  CleanUp(&solutions);
  CleanUp(&tracks);

  last->cd();
  return kTRUE;
}
// *****End of Process()***** //


// *****Terminate()***** //
Bool_t ZH4JAnalysis::Terminate()
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
