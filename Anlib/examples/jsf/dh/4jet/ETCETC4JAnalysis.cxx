//*************************************************************************
//* ========================
//*  ETCETC4JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC eta+/- pair data. 
//* (Requires)
//* 	library Anlib
//* 	library ETCETCStudy
//* (Provides)
//* 	class ETCETC4JAnalysis
//* 	class ETCETC4JAnalysisBuf
//* (Usage)
//*   Take a look at Anl.C.
//* (Update Recored)
//*   2008/11/18  K.Fujii	Original version.
//*
//*************************************************************************

//#define __CHEAT__

#include "ETCETC4JAnalysis.h"
#include "ETCETCSpring.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TH1D.h"
#include "JSFSIMDST.h"
#include "Anlib.h"
#include "ANLVTXTagger.h"
#ifdef __CHEAT__
#include "ANLCheatedJetFinder.h"
#endif

#include <sstream>
#include <iomanip>

using namespace std;

static const Double_t kMassW   = 80.00; // W mass
static const Double_t kMassZ   = 91.19; // Z mass
static const Double_t kSigmaMw =   6.0; // W mass resolution
static const Double_t kSigmaMz =   6.0; // W mass resolution

typedef enum { kElectron = 11, kMuon = 13 } EPID;

Bool_t gDEBUG = kFALSE;
Int_t  gNgoods  = 0;
const Int_t  kMaxCuts = 100;
      stringstream gCutName[kMaxCuts];
      TH1D  *hStat = 0;

//_____________________________________________________________________
//  --------------------
//  ETCETC4JAnalysis Class
//  --------------------
//
//

ClassImp(ETCETC4JAnalysis)

ETCETC4JAnalysis::ETCETC4JAnalysis(const Char_t *name, const Char_t *title)
                : JSFModule(name, title),
                  fNtracksCut(   25),   // No. of tracks
                  fEtrackCut (  0.1),   // track energy
                  fEvisLoCut (  0.0),   // Minimum visible energy
                  fEvisHiCut (500.0),   // Maximum visible energy
                  fPtCut     (  0.0),   // Pt minimum
                  fPlCut     (999.0),   // Pl maximum
                  fElCut     (999.0),   // El maximum
                  fYcutCut   (0.004),   // y_cut to force the event to 4 jets
                  fNjetsCut  (    4),   // No. of jets
                  fEjetCut   (  5.0),   // E_jet minimum
                  fCosjetCut ( 0.99),   // |cos(theta_j)| maximum
                  fM2jCut    ( 18.0),   // |m_jj-m_W| maximum
                  fCoswCut   ( 1.00),   // cos(theta_bw) maximum
                  fMM1Cut    (  70.),   // > mm_WW  
                  fMM2Cut    ( 120.),   // < mm_WW
                  fAcopCut   (  0.0),   // Acoplanarity minimum
                  fBtagNsig  (  1.0),   // Nsig for b-tag
                  fBtagNoffv (    2),   // Noffv for b-tag
                  fBTtagNsig (  3.0),   // Nsig for b-tag (tight)
                  fBTtagNoffv(    2)    // Noffv for b-tag (tight)
{
}

//_____________________________________________________________________
ETCETC4JAnalysis::~ETCETC4JAnalysis()
{
}

//_____________________________________________________________________
Bool_t ETCETC4JAnalysis::Initialize()
{
  //--
  //  Read in Generator info.
  //--
  gJSF->GetInput()->cd("/conf/init");
  ETCETCBases *bsp = static_cast<ETCETCBases *>(gROOT->FindObject("ETCETCBases"));
  cerr << "------------------------------------" << endl
       << " Ecm = " << bsp->GetEcmInit() << " GeV" << endl
       << "------------------------------------" << endl;
  SetEcm(bsp->GetEcmInit()); 

  return 0;
}

//_________________________________________________________
Bool_t ETCETC4JAnalysis::Process(Int_t ev)
{
  //--
  // Remember the previous directory.
  //--

  TDirectory *last = gDirectory;
  gFile->cd("/");
  if (!hStat) hStat = new TH1D("hStat","Cut Statistics", 20, 0., 20.);

  static TNtupleD *hEvt = 0;
  if (!hEvt) {
    stringstream tupstr;
    tupstr << "nev:ecm:ntracks:evis:pt:pl:elmax:ycut:chi2"              << ":"
           << "nsols:ejmin:csjmax"                                      << ":"
           << "csw1:csw2:ew1:ew2"                                       << ":"
           << "mw1:mw2:mm:acop"                                         << ":"
           << "csj11h:fij11h:csj12h:fij12h:csj21h:fij21h:csj22h:fij22h" << ":"
	   << "pj1e:pj1x:pj1y:pj1z:pj2e:pj2x:pj2y:pj2z"                 << ":"
	   << "pj3e:pj3x:pj3y:pj3z:pj4e:pj4x:pj4y:pj4z"                 << ends;

    hEvt = new TNtupleD("hEvt", "", tupstr.str().data());
  }

  // ---------------------
  // Analysis starts here.
  // ---------------------

  Double_t selid = -0.5;
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "No Cuts" << ends;
  }

  //--
  // Get event buffer and make combined tracks accessible.
  //--
 
  JSFSIMDST     *sdsp      = static_cast<JSFSIMDST *>(gJSF->FindModule("JSFSIMDST"));
  JSFSIMDSTBuf  *evtp      = static_cast<JSFSIMDSTBuf *>(sdsp->EventBuf());

  Int_t          ntrks   = evtp->GetNLTKCLTracks(); 	// No. of tracks 
  TObjArray     *trks    = evtp->GetLTKCLTracks(); 	// combined tracks

  //--
  // Select good tracks and store them in "TObjArray tracks".
  //--
  //
  ANL4DVector qsum;
  TObjArray tracks(1000);
  tracks.SetOwner();       // the "tracks" array owns tracks in it
  Double_t elmax   = 0.;
  Int_t    ntracks = 0;
  for (Int_t i = 0; i < ntrks; i++) {
    JSFLTKCLTrack *tp = static_cast<JSFLTKCLTrack*>(trks->UncheckedAt(i));
    if (tp->GetE() > fEtrackCut) {
      ANLTrack *qtp = new ANLTrack(tp);
      tracks.Add(qtp); 		// track 4-momentum
      qsum += *qtp;		// total 4-momentum
      ntracks++;
      if (tp->GetType() == kMuon || tp->GetType() == kElectron) {
        if (tp->GetE() > elmax) elmax = tp->GetE();
      }
    }				// *qt stays.
  }
  if (gDEBUG) cerr << "Ntracks = " << ntracks << endl;

  //--
  // Cut on No. of tracks.
  //--

  if (ntracks < fNtracksCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Ntracks >= " << fNtracksCut << ends;
  }

  Double_t evis = qsum(0);      // E_vis
  Double_t pt   = qsum.GetPt(); // P_t
  Double_t pl   = qsum(3);      // P_l

  if (gDEBUG) cerr << "Evis = " << evis << " Pt = " << pt << " Pl = " << pl << endl;
  
  //--
  // Cut on Evis.
  //--

  if (evis < fEvisLoCut || evis > fEvisHiCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << fEvisLoCut << " < Evis < " << fEvisHiCut << ends;
  }

  //--
  // Cut on Pt.
  //--

  if (pt < fPtCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Pt > " << fPtCut << ends;
  }
 
  //--
  // Cut on Pl.
  //--

  if (TMath::Abs(pl) > fPlCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|Pl| <= " << fPlCut << ends;
  }
 
  //--
  // Cut on Elepton.
  //--

  if (elmax > fElCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "El < " << fElCut << ends;
  }

  //--
  // Find jets.
  //--

#ifndef __CHEAT__
  ANLJadeEJetFinder jclust(fYcutCut);
#else
  ANLCheatedJadeEJetFinder jclust(fYcutCut);
  TIter nexttrk(&tracks);
  ANLTrack *trkp;
  while ((trkp = static_cast<ANLTrack *>(nexttrk()))) {
    trkp->SetColorSingletID();
  }
#endif
  jclust.Initialize(tracks);
  jclust.FindJets();
  Double_t ycut  = jclust.GetYcut();
  Int_t    njets = jclust.GetNjets();

  if (gDEBUG) cerr << "Ycut = " << ycut << " Njets = " << njets << endl;

  //--
  // Cut on No. of jets.
  //--

  if (njets < fNjetsCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Njets >= " << fNjetsCut << " for Ycut = " << fYcutCut << ends;
  }

  //--
  // Now force the event to be fNjetsCut.
  //--

  jclust.ForceNJets(fNjetsCut);
  njets = jclust.GetNjets();
  ycut  = jclust.GetYcut();

  if (gDEBUG) cerr << "Ycut = " << ycut << " Njets = " << njets << endl;

  //--
  // Make sure that No. of jets is fNjetsCut.
  //--

  if (njets != fNjetsCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Njets = " << fNjetsCut << ends;
  }

  //--
  // Loop over jets and decide Ejet_min and |cos(theta_j)|_max.
  //--

  TObjArray &jets = jclust.GetJets();
  TIter nextjet(&jets);
  ANLJet *jetp;
  Double_t ejetmin = 999999.;
  Double_t cosjmax = 0.;
  while ((jetp = static_cast<ANLJet *>(nextjet()))) {
    ANLJet &jet = *jetp;
    if (gDEBUG) jet.DebugPrint();
#ifdef __CHEAT__
    if (gDEBUG) {
      ANLTaggedJet *tjetp = dynamic_cast<ANLTaggedJet *>(jetp);
      cerr << " " << tjetp->GetTag();
    }
#endif
    Double_t ejet = jet.E();
    if (ejet < ejetmin) ejetmin = ejet;
    Double_t cosj = jet.CosTheta();
    if (TMath::Abs(cosj) > TMath::Abs(cosjmax)) cosjmax = cosj;
  }

  //--
  // Cut on Ejet_min.
  //--

  if (ejetmin < fEjetCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Ejet >= " << fEjetCut << ends;
  }

  //--
  // Cut on |cos(theta_j)|_max.
  //--

  if (TMath::Abs(cosjmax) > fCosjetCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_j)| <= " << fCosjetCut << ends;
  }

  //--
  // Find W candidates in a given mass window.
  //--
  
  TObjArray solutions(10);
  solutions.SetOwner();    // The "solutions" array owns solutions in it.
  ANLPairCombiner w1candidates(jets,jets);
  ANLPair *w1p, *w2p;
  while ((w1p = static_cast<ANLPair *>(w1candidates()))) {
    ANLPair &w1 = *w1p;
#ifdef __CHEAT__
    ANLTaggedJet *j1p = dynamic_cast<ANLTaggedJet *>(w1[0]);
    ANLTaggedJet *j2p = dynamic_cast<ANLTaggedJet *>(w1[1]);
    if (!((j1p->GetTag() == - 7 && j2p->GetTag() == - 7) ||
          (j1p->GetTag() == - 7 && j2p->GetTag() == - 8) ||
          (j1p->GetTag() == - 8 && j2p->GetTag() == - 7) ||
          (j1p->GetTag() == - 8 && j2p->GetTag() == - 8) ||
          (j1p->GetTag() == - 9 && j2p->GetTag() == - 9) ||
          (j1p->GetTag() == - 9 && j2p->GetTag() == -10) ||
          (j1p->GetTag() == -10 && j2p->GetTag() == - 9) ||
          (j1p->GetTag() == -10 && j2p->GetTag() == -10))) continue;
#endif
    Double_t w1mass = w1().GetMass();
    if (TMath::Abs(w1mass - kMassW) > fM2jCut) continue; // in the Mw window
    w1.LockChildren();
    ANLPairCombiner w2candidates(w1candidates);
    while ((w2p = static_cast<ANLPair *>(w2candidates()))) {
      ANLPair &w2 = *w2p;
#ifdef __CHEAT__
      j1p = dynamic_cast<ANLTaggedJet *>(w2[0]);
      j2p = dynamic_cast<ANLTaggedJet *>(w2[1]);
      if (!((j1p->GetTag() == - 7 && j2p->GetTag() == - 7) ||
            (j1p->GetTag() == - 7 && j2p->GetTag() == - 8) ||
            (j1p->GetTag() == - 8 && j2p->GetTag() == - 7) ||
            (j1p->GetTag() == - 8 && j2p->GetTag() == - 8) ||
            (j1p->GetTag() == - 9 && j2p->GetTag() == - 9) ||
            (j1p->GetTag() == - 9 && j2p->GetTag() == -10) ||
            (j1p->GetTag() == -10 && j2p->GetTag() == - 9) ||
            (j1p->GetTag() == -10 && j2p->GetTag() == -10))) continue;
#endif
      if (w2.IsLocked()) continue;
      Double_t w2mass = w2().GetMass();
      if (TMath::Abs(w2mass - kMassW) > fM2jCut) continue; // in the Mw window
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

  //--
  // Cut on No. of solutions.
  //--

  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|m_jj - m_W| <= " << fM2jCut << ends;
  }
  
  //--
  // Cut on cos(theta_W).
  //--

  TIter nextsol(&solutions);
  ANLPair *solp;
  while ((solp = static_cast<ANLPair *>(nextsol()))) {
    ANLPair &sol = *solp;
    ANL4DVector  &w1 = *static_cast<ANL4DVector *>(sol[0]);
    ANL4DVector  &w2 = *static_cast<ANL4DVector *>(sol[1]);
    Double_t cosw1 = w1.CosTheta();
    Double_t cosw2 = w2.CosTheta();
    if (TMath::Abs(cosw1) > fCoswCut || TMath::Abs(cosw2) > fCoswCut) {
      solutions.Remove(solp);
    }
  }
  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_w)| < " << fCoswCut << ends;
  }

  //--
  // Cut on missing mass.
  //--

  ANL4DVector qcm(fEcm);
  ANL4DVector qmm = qcm -qsum;
  Double_t mm = qmm.GetMass();

  if (mm > fMM1Cut && mm < fMM2Cut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "mm < " << fMM1Cut << " || mm > " << fMM2Cut << ends;
  }

  //-- 
  // Cut on Acop.
  //-- 

  nextsol.Reset();
  while ((solp = static_cast<ANLPair *>(nextsol()))) {
    ANLPair &sol = *solp;
    ANL4DVector  &w1 = *static_cast<ANL4DVector *>(sol[0]);
    ANL4DVector  &w2 = *static_cast<ANL4DVector *>(sol[1]);
    Double_t acop = w1.Acop(w2);
    if (acop < fAcopCut) {
      solutions.Remove(solp);
    }
  }
  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Acop > " << fAcopCut << ends;
  }

  // ----------------------
  // End of event selection
  // ----------------------

  gNgoods++;

  cerr << "--------------------------------------------------" << endl
       << "Event "                   << gJSF->GetEventNumber()
       << ": Number of solutions = " << solutions.GetEntries() << endl
       << "--------------------------------------------------" << endl;

  //--
  // Sort the solutions in the ascending order of chi2 vlues.
  //--

  solutions.Sort();

  if (gDEBUG) {
    Int_t nj = 0;
    nextjet.Reset();
    while ((jetp = static_cast<ANLJet *>(nextjet()))) {
       cerr << "------" << endl
            << "Jet " << ++nj << endl
            << "------" << endl;
       jetp->DebugPrint();
    }
  }

  //--
  // Select the best solution.
  //--
  solutions.Sort();  
  ANLPair &sol = *static_cast<ANLPair *>(solutions.At(0));
  Int_t nsols = solutions.GetEntries();
  ANLPair      &w1    = *static_cast<ANLPair *>(sol[0]);
  ANLPair      &w2    = *static_cast<ANLPair *>(sol[1]);
  ANLJet       &j11   = *static_cast<ANLJet  *>(w1[0]);
  ANLJet       &j12   = *static_cast<ANLJet  *>(w1[1]);
  ANLJet       &j21   = *static_cast<ANLJet  *>(w2[0]);
  ANLJet       &j22   = *static_cast<ANLJet  *>(w2[1]);

  //--
  // Calculate helicity angles.
  //--
  TVector3    ez     = TVector3(0., 0., 1.);
  TVector3    ew1z   = w1.Vect().Unit();
  TVector3    ew1x   = ew1z.Cross(ez).Unit();
  TVector3    ew1y   = ew1z.Cross(ew1x);

  TVector3    bstw1  = TVector3(0., 0., w1.Vect().Mag()/w1.E());
  ANL4DVector j11h   = ANL4DVector(j11.E(), j11.Vect()*ew1x,
                                            j11.Vect()*ew1y,
                                            j11.Vect()*ew1z);
  j11h.Boost(-bstw1);
  Double_t    csj11h = j11h.CosTheta();
  Double_t    fij11h = j11h.Phi();

  ANL4DVector j12h   = ANL4DVector(j12.E(), j12.Vect()*ew1x,
                                            j12.Vect()*ew1y,
                                            j12.Vect()*ew1z);
  j12h.Boost(-bstw1);
  Double_t    csj12h = j12h.CosTheta();
  Double_t    fij12h = j12h.Phi();


  TVector3    ew2z   = w2.Vect().Unit();
  TVector3    ew2x   = ew2z.Cross(ez).Unit();
  TVector3    ew2y   = ew2z.Cross(ew2x);

  TVector3    bstw2  = TVector3(0., 0., w2.Vect().Mag()/w2.E());
  ANL4DVector j21h   = ANL4DVector(j21.E(), j21.Vect()*ew2x,
                                            j21.Vect()*ew2y,
                                            j21.Vect()*ew2z);
  j21h.Boost(-bstw2);
  Double_t    csj21h = j21h.CosTheta();
  Double_t    fij21h = j21h.Phi();

  ANL4DVector j22h   = ANL4DVector(j22.E(), j22.Vect()*ew2x,
                                            j22.Vect()*ew2y,
                                            j22.Vect()*ew2z);
  j22h.Boost(-bstw2);
  Double_t    csj22h = j22h.CosTheta();
  Double_t    fij22h = j22h.Phi();

  //--
  // Now store this in the Ntuple.
  //--

  Int_t         nev   = gJSF->GetEventNumber();
  Double_t      ecm   = GetEcm();
  Double_t      chi2  = sol.GetQuality();
  Double_t      mw1   = w1.GetMass();
  Double_t      mw2   = w2.GetMass();
  Double_t      csw1  = w1.CosTheta();
  Double_t      csw2  = w2.CosTheta();
  Double_t      ew1   = w1.E();
  Double_t      ew2   = w2.E();
  Double_t      acop  = w1.Acop(w2);
#ifdef __CHEAT__
  ANLTaggedJet *w1j1p = dynamic_cast<ANLTaggedJet *>(w1[0]);
  ANLTaggedJet *w1j2p = dynamic_cast<ANLTaggedJet *>(w1[1]);
  ANLTaggedJet *w2j1p = dynamic_cast<ANLTaggedJet *>(w2[0]);
  ANLTaggedJet *w2j2p = dynamic_cast<ANLTaggedJet *>(w2[1]);
#if 1
  cerr << "(w1j1,w1j2; w2j1,w2j2) = ("
       << w1j1p->GetTag() << ","
       << w1j2p->GetTag() << ";"
       << w2j1p->GetTag() << ","
       << w2j2p->GetTag() << ";"
       << endl;
#endif
#endif

  Double_t data[100];
  data[ 0] = nev;
  data[ 1] = ecm;
  data[ 2] = ntracks;
  data[ 3] = evis;
  data[ 4] = pt;
  data[ 5] = pl;
  data[ 6] = elmax;
  data[ 7] = ycut;
  data[ 8] = chi2;
  data[ 9] = nsols;
  data[10] = ejetmin;
  data[11] = cosjmax;
  data[12] = csw1;
  data[13] = csw2;
  data[14] = ew1;
  data[15] = ew2;
  data[16] = mw1;
  data[17] = mw2;
  data[18] = mm;
  data[19] = acop;
  data[20] = csj11h;
  data[21] = fij11h;
  data[22] = csj12h;
  data[23] = fij12h;
  data[24] = csj21h;
  data[25] = fij21h;
  data[26] = csj22h;
  data[27] = fij22h;
  data[28] = j11.E();
  data[29] = j11.Px();
  data[30] = j11.Py();
  data[31] = j11.Pz();
  data[32] = j11.E();
  data[33] = j12.Px();
  data[34] = j12.Py();
  data[35] = j12.Pz();
  data[36] = j21.E();
  data[37] = j21.Px();
  data[38] = j21.Py();
  data[39] = j21.Pz();
  data[40] = j21.E();
  data[41] = j22.Px();
  data[42] = j22.Py();
  data[43] = j22.Pz();

  hEvt->Fill(data);

  last->cd();

  // ------------------------
  //  That's it.
  // ------------------------
  return kTRUE;
}

//_________________________________________________________
Bool_t ETCETC4JAnalysis::Terminate()
{
  // This function is called at the end of job.
  //--
  // Print out cut statistics
  //--
  cerr << endl
       << "  =============" << endl
       << "   Cut Summary " << endl
       << "  =============" << endl
       << endl
       << "  -----------------------------------------------------------" << endl
       << "   ID   No.Events    Cut Description                         " << endl
       << "  -----------------------------------------------------------" << endl;
  for (int id=0; id<kMaxCuts && gCutName[id].str().data()[0]; id++) {
    cerr << "  " << setw( 3) << id
         << "  " << setw(10) << static_cast<int>(hStat->GetBinContent(id+1))
         << "  : " << gCutName[id].str().data() << endl;
  }
  cerr << "  -----------------------------------------------------------" << endl;
  //--
  // That's it!
  //--
  return 0;
}
