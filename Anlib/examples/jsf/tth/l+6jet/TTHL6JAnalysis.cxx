//*************************************************************************
//* ========================
//*  TTHL6JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC ttbar data to select lepton + 6-jet events.
//* (Requires)
//* 	library Anlib
//* 	library TTHStudy
//* (Provides)
//* 	class TTHL6JAnalysis
//* 	class TTHL6JAnalysisBuf
//* (Usage)
//*   Take a look at anlL6J.C.
//* (Update Recored)
//*   2002/08/15  K.Fujii	Derived from TTL4JAnalysis.cxx.
///*
//*************************************************************************

#define __CHEAT__

#include "TTHL6JAnalysis.h"
#include "TTHSpring.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TH1D.h"
#include "JSFSIMDST.h"
#include "Anlib.h"
#include "ANLVTXTagger.h"
#include "ANLCheatedJetFinder.h"

#include <sstream>
#include <iomanip>

static const Double_t kMassW   = 80.00; 	// W mass
static const Double_t kMassZ   = 91.19; 	// Z mass
static const Double_t kMasst   = 175.0; 	// top mass
static const Double_t kMassH   = 120.0;		// Higgs mass
#if 0
static const Double_t kSigmaMw =   8.0; 	// W mass resolution
static const Double_t kSigmaMz =   8.0; 	// Z mass resolution
static const Double_t kSigmaMt =  15.0; 	// top mass resolution
static const Double_t kSigmaMh =   8.0; 	// H mass resolution
#else
static const Double_t kSigmaMw =   6.0; 	// W mass resolution
static const Double_t kSigmaMz =   6.0; 	// Z mass resolution
static const Double_t kSigmaMt =  14.0; 	// top mass resolution
static const Double_t kSigmaMh =  11.0; 	// H mass resolution
#endif

Bool_t gDEBUG = kFALSE;
Int_t  gNgoods  = 0;
const Int_t  kMaxCuts = 100;
      stringstream gCutName[kMaxCuts];
      TH1D  *hStat = 0;

//_____________________________________________________________________
//  -------------------
//  TTHL6JAnalysis Class
//  -------------------
//
//

ClassImp(TTHL6JAnalysis)

TTHL6JAnalysis::TTHL6JAnalysis(const Char_t *name,
		               const Char_t *title)
              : JSFModule(name, title),
                fNtracksCut(   25),   // No. of tracks
                fEtrackCut (  0.1),   // track energy
                fEvisCut   ( 80.0),   // Minimum visible energy
                fPtCut     (999.0),   // Pt maximum
                fPlCut     (999.0),   // Pl maximum
                fEleptonCut( 10.0),    // Elepton mimimum
                fCosConeCut( 0.94),    // cos(theta_cone)
                fYcutCut   (0.004),   // y_cut to force the event to 4 jets
                fNjetsCut  (    8),   // No. of jets
                fEjetCut   (  5.0),   // E_jet minimum
                fCosjetCut ( 0.99),   // |cos(theta_j)| maximum
                fCosbwCut  ( 1.00),   // cos(theta_bw) maximum
                fM2jCut    ( 18.0),   // |m_jj-m_W| maximum
                fM3jCut    ( 24.0),   // |m_3j-m_t| maximum
                fThrustCut (  0.9),   // Thrust maximum
                fBtagNsig  (  1.0),   // Nsig for b-tag
                fBtagNoffv (    2),   // Noffv for b-tag
                fBTtagNsig (  3.0),   // Nsig for b-tag (tight)
                fBTtagNoffv(    2)    // Noffv for b-tag (tight)
{
}

//_____________________________________________________________________
TTHL6JAnalysis::~TTHL6JAnalysis()
{
}

//_____________________________________________________________________
Bool_t TTHL6JAnalysis::Initialize()
{
  //--
  //  Read in Generator info.
  //--
  gJSF->GetInput()->cd("/conf/init");
  TTHBases *bsp = static_cast<TTHBases *>(gROOT->FindObject("TTHBases"));
  cerr << "------------------------------------" << endl
       << " Ecm = " << bsp->GetRoots() << " GeV" << endl
       << "------------------------------------" << endl;
  SetEcm(bsp->GetRoots());

  return kTRUE;
}


//_________________________________________________________
Bool_t TTHL6JAnalysis::Process(Int_t ev)
{
  using namespace std;

  // Remember the previous directory.

  TDirectory *last = gDirectory;
  gFile->cd("/");
  if (!hStat) hStat = new TH1D("hStat","Cut Statistics", 20, 0., 20.);

  static TNtupleD *hEvt = 0;
  if (!hEvt) {
    stringstream tupstr;
    tupstr << "nev:ecm:ntracks:evis:pt:pl:ycut:chi2"                    << ":"
           << "nsols:ejmin:csjmax"                                      << ":"
           << "csbw1:csbw2"                                             << ":"
           << "mw1:mw2:mt1:mt2:mh:mtt:thrust"                           << ":"
           << ends;

    hEvt = new TNtupleD("hEvt", "", tupstr.str().data());
  }
#if 0
  static TNtupleD *hLep = 0;
  if (!hLep) {
    stringstream tupstr;
    tupstr << "nev:ecm:ntracks:evis:pt:pl"                              << ":"
           << "elep:econe"                                              << ":"
           << ends;

    hLep = new TNtupleD("hLep", "", tupstr.str().data());
  }
#endif

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

  JSFSIMDST     *sdsp    = static_cast<JSFSIMDST *>(gJSF->FindModule("JSFSIMDST"));
  JSFSIMDSTBuf  *evtp    = static_cast<JSFSIMDSTBuf *>(sdsp->EventBuf());

  Int_t          ntrks   = evtp->GetNLTKCLTracks(); 	// No. of tracks 
  TObjArray     *trks    = evtp->GetLTKCLTracks(); 	// combined tracks

  //--
  // Select good tracks and store them in "TObjArray tracks".
  //--

  ANL4DVector qsum;
  TObjArray tracks(1000);
  tracks.SetOwner();
  Int_t ntracks = 0;
  for (Int_t i = 0; i < ntrks; i++) {
    JSFLTKCLTrack *tp = static_cast<JSFLTKCLTrack*>(trks->UncheckedAt(i));
    if (tp->GetE() > fEtrackCut) {
      ANLTrack *qtp = new ANLTrack(tp);
      tracks.Add(qtp); 		// track 4-momentum
      qsum += *qtp;		// total 4-momentum
      ntracks++;
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

  Double_t evis = qsum(0);		// E_vis
  Double_t pt   = qsum.GetPt();		// P_t
  Double_t pl   = qsum(3);		// P_l

  if (gDEBUG) cerr << "Evis = " << evis << " Pt = " << pt << " Pl = " << pl << endl;

  //--
  // Cut on Evis.
  //--

  if (evis < fEvisCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Evis >= " << fEvisCut << ends;
  }

  //--
  // Cut on Pt.
  //--

  if (pt > fPtCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Pt <= " << fPtCut << ends;
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
  // Find Isolated Lepton.
  //--

  TObjArray lptracks(20);
  Int_t nlptracks = 0;
  TIter nexttrk (&tracks);
  ANLTrack *trkp;
  Double_t lpcharge = 0.;
  while ((trkp = static_cast<ANLTrack *>(nexttrk()))) {
    ANLTrack &trk = *trkp;
    if (!trk.IsLepton()) continue;
    Double_t elepton = trk.E();
#if 0
    Int_t    nev   = gJSF->GetEventNumber();
    Double_t ecm   = GetEcm();
    Double_t econe = trk.GetConeEnergy(fCosConeCut, &tracks);
    hLep->Fill(nev,ecm,ntracks,evis,pt,pl,elepton,econe);
    if (elepton < fEleptonCut) continue;
#else
    if (elepton < fEleptonCut) continue;
    Double_t econe = trk.GetConeEnergy(fCosConeCut, &tracks);
#endif
    if (econe <= fEconeCut) {
      lptracks.Add(trkp);
      nlptracks++;
      trk.Lock();
      lpcharge = trk.GetCharge();
    }
  }

  //--
  // Require only one Isolated Lepton.
  //--

  if (nlptracks != 1) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Nlptracks = 1" << ends;
  }

  //--
  // Find jets.
  //--

#ifndef __CHEAT__
  ANLJadeEJetFinder jclust(fYcutCut);
#else
  ANLCheatedJadeEJetFinder jclust(fYcutCut);
  nexttrk.Reset();
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
  // Find W and top candidates in given mass windows.
  //--

  ANL4DVector qcm(GetEcm());
  ANL4DVector qnu = qcm - qsum;

#ifndef __CHEAT__
  ANLJet lpjet(lptracks);
  ANLJet nujet;
#else
  ANLTaggedJet lpjet(lptracks);
  ANLTaggedJet nujet;
#endif
  nujet.Add(&qnu);

  ANLPair w1(&lpjet, &nujet);
  ANLPair *w1p = &w1;

  ANLVTXTagger btag (fBtagNsig,fBtagNoffv);   // (sigma, n-offv tracks) b-tagger
  ANLVTXTagger bttag(fBTtagNsig,fBTtagNoffv); // (sigma, n-offv tracks) tight b-tagger

  TObjArray solutions(10);
  solutions.SetOwner();    // The "solutions" array owns solutions in it.

  ANLPairCombiner w2candidates(jets,jets);
  ANLPair *w2p, *bbp, *hp;
  while ((w2p = static_cast<ANLPair *>(w2candidates()))) {
    ANLPair &w2 = *w2p;
#ifdef __CHEAT__
    ANLTaggedJet *j1p = dynamic_cast<ANLTaggedJet *>(w2[0]);
    ANLTaggedJet *j2p = dynamic_cast<ANLTaggedJet *>(w2[1]);
    if (!((j1p->GetTag() == -8 && j2p->GetTag() == -8) ||
          (j1p->GetTag() == -8 && j2p->GetTag() == -9) ||
          (j1p->GetTag() == -9 && j2p->GetTag() == -8) ||
          (j1p->GetTag() == -9 && j2p->GetTag() == -9) ||
          (j1p->GetTag() == -10 && j2p->GetTag() == -10) ||
          (j1p->GetTag() == -10 && j2p->GetTag() == -11) ||
          (j1p->GetTag() == -11 && j2p->GetTag() == -10) ||
          (j1p->GetTag() == -11 && j2p->GetTag() == -11))) continue;
#endif
    if (bttag(*static_cast<ANLJet *>(w2[0])) || 
        bttag(*static_cast<ANLJet *>(w2[1]))) continue;   // anti-b-tag for W daughters
    Double_t w2mass = w2.GetMass();
    if (TMath::Abs(w2mass - kMassW) > fM2jCut) continue;
    w2.LockChildren();
    ANLPairCombiner bbcandidates(w2candidates);
    bbcandidates.Reset();
    while ((bbp = static_cast<ANLPair *>(bbcandidates()))) {
      ANLPair &bb = *bbp;
      if (bb.IsLocked()) continue;
#ifdef __CHEAT__
      j1p = dynamic_cast<ANLTaggedJet *>(bb[0]);
      j2p = dynamic_cast<ANLTaggedJet *>(bb[1]);
      if (!((j1p->GetTag() == -4 && j2p->GetTag() == -4) ||
            (j1p->GetTag() == -4 && j2p->GetTag() == -6) ||
            (j1p->GetTag() == -6 && j2p->GetTag() == -4) ||
            (j1p->GetTag() == -6 && j2p->GetTag() == -6))) continue;
#endif
      if (!btag(*static_cast<ANLJet *>(bb[0])) || 
          !btag(*static_cast<ANLJet *>(bb[1]))) continue;   // double b-tag for b's from t's
      bb.LockChildren();
      for (Int_t i = 0; i < 2; i++) {
	ANL4DVector *b1p = static_cast<ANL4DVector *>(bb[i]);
	ANL4DVector *b2p = static_cast<ANL4DVector *>(bb[1-i]);
	ANLPair *bw1p = new ANLPair(b1p,w1p);
	ANLPair *bw2p = new ANLPair(b2p,w2p);
	ANLPair &bw2  = *bw2p;
	// Double_t t1mass = bw1.GetMass();
	Double_t t2mass = bw2.GetMass();
	if (TMath::Abs(t2mass - kMasst) > fM3jCut) {
	  delete bw1p;
	  delete bw2p;
	  continue;
	}
	ANLPairCombiner hcandidates(bbcandidates);
	hcandidates.Reset();
	Bool_t ok = kFALSE;
        while ((hp = static_cast<ANLPair *>(hcandidates()))) { 
	  ANLPair &h = *hp;
#ifdef __CHEAT__
          j1p = dynamic_cast<ANLTaggedJet *>(h[0]);
          j2p = dynamic_cast<ANLTaggedJet *>(h[1]);
          if (j1p->GetTag() != j2p->GetTag()) continue;
          if (j1p->GetTag() != -3) continue;
#endif
	  Double_t hmass = h.GetMass();
	  if (h.IsLocked() ||
	     !btag(*static_cast<ANLJet *>(h[0])) || // double b-tag
	     !btag(*static_cast<ANLJet *>(h[1])) || // h daughters (h->bb)
	      TMath::Abs(hmass - kMassH) > fM2jCut) continue; // h in the Mh window
	  Double_t chi2 = TMath::Power((w2mass - kMassW)/kSigmaMw,2.)
                        + TMath::Power((t2mass - kMasst)/kSigmaMt,2.)
                        + TMath::Power((hmass  - kMassH)/kSigmaMh,2.);
	  ANLPair *ttp = new ANLPair(bw1p,bw2p);
	  ANLPair *hpp = new ANLPair(h);
	  solutions.Add(new ANLPair(ttp,hpp,chi2));
	  ok = kTRUE;
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

  //--
  // Cut on No. of solutions.
  //--

  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|m_jj - m_W| <= " << fM2jCut << " && "
	                   << "|m_3j - m_t| <= " << fM3jCut << " && "
	                   << "|m_jj - m_h| <= " << fM2jCut << ends;
  }

  //--
  // Cut on cos(theta_bW).
  //--

  TIter nextsol(&solutions);
  ANLPair *solp;
  while ((solp = static_cast<ANLPair *>(nextsol()))) {
    ANLPair &sol = *solp;
    ANLPair &tt  = *static_cast<ANLPair *>(sol[0]);
    ANLPair &bw1 = *static_cast<ANLPair *>(tt[0]);
    ANLPair &bw2 = *static_cast<ANLPair *>(tt[1]);
    ANLJet  &b1  = *static_cast<ANLJet  *>(bw1[0]);
    ANLPair &w1  = *static_cast<ANLPair *>(bw1[1]);
    ANLJet  &b2  = *static_cast<ANLJet  *>(bw2[0]);
    ANLPair &w2  = *static_cast<ANLPair *>(bw2[1]);
    Double_t cosbw1 = b1.CosTheta(w1);
    Double_t cosbw2 = b2.CosTheta(w2);
    if (cosbw1 > fCosbwCut || cosbw2 > fCosbwCut) {
      solutions.Remove(solp);
      ANLPair &tt = *static_cast<ANLPair *>(sol[0]);
      tt.Delete();
      sol.Delete();
    }
  }
  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_bw)| > " << fCosbwCut << ends;
  }

  //--
  // Cut on Thrust.
  //--

  ANLEventShape eshape;
  eshape.Initialize(tracks);
  Double_t thrust = eshape.GetThrust();
  if (thrust > fThrustCut) {
  	nextsol.Reset();
  	while ((solp = static_cast<ANLPair *>(nextsol()))) {
          ANLPair *ttp = static_cast<ANLPair *>((*solp)[0]);
          ttp ->Delete();
          solp->Delete();
        }
  	return kFALSE;
  }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Thrust <= " << fThrustCut << ends;
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

  Int_t nsols = solutions.GetEntries();
  nextsol.Reset();
  Int_t nsol = 0;
  while ((solp = static_cast<ANLPair *>(nextsol()))) {
    if (nsol++) break;				// choose the best
    ANLPair &sol = *solp;
    ANLPair &tt  = *static_cast<ANLPair *>(sol[0]);
    ANLPair &h   = *static_cast<ANLPair *>(sol[1]);
    ANLPair &bw1 = *static_cast<ANLPair *>(tt[0]);
    ANLPair &bw2 = *static_cast<ANLPair *>(tt[1]);
    ANLJet  &b1  = *static_cast<ANLJet  *>(bw1[0]);
    ANLPair &w1  = *static_cast<ANLPair *>(bw1[1]);
    ANLJet  &b2  = *static_cast<ANLJet  *>(bw2[0]);
    ANLPair &w2  = *static_cast<ANLPair *>(bw2[1]);
    Int_t    nev   = gJSF->GetEventNumber();
    Double_t ecm   = GetEcm();
    Double_t chi2  = sol.GetQuality();
    Double_t mw1   = w1.GetMass();
    Double_t mw2   = w2.GetMass();
    Double_t mt1   = bw1.GetMass();
    Double_t mt2   = bw2.GetMass();
    Double_t mh    = h.GetMass();
    Double_t mtt   = tt.GetMass();
    Double_t csbw1 = b1.CosTheta(w1);
    Double_t csbw2 = b2.CosTheta(w2);
#ifdef __CHEAT__
    ANLTaggedJet *b1p   = dynamic_cast<ANLTaggedJet *>(&b1);
    ANLTaggedJet *w1j1p = dynamic_cast<ANLTaggedJet *>(w1[0]);
    ANLTaggedJet *w1j2p = dynamic_cast<ANLTaggedJet *>(w1[1]);
    ANLTaggedJet *b2p   = dynamic_cast<ANLTaggedJet *>(&b2);
    ANLTaggedJet *w2j1p = dynamic_cast<ANLTaggedJet *>(w2[0]);
    ANLTaggedJet *w2j2p = dynamic_cast<ANLTaggedJet *>(w2[1]);
    ANLTaggedJet *hj1p  = dynamic_cast<ANLTaggedJet *>(h[0]);
    ANLTaggedJet *hj2p  = dynamic_cast<ANLTaggedJet *>(h[1]);
#if 1
    cerr << "(b1,w1j1,w1j2; b2,w2j1,w2j2; hj1, hj2) = ("
         << b1p  ->GetTag() << ","
         << w1j1p->GetTag() << ","
         << w1j2p->GetTag() << ","
         << b2p  ->GetTag() << ","
         << w2j1p->GetTag() << ","
         << w2j2p->GetTag() << ";"
	 << hj1p ->GetTag() << ","
	 << hj2p ->GetTag() << ")"
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
    data[ 6] = ycut;
    data[ 7] = chi2;
    data[ 8] = nsols;
    data[ 9] = ejetmin;
    data[10] = cosjmax;
    data[11] = csbw1;
    data[12] = csbw2;
    data[13] = mw1;
    data[14] = mw2;
    data[15] = mt1;
    data[16] = mt2;
    data[17] = mh;
    data[18] = mtt;
    data[19] = thrust;

    hEvt->Fill(data);
  }

  //--
  // Clean up.
  //--

  nextsol.Reset();
  while ((solp = static_cast<ANLPair *>(nextsol()))) {
    ANLPair *ttp = static_cast<ANLPair *>((*solp)[0]);
    ttp ->Delete();
    solp->Delete();
  }

  last->cd();

  // ------------------------
  //  That's it.
  // ------------------------
  return kTRUE;
}

//_________________________________________________________
Bool_t TTHL6JAnalysis::Terminate()
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
