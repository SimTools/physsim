//*************************************************************************
//* ========================
//*  ZHH4JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC ZHH data. 
//* (Requires)
//* 	library Anlib
//* 	library ZHHStudy
//* (Provides)
//* 	class ZHH4JAnalysis
//* 	class ZHH4JAnalysisBuf
//* (Usage)
//*   Take a look at Anl.C.
//* (Update Recored)
//*   2008/11/18  K.Fujii	Original version.
//*
//*************************************************************************

#define __CHEAT__

#include "ZHH4JAnalysis.h"
#include "ZHHSpring.h"
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
static const Double_t kMassH   =120.00; // H mass
static const Double_t kSigmaMw =   6.0; // W mass resolution
static const Double_t kSigmaMz =   6.0; // Z mass resolution
static const Double_t kSigmaMh =   6.0; // H mass resolution

typedef enum { kElectron = 11, kMuon = 13 } EPID;

Bool_t gDEBUG = kFALSE;
Int_t  gNgoods  = 0;
const Int_t  kMaxCuts = 100;
      stringstream gCutName[kMaxCuts];
      TH1D  *hStat = 0;

//_____________________________________________________________________
//  --------------------
//  ZHH4JAnalysis Class
//  --------------------
//
//

ClassImp(ZHH4JAnalysis)

ZHH4JAnalysis::ZHH4JAnalysis(const Char_t *name, const Char_t *title)
                : JSFModule(name, title),
                  fNtracksCut(   25),   // No. of tracks
                  fEtrackCut (  0.1),   // track energy
                  fEvisLoCut (  0.0),   // Minimum visible energy
                  fEvisHiCut (999.0),   // Maximum visible energy
                  fPtCut     (  0.0),   // Pt minimum
                  fPlCut     (999.0),   // Pl maximum
                  fElCut     (999.0),   // El maximum
                  fYcutCut   (0.004),   // y_cut to force the event to 4 jets
                  fNjetsCut  (    4),   // No. of jets
                  fEjetCut   (  5.0),   // E_jet minimum
                  fCosjetCut ( 0.99),   // |cos(theta_j)| maximum
                  fM2jCut    ( 18.0),   // |m_jj-m_H| maximum
                  fCoshCut   ( 1.00),   // cos(theta_bh) maximum
                  fMM1Cut    (  70.),   // > mm_HH  
                  fMM2Cut    ( 120.),   // < mm_HH
                  fAcopCut   (  0.0),   // Acoplanarity minimum
                  fBtagNsig  (  1.0),   // Nsig for b-tag
                  fBtagNoffv (    2),   // Noffv for b-tag
                  fBTtagNsig (  3.0),   // Nsig for b-tag (tight)
                  fBTtagNoffv(    2)    // Noffv for b-tag (tight)
{
}

//_____________________________________________________________________
ZHH4JAnalysis::~ZHH4JAnalysis()
{
}

//_____________________________________________________________________
Bool_t ZHH4JAnalysis::Initialize()
{
  //--
  //  Read in Generator info.
  //--
  gJSF->GetInput()->cd("/conf/init");
  ZHHBases *bsp = static_cast<ZHHBases *>(gROOT->FindObject("ZHHBases"));
  cerr << "------------------------------------" << endl
       << " Ecm = " << bsp->GetEcmInit() << " GeV" << endl
       << "------------------------------------" << endl;
  SetEcm(bsp->GetEcmInit()); 

  return 0;
}

//_________________________________________________________
Bool_t ZHH4JAnalysis::Process(Int_t ev)
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
           << "csh1:csh2:eh1:eh2"                                       << ":"
           << "mh1:mh2:mm:acop"                                         << ":"
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

  if (pt > fPtCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Pt < " << fPtCut << ends;
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
  // Find H candidates in a given mass window.
  //--
  
  ANLVTXTagger btag (fBtagNsig,fBtagNoffv);   // (sigma, n-offv tracks) b-tagger
  ANLVTXTagger bttag(fBTtagNsig,fBTtagNoffv); // (sigma, n-offv tracks) tight b-tagger

  TObjArray solutions(10);
  solutions.SetOwner();    // The "solutions" array owns solutions in it.
  ANLPairCombiner h1candidates(jets,jets);
  ANLPair *h1p, *h2p;
  while ((h1p = static_cast<ANLPair *>(h1candidates()))) {
    ANLPair &h1 = *h1p;
#ifdef __CHEAT__
    ANLTaggedJet *j1p = dynamic_cast<ANLTaggedJet *>(h1[0]);
    ANLTaggedJet *j2p = dynamic_cast<ANLTaggedJet *>(h1[1]);
    if (!((j1p->GetTag() == -1 && j2p->GetTag() == -1) ||
          (j1p->GetTag() == -2 && j2p->GetTag() == -2))) continue;
#endif
    if (!btag(*static_cast<ANLJet *>(j1p)) || 
        !btag(*static_cast<ANLJet *>(j2p))) continue;   // double b-tag for b's from t's
    Double_t h1mass = h1().GetMass();
    if (TMath::Abs(h1mass - kMassH) > fM2jCut) continue; // in the Mh window
    h1.LockChildren();
    ANLPairCombiner h2candidates(h1candidates);
    while ((h2p = static_cast<ANLPair *>(h2candidates()))) {
      ANLPair &h2 = *h2p;
#ifdef __CHEAT__
      j1p = dynamic_cast<ANLTaggedJet *>(h2[0]);
      j2p = dynamic_cast<ANLTaggedJet *>(h2[1]);
      if (!((j1p->GetTag() == -1 && j2p->GetTag() == -1) ||
            (j1p->GetTag() == -2 && j2p->GetTag() == -2))) continue;
#endif
      if (h2.IsLocked()) continue;
      if (!btag(*static_cast<ANLJet *>(j1p)) || 
          !btag(*static_cast<ANLJet *>(j2p))) continue;   // double b-tag for b's from t's
      Double_t h2mass = h2().GetMass();
      if (TMath::Abs(h2mass - kMassH) > fM2jCut) continue; // in the Mh window
      if (gDEBUG) { 
        cerr << " M_h1 = " << h1mass << " M_h2 = " << h2mass << endl;
        cerr << " h1p  = " << (void *)h1p 
             << " h2p  = " << (void *)h2p << endl;
        cerr << " h1[0] = " << (void *)h1[0] 
             << " h1[1] = " << (void *)h1[1]
             << " h2[0] = " << (void *)h2[0] 
             << " h2[1] = " << (void *)h2[1] << endl;
      }
      Double_t chi2 = TMath::Power((h1mass - kMassH)/kSigmaMh,2.)
                    + TMath::Power((h2mass - kMassH)/kSigmaMh,2.);
      solutions.Add(new ANLPair(h1p,h2p,chi2));
    }
    h1.UnlockChildren();
  }

  //--
  // Cut on No. of solutions.
  //--

  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|m_jj - m_H| <= " << fM2jCut << ends;
  }
  
  //--
  // Cut on cos(theta_H).
  //--

  TIter nextsol(&solutions);
  ANLPair *solp;
  while ((solp = static_cast<ANLPair *>(nextsol()))) {
    ANLPair &sol = *solp;
    ANL4DVector  &h1 = *static_cast<ANL4DVector *>(sol[0]);
    ANL4DVector  &h2 = *static_cast<ANL4DVector *>(sol[1]);
    Double_t cosh1 = h1.CosTheta();
    Double_t cosh2 = h2.CosTheta();
    if (TMath::Abs(cosh1) > fCoshCut || TMath::Abs(cosh2) > fCoshCut) {
      solutions.Remove(solp);
    }
  }
  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_h)| < " << fCoshCut << ends;
  }

  //--
  // Cut on missing mass.
  //--

  ANL4DVector qcm(fEcm);
  ANL4DVector qmm = qcm -qsum;
  Double_t mm = qmm.GetMass();

  if (mm < fMM1Cut || mm > fMM2Cut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << fMM1Cut << " < mm < " << fMM2Cut << ends;
  }

  //-- 
  // Cut on Acop.
  //-- 

  nextsol.Reset();
  while ((solp = static_cast<ANLPair *>(nextsol()))) {
    ANLPair &sol = *solp;
    ANL4DVector  &h1 = *static_cast<ANL4DVector *>(sol[0]);
    ANL4DVector  &h2 = *static_cast<ANL4DVector *>(sol[1]);
    Double_t acop = h1.Acop(h2);
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
  ANLPair      &h1    = *static_cast<ANLPair *>(sol[0]);
  ANLPair      &h2    = *static_cast<ANLPair *>(sol[1]);
  ANLJet       &j11   = *static_cast<ANLJet  *>(h1[0]);
  ANLJet       &j12   = *static_cast<ANLJet  *>(h1[1]);
  ANLJet       &j21   = *static_cast<ANLJet  *>(h2[0]);
  ANLJet       &j22   = *static_cast<ANLJet  *>(h2[1]);

  //--
  // Calculate helicity angles.
  //--
  TVector3    ez     = TVector3(0., 0., 1.);
  TVector3    eh1z   = h1.Vect().Unit();
  TVector3    eh1x   = eh1z.Cross(ez).Unit();
  TVector3    eh1y   = eh1z.Cross(eh1x);

  TVector3    bsth1  = TVector3(0., 0., h1.Vect().Mag()/h1.E());
  ANL4DVector j11h   = ANL4DVector(j11.E(), j11.Vect()*eh1x,
                                            j11.Vect()*eh1y,
                                            j11.Vect()*eh1z);
  j11h.Boost(-bsth1);
  Double_t    csj11h = j11h.CosTheta();
  Double_t    fij11h = j11h.Phi();

  ANL4DVector j12h   = ANL4DVector(j12.E(), j12.Vect()*eh1x,
                                            j12.Vect()*eh1y,
                                            j12.Vect()*eh1z);
  j12h.Boost(-bsth1);
  Double_t    csj12h = j12h.CosTheta();
  Double_t    fij12h = j12h.Phi();


  TVector3    eh2z   = h2.Vect().Unit();
  TVector3    eh2x   = eh2z.Cross(ez).Unit();
  TVector3    eh2y   = eh2z.Cross(eh2x);

  TVector3    bsth2  = TVector3(0., 0., h2.Vect().Mag()/h2.E());
  ANL4DVector j21h   = ANL4DVector(j21.E(), j21.Vect()*eh2x,
                                            j21.Vect()*eh2y,
                                            j21.Vect()*eh2z);
  j21h.Boost(-bsth2);
  Double_t    csj21h = j21h.CosTheta();
  Double_t    fij21h = j21h.Phi();

  ANL4DVector j22h   = ANL4DVector(j22.E(), j22.Vect()*eh2x,
                                            j22.Vect()*eh2y,
                                            j22.Vect()*eh2z);
  j22h.Boost(-bsth2);
  Double_t    csj22h = j22h.CosTheta();
  Double_t    fij22h = j22h.Phi();

  //--
  // Now store this in the Ntuple.
  //--

  Int_t         nev   = gJSF->GetEventNumber();
  Double_t      ecm   = GetEcm();
  Double_t      chi2  = sol.GetQuality();
  Double_t      mh1   = h1.GetMass();
  Double_t      mh2   = h2.GetMass();
  Double_t      csh1  = h1.CosTheta();
  Double_t      csh2  = h2.CosTheta();
  Double_t      eh1   = h1.E();
  Double_t      eh2   = h2.E();
  Double_t      acop  = h1.Acop(h2);
#ifdef __CHEAT__
  ANLTaggedJet *h1j1p = dynamic_cast<ANLTaggedJet *>(h1[0]);
  ANLTaggedJet *h1j2p = dynamic_cast<ANLTaggedJet *>(h1[1]);
  ANLTaggedJet *h2j1p = dynamic_cast<ANLTaggedJet *>(h2[0]);
  ANLTaggedJet *h2j2p = dynamic_cast<ANLTaggedJet *>(h2[1]);
#if 1
  cerr << "(h1j1,h1j2; h2j1,h2j2) = ("
       << h1j1p->GetTag() << ","
       << h1j2p->GetTag() << ";"
       << h2j1p->GetTag() << ","
       << h2j2p->GetTag() << ";"
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
  data[12] = csh1;
  data[13] = csh2;
  data[14] = eh1;
  data[15] = eh2;
  data[16] = mh1;
  data[17] = mh2;
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
  data[32] = j12.E();
  data[33] = j12.Px();
  data[34] = j12.Py();
  data[35] = j12.Pz();
  data[36] = j21.E();
  data[37] = j21.Px();
  data[38] = j21.Py();
  data[39] = j21.Pz();
  data[40] = j22.E();
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
Bool_t ZHH4JAnalysis::Terminate()
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
