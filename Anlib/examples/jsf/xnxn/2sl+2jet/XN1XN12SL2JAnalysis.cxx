//***************************************************************************
//*  ===========================
//*  XN1XN12SL2JAnalysis Classes
//*  ===========================
//*
//*  (Description)
//*	A user analysis class for JLC analyses.
//*	This reads and analyzes MC e+e- -> XN1XN1 -> (stau tau) (stau tau)
//*  (Requires)
//*	library Anlib (in physsim-99a-1 if K.Fujii)
//*	library XN1XN1Study (also in physsim)
//*  (Provides)
//*	class XN1XN12SL2JAnalysis
//*  (Usage)
//*	...
//*  (Update Record)
//*    2010/10/10  K.Fujii	Original version
//***************************************************************************

#include "XN1XN12SL2JAnalysis.h"
#include "XN1XN1Spring.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TH1D.h"
#include "JSFSIMDST.h"
#include "Anlib.h"
#include "ANLVTXTagger.h"

#include <sstream>
#include <iomanip>

using namespace std;

static const Double_t kMassTau  =   1.7; // tau mass
static const Double_t kMassStau = 150.0; // stau mass
static const Double_t kMassX    = 200.0; // X mass
static const Double_t kSigmaMx  =  40.0; // X mass resolution

typedef enum { kElectron = 11, kMuon = 13 } EPID;

Bool_t gDEBUG = kFALSE;
//Bool_t gDEBUG = kTRUE;
Int_t  gNgoods  = 0;
const Int_t  kMaxCuts = 100;
      stringstream gCutName[kMaxCuts];
      TH1D  *hStat = 0;

//_____________________________________________________________________
//  ------------------------
//  XN1XN12SL2JAnalysis Class
//  ------------------------

ClassImp(XN1XN12SL2JAnalysis)

XN1XN12SL2JAnalysis::XN1XN12SL2JAnalysis(const Char_t *name, const Char_t *title)
               : JSFModule(name, title),
                  fNtracksCut(    3),   // No. of tracks
                  fEtrackCut (  0.1),   // track energy
                  fEvisLoCut (  0.0),   // Minimum visible energy
                  fEvisHiCut (999.0),   // Maximum visible energy
                  fPtCut     (999.0),   // Pt minimum
                  fPlCut     (999.0),   // Pl maximum
                  fEleptonCut( 50.0),   // El maximum
		  fCosConeCut( 0.94),   // Cos(cone) cut
		  fEconeCut  ( 2.00),   // Econe cut
                  fYcutCut   (0.004),   // y_cut to force the event to 4 jets
                  fNjetsCut  (    2),   // No. of jets
                  fEjetCut   (  5.0),   // E_jet minimum
                  fCosjetCut ( 0.99),   // |cos(theta_j)| maximum
                  fM2jCut    ( 40.0),   // |m_jj-m_H| maximum
                  fCosxCut   ( 1.00),   // cos(theta_bh) maximum
                  fMM1Cut    (   0.),   // > mm_HH  
                  fMM2Cut    ( 999.),   // < mm_HH
                  fAcopCut   (  0.0)    // Acoplanarity minimum
{
}

XN1XN12SL2JAnalysis::~XN1XN12SL2JAnalysis()
{
}

Bool_t XN1XN12SL2JAnalysis::Initialize()
{
#if 0
  //--
  //  Read in Generator info.
  //--
  gJSF->GetInput()->cd("/conf/init");
  XN1XN1Bases *bsp = static_cast<XN1XN1Bases *>(gROOT->FindObject("XN1XN1HBases"));
  SetEcm(bsp->GetRoots()); 
#endif
  cerr << "------------------------------------" << endl
       << " Ecm = " << GetEcm() << " GeV"        << endl
       << "------------------------------------" << endl;

  return 0;
}


Bool_t XN1XN12SL2JAnalysis::Process(Int_t ev)
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
           << "csx1:csx2:ex1:ex2"                                       << ":"
           << "mx1:mx2:mm:acop"                                         << ":"
	   << "pj1e:pj1x:pj1y:pj1z:pj2e:pj2x:pj2y:pj2z"                 << ":"
	   << "pj3e:pj3x:pj3y:pj3z:pj4e:pj4x:pj4y:pj4z"                 << ":"
	   << "cstau1h:fitau1h:csstau1h:fistau1h"                       << ":"
	   << "cstau2h:fitau2h:csstau2h:fistau2h"                       << ends;

    hEvt = new TNtupleD("hEvt", "", tupstr.str().data());
  }

  //--
  // Analysis starts here.
  //--
  Float_t selid = -0.5;
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
  // Select good tracks
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
    }
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
  // Find Isolated sLepton.
  //--
  TObjArray lptracks(20);
  Int_t nlptracks = 0;
  TIter nexttrk (&tracks);
  ANLTrack *trkp;
  Double_t lpcharge = 0.;
  while ((trkp = static_cast<ANLTrack *>(nexttrk()))) {
    ANLTrack &trk = *trkp;
    if (! trk.IsLepton()) continue;
    Double_t elepton = trk.E();
    if (elepton < fEleptonCut) continue;
    Double_t econe = trk.GetConeEnergy(fCosConeCut, &tracks);
    if (econe <= fEconeCut) {
      lptracks.Add(trkp);
      trk.SetE(TMath::Sqrt(trk.GetMag2()+kMassStau*kMassStau));
      nlptracks++;
      trk.Lock();
      lpcharge += trk.GetCharge();
    }
  }

  if (gDEBUG) cerr << "nlptracks = " << nlptracks << endl;

  //--
  // Require only two Isolated sLeptons.
  //--
  if (nlptracks != 2 || lpcharge) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Nslptracks = 2" << ends;
  }

  //--
  // Find jets.
  //--
  ANLJadeEJetFinder jclust(fYcutCut);
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
  // Make sure the No. of Jets is fNjetsCut.
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
  // Unlock lepton tracks.
  //--
  TIter nextlptrk(&lptracks);
  while ((trkp = static_cast<ANLTrack *>(nextlptrk()))) trkp->Unlock();

  //--
  // Find XN1XN1 candidates in given mass window.
  //--
  TObjArray solutions(10);
  ANLPairCombiner x1candidates(lptracks,jets);
  ANLPair *x1p, *x2p;
  while ((x1p = static_cast<ANLPair *>(x1candidates()))) {
    ANLPair &x1 = *x1p;
#if 0
    Double_t x1mass = x1().GetMass();
#else
    ANL4DVector pvstau1 = x1(0);
    ANL4DVector pvtau1  = x1(1);
    Double_t etau1   = fEcm/2 - pvstau1(0);
    Double_t apvtau1 = TMath::Sqrt((etau1-kMassTau)*(etau1+kMassTau));
    Double_t apv1    = pvtau1.GetMag();
    pvtau1(0)  = etau1;
    pvtau1(1) *= apvtau1/apv1;
    pvtau1(2) *= apvtau1/apv1;
    pvtau1(3) *= apvtau1/apv1;
    Double_t x1mass = (pvtau1 + pvstau1).GetMass();
#endif
    if (TMath::Abs(x1mass - kMassX) > fM2jCut) continue;
    x1.LockChildren();
    ANLPairCombiner x2candidates(x1candidates);
    while ((x2p = static_cast<ANLPair*>(x2candidates()))) {
      ANLPair &x2 = *x2p;
      if (x2.IsLocked()) continue;
#if 0
      Double_t x2mass = x2().GetMass();
#else
      ANL4DVector pvstau2 = x2(0);
      ANL4DVector pvtau2  = x2(1);
      Double_t etau2   = fEcm/2 - pvstau2(0);
      Double_t apvtau2 = TMath::Sqrt((etau2-kMassTau)*(etau2+kMassTau));
      Double_t apv2    = pvtau2.GetMag();
      pvtau2(0)  = etau2;
      pvtau2(1) *= apvtau2/apv2;
      pvtau2(2) *= apvtau2/apv2;
      pvtau2(3) *= apvtau2/apv2;
      Double_t x2mass = (pvtau2 + pvstau2).GetMass();
#endif
      if (TMath::Abs(x2mass - kMassX) > fM2jCut) continue;
      if (gDEBUG) {
        cerr << " M_X1 = " << x1mass << " M_X2 = " << x2mass << endl;
        cerr << " x1p  = " << (void*)x1p
             << " x2p  = " << (void*)x2p << endl;
        cerr << " x1[0] = " << (void*)x1[0]
             << " x1[1] = " << (void*)x1[1]
             << " x2[0] = " << (void*)x2[0]
             << " x2[1] = " << (void*)x2[1] << endl;
      }
      Double_t chi2 = TMath::Power((x1mass - kMassX)/kSigmaMx,2.)
                    + TMath::Power((x2mass - kMassX)/kSigmaMx,2.);
      solutions.Add(new ANLPair(x1p,x2p,chi2));
    }
    x1.UnlockChildren();
  }

  //--
  // Cut on number of solutions.
  //--
  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|m_jj - m_H| <= " << fM2jCut << ends;
  }

  //--
  // Cut on cos(theta_XN1XN1).
  //--
  TIter nextsol(&solutions);
  ANLPair *solp;
  while ((solp = static_cast<ANLPair *>(nextsol()))) {
    ANLPair &sol = *solp;
    ANL4DVector  &x1 = *static_cast<ANL4DVector *>(sol[0]);
    ANL4DVector  &x2 = *static_cast<ANL4DVector *>(sol[1]);
    Double_t cosx1 = x1.CosTheta();
    Double_t cosx2 = x2.CosTheta();
    if (TMath::Abs(cosx1) > fCosxCut || TMath::Abs(cosx2) > fCosxCut) {
      solutions.Remove(solp);
    }
  }
  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_h)| < " << fCosxCut << ends;
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
    ANL4DVector  &x1 = *static_cast<ANL4DVector *>(sol[0]);
    ANL4DVector  &x2 = *static_cast<ANL4DVector *>(sol[1]);
    Double_t acop = x1.Acop(x2);
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
  // Sort the solutions in the ascending order of chi2 values.
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
  ANLPair      &x1    = *static_cast<ANLPair *>(sol[0]);
  ANLPair      &x2    = *static_cast<ANLPair *>(sol[1]);
  ANLJet       &j11   = *static_cast<ANLJet  *>(x1[0]);
  ANLJet       &j12   = *static_cast<ANLJet  *>(x1[1]);
  ANLJet       &j21   = *static_cast<ANLJet  *>(x2[0]);
  ANLJet       &j22   = *static_cast<ANLJet  *>(x2[1]);

  ANL4DVector   pvstau1 = j11;
  ANL4DVector   pvtau1  = j12;
  Double_t      etau1   = fEcm/2 - pvstau1(0);
  Double_t      apvtau1 = TMath::Sqrt((etau1-kMassTau)*(etau1+kMassTau));
  Double_t      apv1    = pvtau1.GetMag();
  pvtau1(0)  = etau1;
  pvtau1(1) *= apvtau1/apv1;
  pvtau1(2) *= apvtau1/apv1;
  pvtau1(3) *= apvtau1/apv1;
  ANL4DVector pvx1 = pvstau1 + pvtau1;

  ANL4DVector   pvstau2 = j21;
  ANL4DVector   pvtau2  = j22;
  Double_t      etau2   = fEcm/2 - pvstau2(0);
  Double_t      apvtau2 = TMath::Sqrt((etau2-kMassTau)*(etau2+kMassTau));
  Double_t      apv2    = pvtau2.GetMag();
  pvtau2(0)  = etau2;
  pvtau2(1) *= apvtau2/apv2;
  pvtau2(2) *= apvtau2/apv2;
  pvtau2(3) *= apvtau2/apv2;
  ANL4DVector pvx2 = pvstau2 + pvtau2;

  //--
  // Now store this in the Ntuple.
  //--

  Int_t         nev   = gJSF->GetEventNumber();
  Double_t      ecm   = GetEcm();
  Double_t      chi2  = sol.GetQuality();
  Double_t      mx1   = pvx1.GetMass();
  Double_t      mx2   = pvx2.GetMass();
  Double_t      csx1  = pvx1.CosTheta();
  Double_t      csx2  = pvx2.CosTheta();
  Double_t      ex1   = pvx1.E();
  Double_t      ex2   = pvx2.E();
  Double_t      acop  = pvx1.Acop(x2);
  //--
  // Calculate helicity angles.
  //--
  TVector3    ez     = TVector3(0., 0., 1.);
  TVector3    ex1z   = pvx1.Vect().Unit();
  TVector3    ex1x   = ex1z.Cross(ez).Unit();
  TVector3    ex1y   = ex1z.Cross(ex1x);

  TVector3    bstx1 = TVector3(0., 0., pvx1.Vect().Mag()/pvx1.E());
  ANL4DVector k1h   = ANL4DVector(pvtau1.E(), pvtau1.Vect()*ex1x,
                                              pvtau1.Vect()*ex1y,
                                              pvtau1.Vect()*ex1z);
  k1h.Boost(-bstx1);
  Double_t    cstau1h = k1h.CosTheta();
  Double_t    fitau1h = k1h.Phi();

  ANL4DVector k2h   = ANL4DVector(pvstau1.E(), pvstau1.Vect()*ex1x,
                                               pvstau1.Vect()*ex1y,
                                               pvstau1.Vect()*ex1z);
  k2h.Boost(-bstx1);
  Double_t    csstau1h = k2h.CosTheta();
  Double_t    fistau1h = k2h.Phi();

  TVector3    ex2z   = pvx2.Vect().Unit();
  TVector3    ex2x   = ex2z.Cross(ez).Unit();
  TVector3    ex2y   = ex2z.Cross(ex2x);
  TVector3    bstx2  = TVector3(0., 0., pvx2.Vect().Mag()/pvx2.E());
  ANL4DVector p1h    = ANL4DVector(pvtau2.E(), pvtau2.Vect()*ex2x,
                                                pvtau2.Vect()*ex2y,
                                                pvtau2.Vect()*ex2z);
  p1h.Boost(-bstx2);
  Double_t    cstau2h = p1h.CosTheta();
  Double_t    fitau2h = p1h.Phi();

  ANL4DVector p2h   = ANL4DVector(pvstau2.E(), pvstau2.Vect()*ex2x,
                                               pvstau2.Vect()*ex2y,
                                               pvstau2.Vect()*ex2z);
  p2h.Boost(-bstx2);
  Double_t    csstau2h = p2h.CosTheta();
  Double_t    fistau2h = p2h.Phi();

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
  data[12] = csx1;
  data[13] = csx2;
  data[14] = ex1;
  data[15] = ex2;
  data[16] = mx1;
  data[17] = mx2;
  data[18] = mm;
  data[19] = acop;
  data[20] = pvstau1.E();
  data[21] = pvstau1.Px();
  data[22] = pvstau1.Py();
  data[23] = pvstau1.Pz();
  data[24] = pvtau1.E();
  data[25] = pvtau1.Px();
  data[26] = pvtau1.Py();
  data[27] = pvtau1.Pz();
  data[28] = pvstau2.E();
  data[29] = pvstau2.Px();
  data[30] = pvstau2.Py();
  data[31] = pvstau2.Pz();
  data[32] = pvtau2.E();
  data[33] = pvtau2.Px();
  data[34] = pvtau2.Py();
  data[35] = pvtau2.Pz();

  data[36] = cstau1h;
  data[37] = fitau1h;
  data[38] = csstau1h;
  data[39] = fistau1h;
  data[40] = cstau2h;
  data[41] = fitau2h;
  data[42] = csstau2h;
  data[43] = fistau2h;

  hEvt->Fill(data);

  last->cd();

  // ------------------------
  //  That's it.
  // ------------------------
  return kTRUE;
}


Bool_t XN1XN12SL2JAnalysis::Terminate()
{
  // This function is called at the end of job.
  cout << endl
       << "  =============" << endl
       << "   Cut Summary " << endl
       << "  =============" << endl
       << endl
       << "  ---------------------------------------------------------" << endl
       << "  ID   No.Events     Cut Description                       " << endl
       << "  ---------------------------------------------------------" << endl;
  for (int id=0; id<kMaxCuts && gCutName[id].str().data()[0]; id++) {
    cout << "  " << setw( 3) << id
         << "  " << setw(10) << static_cast<int>(hStat->GetBinContent(id+1))
         << "  : " << gCutName[id].str().data() << endl;
  }
  cout << "  ---------------------------------------------------------" << endl;

  //--
  // That's it!
  //--
  return 0;
}
