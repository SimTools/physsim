//*************************************************************************
//* ========================
//*  ETRETI4JAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyze MC eta_R eta_I data. 
//* (Requires)
//* 	library Anlib
//* 	library ETRETIStudy
//* (Provides)
//* 	class ETRETI4JAnalysis
//* 	class ETRETI4JAnalysisBuf
//* (Usage)
//*   Take a look at Anl.C.
//* (Update Recored)
//*   2008/11/18  K.Fujii	Original version.
//*
//*************************************************************************

//#define __CHEAT__

#include "ETRETI4JAnalysis.h"
#include "ETRETISpring.h"
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
//  ETRETI4JAnalysis Class
//  --------------------
//
//

ClassImp(ETRETI4JAnalysis)

ETRETI4JAnalysis::ETRETI4JAnalysis(const Char_t *name, const Char_t *title)
                : JSFModule(name, title),
                  fNtracksCut(   25),   // No. of tracks
                  fEtrackCut (  0.1),   // track energy
                  fEvisLoCut (  0.0),   // Minimum visible energy
                  fEvisHiCut (500.0),   // Maximum visible energy
                  fPtCut     (  0.0),   // Pt minimum
                  fPlCut     (999.0),   // Pl maximum
                  fElCut     (999.0),   // El maximum
                  fYcutCut   (0.004),   // y_cut to force the event to 4 jets
                  fNjetsCut  (    2),   // No. of jets
                  fEjetCut   (  5.0),   // E_jet minimum
                  fCosjetCut ( 0.99),   // |cos(theta_j)| maximum
                  fM2jCut    ( 18.0),   // |m_jj-m_W| maximum
                  fCoszCut   ( 1.00),   // cos(theta_bw) maximum
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
ETRETI4JAnalysis::~ETRETI4JAnalysis()
{
}

//_____________________________________________________________________
Bool_t ETRETI4JAnalysis::Initialize()
{
  //--
  //  Read in Generator info.
  //--
  gJSF->GetInput()->cd("/conf/init");
  ETRETIBases *bsp = static_cast<ETRETIBases *>(gROOT->FindObject("ETRETIBases"));
  cerr << "------------------------------------" << endl
       << " Ecm = " << bsp->GetEcmInit() << " GeV" << endl
       << "------------------------------------" << endl;
  SetEcm(bsp->GetEcmInit()); 

  return 0;
}

//_________________________________________________________
Bool_t ETRETI4JAnalysis::Process(Int_t ev)
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
           << "ejmin:csjmax"                                            << ":"
           << "csz:ez"                                                  << ":"
           << "mz:mm"                                                   << ":"
           << "csj1h:fij1h:csj2h:fij2h"                                 << ":"
	   << "pj1e:pj1x:pj1y:pj1z:pj2e:pj2x:pj2y:pj2z"                 << ends;

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
  // Find Z candidates in a given mass window.
  //--
  
  ANLPairCombiner zcandidates(jets,jets);
  ANLPair z(jets[0], jets[1]);
  Double_t zmass = z().GetMass();
  Double_t chi2 = TMath::Power((zmass - kMassW)/kSigmaMw,2.);
  if (TMath::Abs(zmass - kMassZ) > fM2jCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|m_jj - m_Z| <= " << fM2jCut << ends;
  }

  //--
  // Cut on cos(theta_Z).
  //--

  Double_t cosz = z.CosTheta();
  if (TMath::Abs(cosz) > fCoszCut ) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_z)| < " << fCoszCut << ends;
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

  // ----------------------
  // End of event selection
  // ----------------------

  gNgoods++;

  cerr << "--------------------------------------------------" << endl
       << "Event "                   << gJSF->GetEventNumber() << endl
       << "--------------------------------------------------" << endl;

  ANLJet       &j1   = *static_cast<ANLJet  *>(z[0]);
  ANLJet       &j2   = *static_cast<ANLJet  *>(z[1]);

  //--
  // Calculate helicity angles.
  //--
  TVector3    ez    = TVector3(0., 0., 1.);
  TVector3    ezz   = z.Vect().Unit();
  TVector3    ezx   = ezz.Cross(ez).Unit();
  TVector3    ezy   = ezz.Cross(ezx);

  TVector3    bstz  = TVector3(0., 0., z.Vect().Mag()/z.E());
  ANL4DVector j1h   = ANL4DVector(j1.E(), j1.Vect()*ezx,
                                          j1.Vect()*ezy,
                                          j1.Vect()*ezz);
  j1h.Boost(-bstz);
  Double_t    csj1h = j1h.CosTheta();
  Double_t    fij1h = j1h.Phi();

  ANL4DVector j2h   = ANL4DVector(j2.E(), j2.Vect()*ezx,
                                          j2.Vect()*ezy,
                                          j2.Vect()*ezz);
  j2h.Boost(-bstz);
  Double_t    csj2h = j2h.CosTheta();
  Double_t    fij2h = j2h.Phi();

  //--
  // Now store this in the Ntuple.
  //--

  Double_t data[100];
  data[ 0] = gJSF->GetEventNumber();
  data[ 1] = GetEcm();
  data[ 2] = ntracks;
  data[ 3] = evis;
  data[ 4] = pt;
  data[ 5] = pl;
  data[ 6] = elmax;
  data[ 7] = ycut;
  data[ 8] = chi2;
  data[ 9] = ejetmin;
  data[10] = cosjmax;
  data[11] = z.CosTheta();
  data[12] = z.E();
  data[13] = z.GetMass();
  data[14] = mm;
  data[15] = csj1h;
  data[16] = fij1h;
  data[17] = csj2h;
  data[18] = fij2h;
  data[19] = j1.E();
  data[20] = j1.Px();
  data[21] = j1.Py();
  data[22] = j1.Pz();
  data[23] = j2.E();
  data[24] = j2.Px();
  data[25] = j2.Py();
  data[26] = j2.Pz();

  hEvt->Fill(data);

  last->cd();

  // ------------------------
  //  That's it.
  // ------------------------
  return kTRUE;
}

//_________________________________________________________
Bool_t ETRETI4JAnalysis::Terminate()
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
