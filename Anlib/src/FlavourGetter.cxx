//*************************************************************************
//* =====================
//*  FlavourGetter Class
//* =====================
//*
//* (Description)
//*    A flavour of jets getting class
//* (Requires)
//*     class ANLJetFinder
//*     class ANLTrack
//*     class JSFSIMDST, etc
//* (Provides)
//*     class FlavourGetter
//* (Update Recored)
//*    2000/02/13  K.Ikematsu   Original version
//*    2000/04/28  K.Ikemtasu   Modified GetPrimaryHadronPID method
//*    2001/07/20  K.Ikemtasu   Added GetPrimaryHadronSN method
//*    2001/07/20  K.Ikemtasu   Added GetPartonID method
//*    2001/07/26  K.Ikematsu   Added TTL4JFlavourGetter class
//*
//* $Id$
//*************************************************************************
//
#include "FlavourGetter.h"
//_____________________________________________________________________
//  -------------------
//  FlavourGetter Class
//  -------------------
//
ClassImp(FlavourGetter)

//_____________________________________________________________________
//*--
//*  Getters
//*--
Int_t   FlavourGetter::operator()(const ANLJet &jet){

  Int_t nbhad = 0;
  Int_t nchad = 0;
  Int_t nhad = 0;

  TIter next(&jet.GetParticlesInJet());
  ANLTrack *tp;
  while ((tp = (ANLTrack *)next())) {
    Int_t hadpid = TMath::Abs(GetPrimaryHadronPID(*tp));
    if ( hadpid != 0 ) ++nhad;
    Int_t partonid = 0;
    if ( hadpid/1000 > 0 && hadpid/10000 == 0 ) { // Meson-Baryon branch
      partonid = hadpid/1000;
    } else {
      partonid = (hadpid/100)%10;
    }
    if ( partonid == 5 ) ++nbhad;
    if ( partonid == 4 ) ++nchad;
  }

  Double_t blike = (Double_t)nbhad/(Double_t)nhad;
  Double_t clike = (Double_t)nchad/(Double_t)nhad;
  if ( blike >= 0.2 && blike >= clike ) return 3;
  else if ( clike >= 0.2 && clike > blike ) return 2;
  else return 1;
}

//_____________________________________________________________________
Int_t FlavourGetter::GetPrimaryHadronPID(const ANLTrack &t){
  SearchPrimaryHadron(t);
  return fGpid;
}

Int_t FlavourGetter::GetPrimaryHadronSN(const ANLTrack &t){
  SearchPrimaryHadron(t);
  return fGsn;
}

Int_t FlavourGetter::GetPartonID(const ANLTrack &t){
  SearchPrimaryHadron(t);
  return fGmsn;
}

//_____________________________________________________________________
void FlavourGetter::SearchPrimaryHadron(const ANLTrack &t){

  JSFCDCTrack *cdctp = t.GetLTKCLTrack()->GetCDC();
  if (!cdctp) {
    //cerr << "No pointer to JSFCDCTrack !" << endl;
    fGsn  = 0;
    fGpid = 0;
    fGmsn = 0;
    return;
  }

  fGsn = 0;
  fGsn = t.GetLTKCLTrack()->GetCDC()->GetGenID();
  fGpid = 0;
  while ( fGsn > 0 ) {
    JSFGeneratorParticle *g = (JSFGeneratorParticle *)fGen->UncheckedAt(fGsn-1);
    fGpid = g->GetID();
    fGmsn = g->GetMother();

    /*
    cerr << "(PID, S.N, M.S.N) = ("
         << fGpid << ","
         << fGsn  << ","
         << fGmsn << ")" << endl;
    */
    fGsn = fGmsn;
  }
  //cerr << "-- Search ended --" << endl;
}

//_____________________________________________________________________
//  ------------------------
//  TTL4JFlavourGetter Class
//  ------------------------
//
ClassImp(TTL4JFlavourGetter)

//_____________________________________________________________________
//*--
//*  Getters
//*--
Int_t   TTL4JFlavourGetter::operator()(const ANLJet &jet){

  static const Double_t kThetaCut = 10.0;
  static const Double_t kCHadRatioCut = 0.2;

  Int_t ntrkinjet = 0;
  Int_t ntrkfromb = 0;
  Int_t ntrkfromw = 0;
  Int_t ntrkfromchad = 0;

  TIter next(&jet.GetParticlesInJet());
  ANLTrack *tp;
  while ((tp = (ANLTrack *)next())) {
    Int_t hadpid   = TMath::Abs(GetPrimaryHadronPID(*tp));
    Int_t partonid = TMath::Abs(GetPartonID(*tp));
    /*
    cerr << "(HadPID,PartonID) = (" << hadpid << ","
	 << partonid << ")" << endl;
    */
    if ( hadpid != 0 ) ntrkinjet++;

    Int_t flavour = 0;
    if ( hadpid/1000 > 0 && hadpid/10000 == 0 ) { // Meson-Baryon branch
      flavour = hadpid/1000;
    } else {
      flavour = (hadpid/100)%10;
    }

    if ( partonid == 3 ) {
      ntrkfromb++;
    } else if ( partonid == 7 || partonid == 9 ) {
      ntrkfromw++;
    }
    if ( flavour == 4 ) ntrkfromchad++;
  }

  cerr << "R(Fromb,FromW,FromD) = (" << (Double_t)ntrkfromb/ntrkinjet << ","
       << (Double_t)ntrkfromw/ntrkinjet << ","
       << (Double_t)ntrkfromchad/ntrkinjet << ")" << endl;

  JSFSpringParton *spbbar  = (JSFSpringParton *)fSpgen->UncheckedAt(2);
  JSFSpringParton *spb     = (JSFSpringParton *)fSpgen->UncheckedAt(4);
  JSFSpringParton *spupbar = (JSFSpringParton *)fSpgen->UncheckedAt(6);
  JSFSpringParton *spdn    = (JSFSpringParton *)fSpgen->UncheckedAt(7);
  JSFSpringParton *spup    = (JSFSpringParton *)fSpgen->UncheckedAt(8);
  JSFSpringParton *spdnbar = (JSFSpringParton *)fSpgen->UncheckedAt(9);

  ANL4DVector qspbbar(spbbar->GetPV());
  ANL4DVector qspb(spb->GetPV());
  ANL4DVector qspupbar(spupbar->GetPV());
  ANL4DVector qspdn(spdn->GetPV());
  ANL4DVector qspup(spup->GetPV());
  ANL4DVector qspdnbar(spdnbar->GetPV());

  Double_t thetaspbbarj  = qspbbar.GetTheta(jet);
  Double_t thetaspbj     = qspb.GetTheta(jet);
  Double_t thetaspupbarj = qspupbar.GetTheta(jet);
  Double_t thetaspdnj    = qspdn.GetTheta(jet);
  Double_t thetaspupj    = qspup.GetTheta(jet);
  Double_t thetaspdnbarj = qspdnbar.GetTheta(jet);

  if ( spdn->GetID() == 11 || spdn->GetID() == 13 ||  spdn->GetID() == 15 ) {
  cerr << "Th(bbar_j,b_j,up_j,dnbar_j) = ("
       << thetaspbbarj  << ","
       << thetaspbj     << ","
       << thetaspupj    << ","
       << thetaspdnbarj << ")" << endl;
  } else {
  cerr << "Th(bbar_j,b_j,upbar_j,dn_j) = ("
       << thetaspbbarj  << ","
       << thetaspbj     << ","
       << thetaspupbarj << ","
       << thetaspdnj    << ")" << endl;
  }

  if ( (Double_t)ntrkfromb/ntrkinjet > (Double_t)ntrkfromw/ntrkinjet ) {
    if ( TMath::Min(thetaspbbarj,thetaspbj) < kThetaCut ) {
      if ( thetaspbbarj < thetaspbj ) return -3;
      else return 3;
    } else return TMath::Min(thetaspbbarj,thetaspbj);
  } else {
    if ( spdn->GetID() == 11 || spdn->GetID() == 13 ||  spdn->GetID() == 15 ) {
      //cerr << "W- decayed leptonically" << endl;
      if ( (Double_t)ntrkfromchad/ntrkinjet > kCHadRatioCut ) {
	if ( thetaspupj < kThetaCut ) return 2;
	else return thetaspupj;
      } else {
	if ( TMath::Min(thetaspupj,thetaspdnbarj) < kThetaCut ) {
	  if ( thetaspupj < thetaspdnbarj ) return 1;
	  else return -1;
	} else return TMath::Min(thetaspupj,thetaspdnbarj);
      }
    } else {
      //cerr << "W+ decayed leptonically" << endl;
      if ( (Double_t)ntrkfromchad/ntrkinjet > kCHadRatioCut ) {
	if ( thetaspupbarj < kThetaCut ) return -2;
	else return thetaspupbarj;
      } else {
	if ( TMath::Min(thetaspupbarj,thetaspdnj) < kThetaCut ) {
	  if ( thetaspupbarj < thetaspdnj ) return -1;
	  else return 1;
	} else return TMath::Min(thetaspupbarj,thetaspdnj);
      }
    }
  }
}
