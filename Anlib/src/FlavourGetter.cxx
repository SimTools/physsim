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
//*  Setters
//*--
void FlavourGetter::SetDebug(Bool_t flag) {
  if ( flag != fDEBUG ) {
    fDEBUG = flag;
  }
}

//_____________________________________________________________________
//*--
//*  Getters
//*--
Int_t FlavourGetter::operator()(const ANLJet &jet){

  static const Double_t kBHadRatioCut   = 0.3;
  static const Double_t kCHadRatioCut   = 0.2;
  static const Double_t kPriHadRatioCut = 0.9;

  Int_t ntrkinjet = 0;
  Int_t ntrkfrombhad = 0;
  Int_t ntrkfromchad = 0;

  TIter next(&jet.GetParticlesInJet());
  ANLTrack *tp;
  while ((tp = (ANLTrack *)next())) {
    Int_t hadpid = TMath::Abs(GetPrimaryHadronPID(*tp));
    if ( hadpid != 0 ) ntrkinjet++;

    Int_t flavour = 0;
    if ( hadpid/1000 > 0 && hadpid/10000 == 0 ) { // Meson-Baryon branch
      flavour = hadpid/1000;
    } else {
      flavour = (hadpid/100)%10;
    }
    if ( flavour == 5 ) ntrkfrombhad++;
    if ( flavour == 4 ) ntrkfromchad++;
  }

  Double_t bratio   = (Double_t)ntrkfrombhad/(Double_t)ntrkinjet;
  Double_t cratio   = (Double_t)ntrkfromchad/(Double_t)ntrkinjet;
  Double_t priratio = (Double_t)(ntrkinjet-ntrkfrombhad-ntrkfromchad)
                      /(Double_t)ntrkinjet;

  if (fDEBUG) cerr << "R(Frombhad,Fromchad,Fromprihad) = ("
		   << bratio << ","
		   << cratio << ","
		   << priratio << ")" << endl;

  if ( bratio >= kBHadRatioCut && cratio < kCHadRatioCut ) return 3;
  else if ( bratio < kBHadRatioCut && cratio >= kCHadRatioCut ) return 2;
  else if ( priratio >= kPriHadRatioCut ) return 1;
  else return 0;
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

    if (fDEBUG) cerr << "(PID, S.N, M.S.N) = ("
		     << fGpid << ","
		     << fGsn  << ","
		     << fGmsn << ")" << endl;

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
//*  Setters
//*--
void TTL4JFlavourGetter::SetDebug(Bool_t flag) {
  if ( flag != fDEBUG ) {
    fDEBUG = flag;
  }
}

//_____________________________________________________________________
//*--
//*  Getters
//*--
Int_t TTL4JFlavourGetter::operator()(const ANLJet &jet){

  static const Double_t kThetaCut = 10.0;
  //static const Double_t kCHadRatioCut = 0.25;
  //temporary treatment
  static const Double_t kCHadRatioCut = 0.1;

  Int_t ntrkinjet = 0;
  Int_t ntrkfromb = 0;
  Int_t ntrkfromw = 0;
  Int_t ntrkfromchad = 0;

  TIter next(&jet.GetParticlesInJet());
  ANLTrack *tp;
  while ((tp = (ANLTrack *)next())) {
    Int_t hadpid   = TMath::Abs(GetPrimaryHadronPID(*tp));
    Int_t partonid = TMath::Abs(GetPartonID(*tp));
    if ( hadpid != 0 ) ntrkinjet++;

    Int_t flavour = 0;
    if ( hadpid/1000 > 0 && hadpid/10000 == 0 ) { // Meson-Baryon branch
      flavour = hadpid/1000;
    } else {
      flavour = (hadpid/100)%10;
    }

    if ( partonid == 3 ) ntrkfromb++;
    else if ( partonid == 7 || partonid == 9 ) ntrkfromw++;
    if ( flavour == 4 ) ntrkfromchad++;
  }

  if (fDEBUG)
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
    if (fDEBUG)
      cerr << "W- decayed leptonically" << endl
	   << "Th(bbar_j,b_j,up_j,dnbar_j) = ("
	   << thetaspbbarj  << ","
	   << thetaspbj     << ","
	   << thetaspupj    << ","
	   << thetaspdnbarj << ")" << endl;
  } else {
    if (fDEBUG)
      cerr << "W+ decayed leptonically" << endl
	   << "Th(bbar_j,b_j,upbar_j,dn_j) = ("
	   << thetaspbbarj  << ","
	   << thetaspbj     << ","
	   << thetaspupbarj << ","
	   << thetaspdnj    << ")" << endl;
  }

  if ( (Double_t)ntrkfromb/ntrkinjet > (Double_t)ntrkfromw/ntrkinjet &&
       (Double_t)ntrkfromchad/ntrkinjet < kCHadRatioCut ) {
    if ( TMath::Min(thetaspbbarj,thetaspbj) < kThetaCut ) {
      if ( thetaspbbarj < thetaspbj ) return -3;
      else return 3;
    } else {
      if (fDEBUG) cerr << "Th(bottom_j) = "
		       << TMath::Min(thetaspbbarj,thetaspbj) << endl;
      return 0;
    }
  } else {
    if ( spdn->GetID() == 11 || spdn->GetID() == 13 ||  spdn->GetID() == 15 ) {
      if ( (Double_t)ntrkfromchad/ntrkinjet >= kCHadRatioCut ) {
	if ( thetaspupj < kThetaCut ) return 2;
	else {
	  if (fDEBUG) cerr << "Th(charm_j) = " << thetaspupj << endl;
	  return 0;
	}
      } else {
	if ( TMath::Min(thetaspupj,thetaspdnbarj) < kThetaCut ) {
	  if ( thetaspupj < thetaspdnbarj ) return 1;
	  else return -1;
	} else {
	  if (fDEBUG) cerr << "Th(light_j) = "
			   << TMath::Min(thetaspupj,thetaspdnbarj) << endl;
	  return 0;
	}
      }
    } else {
      if ( (Double_t)ntrkfromchad/ntrkinjet > kCHadRatioCut ) {
	if ( thetaspupbarj < kThetaCut ) return -2;
	else {
	  if (fDEBUG) cerr << "Th(charm_j) = " << thetaspupbarj << endl;
	  return 0;
	}
      } else {
	if ( TMath::Min(thetaspupbarj,thetaspdnj) < kThetaCut ) {
	  if ( thetaspupbarj < thetaspdnj ) return -1;
	  else return 1;
	} else {
	  if (fDEBUG) cerr << "Th(light_j) = "
			   << TMath::Min(thetaspupbarj,thetaspdnj) << endl;
	  return 0;
	}
      }
    }
  }
}
