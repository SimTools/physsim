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
//*    2001/07/31  K.Ikematsu   Supported to search generator particles
//*                             contributing to the EM cluster (gamma)
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
void FlavourGetter::SetData(const ANLJet &jet) {

  fPIDGen.Clear();
  fSNGen.Clear();
  fMSNGen.Clear();

  TIter next(&jet.GetParticlesInJet());
  ANLTrack *tp;
  while ((tp = (ANLTrack *)next())) {
    SearchPrimaryHadron(*tp);
  }

  fPIDGen.SetOwner();  // SetOwner() method only enabled
  fSNGen.SetOwner();   // after adding contents
  fMSNGen.SetOwner();
}

void FlavourGetter::SetDebug(Bool_t flag) {
  if ( flag != fDEBUG ) {
    fDEBUG = flag;
  }
}

//_____________________________________________________________________
//*--
//*  Getters
//*--
Int_t FlavourGetter::operator()(const ANLJet &jet) {

  static const Double_t kBHadRatioCut   = 0.25;
  static const Double_t kCHadRatioCut   = 0.15;
  static const Double_t kPriHadRatioCut = 0.9;

  Int_t ntrkinjet = 0;
  Int_t ntrkfrombhad = 0;
  Int_t ntrkfromchad = 0;

  SetData(jet);

  TIter nexthadpid(&fPIDGen);
  TObjNum *pidp;
  while ((pidp = (TObjNum *)nexthadpid())) {
    // fPIDGen : Not Yet contained primary hadron's PIDs contributing to
    //           the HDC cluster from neutral combined tracks.
    Int_t hadpid = TMath::Abs(pidp->GetNum());
    if (fDEBUG) cerr << "FlavourGetter::operator() : PID = " << hadpid << endl;
    ntrkinjet++;

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
TObjArray & FlavourGetter::GetPrimaryHadronPID() { return fPIDGen; }

TObjArray & FlavourGetter::GetPrimaryHadronSN() { return fSNGen; }

TObjArray & FlavourGetter::GetPartonID() { return fMSNGen; }

//_____________________________________________________________________
void FlavourGetter::SearchPrimaryHadron(const ANLTrack &t) {

  JSFLTKCLTrack *ctp = t.GetLTKCLTrack();

  if (fDEBUG) {
    if ( ctp->GetType() == 1 ) {
      cerr << "Combined track type : pure gamma" << endl;
    } else if ( ctp->GetType() == 2 ) {
      cerr << "Combined track type : gamma in mixed EMC cluster" << endl;
    } else if ( ctp->GetType() == 3 ) {
      cerr << "Combined track type : pure nutral hadron" << endl;
    } else if ( ctp->GetType() == 4 ) {
      cerr << "Combined track type : hadron in mixed HDC cluster" << endl;
    } else if ( ctp->GetType() == 5 ) {
      cerr << "Combined track type : pure charged hadron" << endl;
    } else if ( ctp->GetType() == 11 ) {
      cerr << "Combined track type : electron candidate" << endl;
    } else if ( ctp->GetType() == 13 ) {
      cerr << "Combined track type : muon candidate" << endl;
    } else if ( ctp->GetType() == 6 ) {
      cerr << "Combined track type : unmatched track" << endl;
    } else {
      cerr << "illegal track type !!!" << endl;
    }
    cerr << "  CDCEntries   = " << ctp->GetCDCEntries() << endl;
    cerr << "  EMGenEntries = " << ctp->GetEMGenEntries() << endl;
  }

  Int_t nemgen = ctp->GetEMGenEntries();
  if ( nemgen > 0 ) {
    if (fDEBUG) cerr << "# of neutral tracks in EMC = " << nemgen << endl;
    for (Int_t i = 0; i < nemgen; i++ ) {
      if (fDEBUG) cerr << "i = " << i << endl;
      Int_t gpid = 0;
      Int_t gsn  = 0;
      Int_t gmsn = ctp->GetEMGenAt(i)->GetSerial();
      while ( gmsn >= 0 ) {
	JSFGeneratorParticle *g = (JSFGeneratorParticle *)fGen->UncheckedAt(gmsn-1);
	gpid = g->GetID();
	gsn  = g->GetSerial();
	gmsn = g->GetMother();

	if (fDEBUG) cerr << "(PID, S.N, M.S.N) = ("
			 << gpid << ","
			 << gsn  << ","
			 << gmsn << ")" << endl;
      }
      if (fDEBUG) cerr << "-- Search ended --" << endl;
      TObjNum *gpidp = new TObjNum(gpid);
      TObjNum *gsnp  = new TObjNum(gsn);
      TObjNum *gmsnp = new TObjNum(gmsn);
      fPIDGen.Add(gpidp);  // *gpidp, *gsnp and *gmsnp stays
      fSNGen.Add(gsnp);    // but (TObjArray *)obj->SetOwner() deletes
      fMSNGen.Add(gmsnp);  // its elements.
    }
  }

  Int_t ncdctrk = ctp->GetCDCEntries();
  if ( ncdctrk > 0 ) {
    if (fDEBUG) cerr << "# of charged tracks = " << ncdctrk << endl;
    if ( (ctp->GetType()==2||ctp->GetType()==4) && ctp->GetCharge() == 0 ) {
      if (fDEBUG) cerr << "Charge of this track in mixed CAL cluster = "
		       << ctp->GetCharge() << endl
		       << "Don't use these CDCTrack pointer !" << endl;
      return;
    }
#if 0
    for (Int_t i = 0; i < ncdctrk; i++ ) {
#else // temporary treatment
    for (Int_t i = 0; i < 1; i++ ) {
#endif
      if (fDEBUG) cerr << "i = " << i << endl;
      Int_t gpid = 0;
      Int_t gsn  = 0;
      Int_t gmsn = ctp->GetCDCTrackAt(i)->GetGenID();
      while ( gmsn >= 0 ) {
        JSFGeneratorParticle *g = (JSFGeneratorParticle *)fGen->UncheckedAt(gmsn-1);
        gpid = g->GetID();
        gsn  = g->GetSerial();
        gmsn = g->GetMother();

        if (fDEBUG) cerr << "(PID, S.N, M.S.N) = ("
                         << gpid << ","
                         << gsn  << ","
                         << gmsn << ")" << endl;
      }
      if (fDEBUG) cerr << "-- Search ended --" << endl;
      TObjNum *gpidp = new TObjNum(gpid);
      TObjNum *gsnp  = new TObjNum(gsn);
      TObjNum *gmsnp = new TObjNum(gmsn);
      fPIDGen.Add(gpidp);  // *gpidp, *gsnp and *gmsnp stays
      fSNGen.Add(gsnp);    // but (TObjArray *)obj->SetOwner() deletes
      fMSNGen.Add(gmsnp);  // its elements.
    }
  }
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
Int_t TTL4JFlavourGetter::operator()(const ANLJet &jet) {

  static const Double_t kThetaCut = 10.0;
  static const Double_t kCHadRatioCut = 0.15;

  Int_t ntrkinjet = 0;
  Int_t ntrkfromb = 0;
  Int_t ntrkfromw = 0;
  Int_t ntrkfromchad = 0;

  SetData(jet);

  TIter nextptnid(&GetPartonID());
  TObjNum *ptnidp;
  while ((ptnidp = (TObjNum *)nextptnid())) {
    // fMSNGen : Not Yet contained primary hadron's Mother S.N contributing to
    //           the HDC cluster from neutral combined tracks.
    Int_t partonid = TMath::Abs(ptnidp->GetNum());
    if (fDEBUG) cerr << "TTL4JFlavourGetter::operator() : MSN = " << partonid << endl;

    if ( partonid == 3 ) ntrkfromb++;
    else if ( partonid == 7 || partonid == 9 ) ntrkfromw++;
  }

  TIter nexthadpid(&GetPrimaryHadronPID());
  TObjNum *pidp;
  while ((pidp = (TObjNum *)nexthadpid())) {
    // fPIDGen : Not Yet contained primary hadron's PIDs contributing to
    //           the HDC cluster from neutral combined tracks.
    Int_t hadpid = TMath::Abs(pidp->GetNum());
    if (fDEBUG) cerr << "TTL4JFlavourGetter::operator() : PID = " << hadpid << endl;
    ntrkinjet++;

    Int_t flavour = 0;
    if ( hadpid/1000 > 0 && hadpid/10000 == 0 ) { // Meson-Baryon branch
      flavour = hadpid/1000;
    } else {
      flavour = (hadpid/100)%10;
    }

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
	else if ( thetaspdnbarj < kThetaCut ) return -1;
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
	else if ( thetaspdnj < kThetaCut ) return 1;
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
