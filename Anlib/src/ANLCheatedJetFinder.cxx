//*************************************************************************
//* =============================
//*  ANLCheatedJetFinder Classes
//* =============================
//*
//* $Id$
//*
#include "ANLCheatedJetFinder.h"
//_____________________________________________________________________
//  ------------------
//  ANLTaggedJet Class
//  ------------------
//
ClassImp(ANLTaggedJet)

//*--
//*  Setters
//*--
void ANLTaggedJet::Merge(TObject *part) {
  Add(part);
  fTag = ((ANLTrack *)part)->GetColorSingletID();
}

void ANLTaggedJet::Merge(ANLJet *jet) {
#ifdef __DEBUG__
  cerr << "ANLTaggedJet::Merge() is called ..." << endl;
  cerr << "   fTag = " << fTag << endl;
  cerr << "ANLTaggedJet::Merge ; ANLJet::Merge(jet); is called ..." << endl;
#endif
  ANLJet::Merge(jet);
  if ( fTag == 9999 ) fTag = ((ANLTaggedJet *)jet)->GetTag();
}

//_____________________________________________________________________
//  -------------------------
//  ANLCheatedJetFinder Class
//  -------------------------
//
ClassImp(ANLCheatedJetFinder)

//*--
//*  Basic Services
//*--
Double_t ANLCheatedJetFinder::GetYmass(const ANL4DVector &p1, const ANL4DVector &p2) const {
  ANLTaggedJet *j1 = (ANLTaggedJet *)&p1;
  ANLTaggedJet *j2 = (ANLTaggedJet *)&p2;

  Int_t j1tag = j1->GetTag();
  Int_t j2tag = j2->GetTag();

  if ( j1tag != 9999 && j2tag != 9999 && j1tag != j2tag ) {
#ifdef __DEBUG__
    cerr << "ANLCheatedJetFinder::GetYmass ; (Ymass,j1tag,j2tag) = ("
	 << fEvis*fEvis << "," << j1tag << "," << j2tag << ")" << endl;
#endif
    return fEvis*fEvis;
  } else {
    Double_t minE = TMath::Min(p1.E(),p2.E());
#ifdef __DEBUG__
    cerr << "ANLCheatedJetFinder::GetYmass ; (Ymass,j1tag,j2tag) = ("
	 << 2 * minE * minE * ( 1 - p1.CosTheta(p2) )
	 << "," << j1tag << "," << j2tag << ")" << endl;
#endif
    return 2 * minE * minE * ( 1 - p1.CosTheta(p2) );
  }
}
