//*************************************************************************
//* =============================
//*  ANLCheatedJetFinder Classes
//* =============================
//*
#include "ANLCheatedJetFinder.h"
//_____________________________________________________________________
//  ------------------
//  ANLTaggedJet Class
//  ------------------
//
ClassImp(ANLTaggedJet)

//_____________________________________________________________________
//  -------------------------
//  ANLCheatedJetFinder Class
//  -------------------------
//
Double_t ANLCheatedJetFinder::GetYmass(const ANL4DVector &p1, const ANL4DVector &p2) const {
  ANLTaggedJet *j1 = (ANLTaggedJet *)&p1;
  ANLTaggedJet *j2 = (ANLTaggedJet *)&p2;

  if ( j1->fTag != j2->fTag ) {
    return fEvis;
  } else {
    Double_t minE = TMath::Min(p1.E(),p2.E());
    return 2 * minE * minE * ( 1 - p1.CosTheta(p2) );
  }
}

void ANLCheatedJetFinder::Initialize(const TObjArray &parts) {
  fDone = kFALSE;
  DeleteJets();
  fEvis = 0.;
  //
  // Store each unlocked object as a single jet.
  //   
  TIter next(&parts);
  TObject *obj;
  while ((obj = next())) {
    if (!obj->IsA()->InheritsFrom("ANL4DVector")) {
      cerr << "ANLJetFinder::Initialize input is not a ANL4DVector" << endl;
      continue;
    }
    if (((ANL4DVector *)obj)->IsLocked()) continue;
    ANLTaggedJet *jet = new ANLTaggedJet();
    jet->Merge(obj);          // obj can be a jet, instead of being a track.
    fEvis += jet->E();
    fJets.Add(jet);
  }
}

void ANLCheatedJetFinder::CopyJets(const TObjArray &jets) {
  DeleteJets();
  TIter next(&jets);
  TObject *obj;
  while ((obj = next())) {
    ANLTaggedJet *jet = new ANLTaggedJet();
    jet->Merge((ANLJet *)obj);
    fJets.Add(jet);
  }
}

ClassImp(ANLCheatedJetFinder) 
