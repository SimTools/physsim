#ifndef __ANLCHEATEDJETFINDER__
#define __ANLCHEATEDJETFINDER__
//*************************************************************************
//* $Id$
//*
#include "ANLJetFinder.h"

//_____________________________________________________________________
//  ------------------
//  ANLTaggedJet Class
//  ------------------
//
class ANLTaggedJet : public ANLJet {
friend class ANLCheatedJetFinder;
public:
   ANLTaggedJet() {} 
   ANLTaggedJet(const TObjArray &parts) : ANLJet(parts), fTag(-1) {}
   virtual ~ANLTaggedJet() {}

   inline Int_t GetTag() const { return fTag; }
private:
   Int_t fTag;  // Color singlet ID

   ClassDef(ANLTaggedJet,1)  // ANLTaggedJet class
};

//_____________________________________________________________________
//  -------------------------
//  ANLCheatedJetFinder Class
//  -------------------------
//
class ANLCheatedJetFinder : public ANLJetFinder {
public:
   ANLCheatedJetFinder(Double_t y = 0.) : ANLJetFinder(y) {}  // default constructor

   Double_t GetYmass(const ANL4DVector &p1, const ANL4DVector &p2) const;
   void Initialize(const TObjArray &parts);
   void CopyJets(const TObjArray &jets);

   ClassDef(ANLCheatedJetFinder,1)  // Cheated Jet Finder class
};

#endif
