#ifndef __ANLCHEATEDJETFINDER__
#define __ANLCHEATEDJETFINDER__
//*************************************************************************
//* $Id$
//*
#include "ANLTrack.h"
#include "ANLJetFinder.h"
//_____________________________________________________________________
//  ------------------
//  ANLTaggedJet Class
//  ------------------
//
class ANLTaggedJet : public ANLJet {
friend class ANLCheatedJetFinder;
public:
   ANLTaggedJet() : fTag(9999) {}
   ANLTaggedJet(const TObjArray &parts) : ANLJet(parts), fTag(9999) {}
   virtual ~ANLTaggedJet() {}

   void Merge(TObject *part);
   void Merge(ANLJet *jet);
   void SetTag(Int_t tag) { fTag = tag; }
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
   ANLCheatedJetFinder(const ANLCheatedJetFinder &jf) : ANLJetFinder(jf) {
#ifdef __DEBUG__
     cerr << "  copy constructor is called ..." << endl;
#endif
     ANLTaggedJet *jetsrc;
     ANLTaggedJet *jet;
     TIter nextsrc( & (((ANLCheatedJetFinder &)jf).GetJets()));
     TIter next(&GetJets());
     while ((jet = (ANLTaggedJet *)next())) {
       jetsrc = (ANLTaggedJet *)nextsrc();
       jet->SetTag(jetsrc->GetTag());
#ifdef __DEBUG__
       cerr << "  fTag = " << jet->GetTag() << endl;
#endif
     }
   }

   ANLJet *NewJet() {               // make new ANLTaggedJet object
     return new ANLTaggedJet();
   }
   ANLJetFinder *NewJetFinder(ANLJetFinder *jf) {
#ifdef __DEBUG__
     cerr << "ANLCheatedJetFinder::NewJetFinder() is called ..." << endl;
#endif
     return new ANLCheatedJetFinder(*(ANLCheatedJetFinder *)jf);
                                    // make new ANLCheatedJetFinder object
   }



   Double_t GetYmass(const ANL4DVector &p1, const ANL4DVector &p2) const;

   ClassDef(ANLCheatedJetFinder,1)  // Cheated Jet Finder class
};

#endif
