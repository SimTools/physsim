#ifndef __ANLJETFINDER__
#define __ANLJETFINDER__
//*************************************************************************
//* ======================
//*  ANLJetFinder Classes
//* ======================
//*
//* (Description)
//*    Jet finder classes for JLC analyses
//* (Requires)
//* 	class TLorentzVector
//* 	class Lockable
//* 	class ANL4DVector
//* (Provides)
//* 	class ANLJet
//* 	class ANLJetFinder
//* 	class ANLJadeJetFinder
//* 	class ANLJadeEJetFinder
//* 	class ANLDurhamJetFinder
//* (Usage)
//*     // Example
//*     Double_t ycut = 0.01;
//*	ANLJadeEJetFinder jclust(ycut);
//*     jclust.Initialize(tracks); // tracks: TObjArray of LVector derivatives.
//*     jclust.FindJets();	   // finds jets with ycut = 0.01.
//*     ycut = 0.015;		   
//*	jclust.SetYcut(ycut);	   // One can make the ycut "bigger"
//*	jclust.FindJets();	   // and resume with the new ycut.
//*     ycut = 0.05;
//*	ANLJadeEJetFinder jclust2(jclust);
//*	jclust2.SetYcut(ycut);	   // One can copy the previous jet finder
//*	jclust2.FindJets();	   // and resume with the new "bigger" ycut.
//*	Int_t njets = 2;	   // One can also force the event to be
//* 	jclust2.ForceNJets(njets); // "njets" jets.
//* (Caution)
//*	One CANNOT decrease the ycut without restoring track information
//*	by reinvoking the "Initialize" member function.
//*	On the other hand, one can increase the ycut and resume jet
//*	finding starting from the last result, thereby saving CPU time.
//* (Update Recored)
//*    1999/07/30  K.Fujii	Original version derived from Rob Shanks's
//*				LCD version. This version heavily uses
//*				ROOT facilities to simplify algorithm.
//*				It allows one to "increase" the last ycut
//*				and resume jet finding. It also allows one
//*				to force an event to be "n" jets.
//*    1999/08/14  K.Fujii	Modified GetYmass as suggested by M.Iwasaki.
//*    1999/09/05  K.Ikematsu   Replaced LockableLVector with ANL4DVector.
//*    1999/09/14  K.Ikematsu   Added GetYmax method.
//*    2000/03/28  K.Ikematsu   Bug fixed about fYmassMax
//*
//*************************************************************************
//
#include <iostream.h>
#include "TClass.h"
#include "TMatrix.h"
#include "TObjArray.h"
#include "ANL4DVector.h"
//_____________________________________________________________________
//  ------------
//  ANLJet Class
//  ------------
//
//  For the moment, the ANLJet class only accepts, as its elements,
//  objects derived from ANL4DVector's or ANLJet's themselves, 
//  although it is more desirable that it accepts any TObjects
//  that contain ANL4DVector's. 
//

class ANLJet : public ANL4DVector {
friend class ANLJetFinder;
public:
   ANLJet();
   ANLJet(const TObjArray &parts);
   virtual ~ANLJet();

   Int_t	          GetNparticles() const;
   const TObjArray        &GetParticlesInJet() const;
   const ANL4DVector  &operator()() const;

   void	Add(TObject *part);
   void Merge(TObject *part);
   void Merge(ANLJet *jet);
   void	Remove(TObject *part);

   void DebugPrint(const Char_t *opt = "Brief") const;

private:
   TObjArray fParts;

   ClassDef(ANLJet,1)  // ANLJet class
};

//_____________________________________________________________________
//  ------------------
//  ANLJetFinder Class
//  ------------------
//
//  For the moment, the ANLJet class only accepts, as its elements,
//  objects derived from ANL4DVector's or ANLJet's themselves, 
//  although it is more desirable that it accepts any TObjects
//  that contain ANL4DVector's. 
//
class ANLJetFinder {
public:
   ANLJetFinder(Double_t y = 0.);		// default constructor
   ANLJetFinder(const ANLJetFinder &jf);	// copy constructor
   virtual ~ANLJetFinder();			// default destructor
   
   Bool_t     IsInitialized() const;		// inquires if fJets set
   Double_t   GetYcut() const;			// returns ycut
   Double_t   GetYmax();                        // returns ymax
   Int_t      GetNjets();			// returns No. of jets   
   TObjArray &GetJets();			// returns array of jets

   void SetYcut(Double_t ycut);			// sets ycut
   void Initialize(const TObjArray &parts);	// sets each particle as a jet
   void FindJets();				// finds jets
   void ForceNJets(Int_t njets);		// forces the event to njets

   ANLJetFinder & operator=(const ANLJetFinder & jf); 	// copys a jet finder

   virtual Double_t GetYmass(const ANL4DVector &p1, 
   			     const ANL4DVector &p2) const; 
private:
   void CopyJets(const TObjArray &jets);
   void DeleteJets();

private:
   Bool_t	   fDone;	//  (1,0) = (jet finding done, not yet)
   Double_t 	   fYcut;	//  ycut
   Double_t	   fEvis;	//  Evis
   TObjArray	   fJets;	//  Jets in the event
   TMatrix        *fYmass;	//! pointer to pair mass array
   Double_t        fYmassMax;   //  Ymass_max

   ClassDef(ANLJetFinder,1)  	// ANLJetFinder class
};


//_____________________________________________________________________
//  ----------------------
//  ANLJadeJetFinder Class
//  ----------------------
//
class ANLJadeJetFinder : public ANLJetFinder {
public:
   ANLJadeJetFinder(Double_t y = 0.);

   Double_t GetYmass(const ANL4DVector &p1, const ANL4DVector &p2) const;
   
   ClassDef(ANLJadeJetFinder,1) 	// Jade Jet Finder class
};

//_____________________________________________________________________
//  -----------------------
//  ANLJadeEJetFinder Class
//  -----------------------
//
class ANLJadeEJetFinder : public ANLJetFinder {
public:
   ANLJadeEJetFinder(Double_t y = 0.);

   Double_t GetYmass(const ANL4DVector &p1, const ANL4DVector &p2) const;
   
   ClassDef(ANLJadeEJetFinder,1) 	// Jade EJet Finder class
};

//_____________________________________________________________________
//  ------------------------
//  ANLDurhamJetFinder Class
//  ------------------------
//
class ANLDurhamJetFinder : public ANLJetFinder {
public:
   ANLDurhamJetFinder(Double_t y = 0.);

   Double_t GetYmass(const ANL4DVector &p1, const ANL4DVector &p2) const;
   
   ClassDef(ANLDurhamJetFinder,1) 	// Durham Jet Finder class
};

#endif
