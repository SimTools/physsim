#ifndef __ANLPAIRCOMBINER__
#define __ANLPAIRCOMBINER__
//*************************************************************************
//* ======================
//*  PairCombiner Classes
//* ======================
//*
//* (Description)
//*    A very primitive class library for pair combiners.
//* (Requires)
//* 	class Lockable
//*	class TLorentzVector
//*	class ANL4DVector
//* (Provides)
//* 	class ANLPair
//*	class ANLPairCombiner
//* (Update Recored)
//*    1999/06/05  K.Fujii	Original very primitive version.
//*    1999/08/08  K.Fujii	new members to handle quality.
//*    1999/09/05  K.Ikematsu   Replaced LockableLVector with ANL4DVector.
//*
//*************************************************************************
//
#include "TClass.h"
#include "TObjArray.h"
#include <iostream.h>
#include "ANL4DVector.h"
//_____________________________________________________________________
//  -------------
//  ANLPair Class
//  -------------
//
//  In this implementaion the oder does not matter.
//  The elements are TObject's that should inherit from ANL4DVector
//  for the moment: they can actually be made into any "Lockable" that 
//  supports an operator "+" ("LockableAddable"). Whether we need to
//  prepare a "LockableAddable" class as a (pure) virtual base class
//  is not clear at the moment. 
// 
//
class ANLPair : public ANL4DVector {
friend class ANLPairCombiner;
public:
   ANLPair() : ANL4DVector(0.), fQuality(0.) { fP[0] = 0; fP[1] = 0; }
   ANLPair(TObject *q1p, TObject *q2p, Double_t quality = 0.) 
   	: ANL4DVector(*(ANL4DVector *)q1p + *(ANL4DVector *)q2p),
   	  fQuality(quality)
   {
   	fP[0] = q1p;
   	fP[1] = q2p;
   }
   ANLPair(const ANLPair &p) 
   	: ANL4DVector(p), fQuality(p.fQuality)
   {
   	fP[0] = p.fP[0];
   	fP[1] = p.fP[1];
   }
   virtual ~ANLPair() {}
   
   TObject         *operator[](Int_t i) const { return  fP[i]; }
   ANL4DVector  operator()(Int_t i) const { return *(ANL4DVector *)fP[i]; }
   ANL4DVector  operator()()        const { return *(ANL4DVector *)this;  }

   void  LockChildren()
   {
   	((ANL4DVector *)fP[0])->Lock();
   	((ANL4DVector *)fP[1])->Lock();
   }
   void  UnlockChildren()
   {
   	((ANL4DVector *)fP[0])->Unlock();
   	((ANL4DVector *)fP[1])->Unlock();
   }
   Bool_t IsLocked() const
   {
   	return ((ANL4DVector *)fP[0])->IsLocked() ||
   	       ((ANL4DVector *)fP[1])->IsLocked() ||
   	       ANL4DVector::IsLocked();
   }

   inline void SetQuality(Double_t q) { fQuality = q; }
   inline Double_t GetQuality() const { return fQuality; }
   
   Bool_t IsSortable() const { return kTRUE; }
   Int_t  Compare(const TObject *obj) const
   {
   	if (fQuality < ((ANLPair *)obj)->fQuality) return -1;
   	else if (fQuality > ((ANLPair *)obj)->fQuality) return 1;
   	else return 0;
   }

   virtual void Delete(Option_t *opt="") { delete fP[0]; delete fP[1]; }
   
   void DebugPrint(const Char_t *opt = "Brief") const
   {
   	if (opt == "Brief") {
           cerr << (void *)this << ": ";
 	   cerr << " item 0: " << (void *)fP[0]
  	        << (((ANL4DVector *)fP[0])->IsLocked() ? ":locked" : "       ")
  	        << " item 1: " << (void *)fP[1] 
  	        << (((ANL4DVector *)fP[1])->IsLocked() ? ":locked" : "       ")
  	        << " mass = " << GetMass() << endl;
   	} else {
           cerr << (void *)this << ": " <<  endl;
  	   cerr << "item 0: " << (void *)fP[0] 
  	        << (((ANL4DVector *)fP[0])->IsLocked() ? ":locked" : "       ")
  	        << ":" << endl; 
   	      ((ANL4DVector *)fP[0])->DebugPrint(opt); 
  	   cerr << "item 1: " << (void *)fP[1]
  	        << (((ANL4DVector *)fP[1])->IsLocked() ? ":locked" : "       ")
  	        << ":" << endl; 
   	      ((ANL4DVector *)fP[1])->DebugPrint(opt);
   	} 
   }
private:
   TObject *fP[2];	// pointers to 4-vectors in the pair
   Double_t fQuality;	// quality

   ClassDef(ANLPair,1)  // Lorentz vector pair class
};


//_____________________________________________________________________
//  ---------------------
//  ANLPairCombiner Class
//  ---------------------
//
//  In this implementaion the oder does not matter.
//  stack1 and stack2 are Stack's of lockable.
// 
//
class ANLPairCombiner {
public:
   ANLPairCombiner() : fToDelete(kFALSE), fPairs(1), fPair(&fPairs) 
   	{ fStacks[0] = 0; fStacks[1] = 0; }
   ANLPairCombiner(ANLPairCombiner & pair) 
   	: fToDelete(kFALSE), fPairs(pair.fPairs), fPair(pair.fPair)
   {
   	fStacks[0] = pair.fStacks[0]; fStacks[1] = pair.fStacks[1];
   }
   ANLPairCombiner(TObjArray &stack1, TObjArray &stack2) 
   	: fToDelete(kTRUE), fPair(&fPairs)
   {
   	fStacks[0] = &stack1; fStacks[1] = &stack2;
   	Initialize(*fStacks[0],*fStacks[1]);
   }
   virtual ~ANLPairCombiner()
   {
   	if (fToDelete) fPairs.Delete();
   }

   void      CleanUp();
   TObject * operator()()  { return fPair.Next(); }
   TObject * Next ()       { return fPair.Next(); }
   void      Reset()       { fPair.Reset();  }
   
   void DebugPrint(const Char_t *opt = "Brief") const
   {
   	TObjArrayIter next(&fPairs);
   	ANLPair *p;
   	while ((p = (ANLPair *)next())) {
   		p->DebugPrint(opt);
   	}
   }

private:
   void Initialize(TObjArray &stack1, TObjArray &stack2);

private:
   Bool_t        fToDelete;	// status flag
   TObjArray    *fStacks[2];	// pointers to 4-vector stacks to combine
   TObjArray     fPairs;	// array of ANLPair's
   TObjArrayIter fPair;		// array Iterator

   ClassDef(ANLPairCombiner,1)  // Pair combiner class
};

#endif
