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
#include "ANLPairCombiner.h"
//_____________________________________________________________________
//  -------------
//  ANLPair Class
//  -------------
//
//  In this implementaion the oder does not matter.
//  stack1 and stack2 are Stack's of lockable.
// 
//

ClassImp(ANLPair)

//_____________________________________________________________________
//  ---------------------
//  ANLPairCombiner Class
//  ---------------------
//
//  In this implementaion the oder does not matter.
//  stack1 and stack2 are Stack's of lockable.
// 
//
void ANLPairCombiner::Initialize(TObjArray &stack1, TObjArray &stack2)
{
   if (&stack1 == &stack2) {
      for (Int_t i = 0; i < (stack1.Capacity()-1); i++) {
         if (stack1[i] == 0 || 
             ((ANL4DVector *)stack1[i])->IsLocked()) continue;
         for (Int_t j = i + 1; j < stack2.Capacity(); j++) {
   	    if (stack2[j] == 0 ||
   	        ((ANL4DVector *)stack2[j])->IsLocked()) continue;
   	    fPairs.Add( new ANLPair(stack1[i],stack2[j]) );
   	 }
      }
   } else {
      for (Int_t i = 0; i < (stack1.Capacity()); i++) {
         if (stack1[i] == 0 || 
             ((ANL4DVector *)stack1[i])->IsLocked()) continue;
         for (Int_t j = 0; j < stack2.Capacity(); j++) {
   	    if (stack2[j] == 0 || 
   	        ((ANL4DVector *)stack2[j])->IsLocked()) continue;
   	    fPairs.Add( new ANLPair(stack1[i],stack2[j]) );
   	 }
      }
   }
}

void ANLPairCombiner::CleanUp()
{
   fPair.Reset(); 
   ANLPair *p; 
   while ( (p = (ANLPair *)fPair.Next()) ) { fPairs.Remove(p); }
}

ClassImp(ANLPairCombiner)

