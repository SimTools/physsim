#ifndef __LOCKABLE__
#define __LOCKABLE__
//*************************************************************************
//* ================
//*  Lockable Class
//* ================
//*
//* (Description)
//*    Lockable class supplies basic functions for lockable objects.
//* (Requires)
//* 	none
//* (Provides)
//* 	class Lockable
//* (Update Recored)
//*    1999/06/05  K.Fujii	Original very primitive version.
//*
//*************************************************************************
//
#include <Rtypes.h>
//_____________________________________________________________________
//  ------------------------------
//  Base Class for Lockale Objects
//  ------------------------------
//
class Lockable {
public:
   Lockable();
   virtual ~Lockable();
   virtual Bool_t IsLocked() const;
   virtual void   Lock();
   virtual void   Unlock();
private:
   Bool_t fStatus;	 // lock byte

   ClassDef(Lockable,1)  // Base class for lockable objects
};

#endif
