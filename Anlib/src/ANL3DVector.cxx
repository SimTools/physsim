//*************************************************************************
//* ===================
//*  ANL3DVector Class
//* ===================
//*
//* (Description)
//*    A very primitive lockable 3D vector class.
//* (Requires)
//*	class TVector3
//* 	class Lockable
//* 	class ANL2DVector
//* (Provides)
//* 	class ANL3DVector
//* (Update Recored)
//*    1999/09/13  K.Ikematsu	Original version.
//*    2000/03/23  K.Ikematsu	Added GetNorm method.
//*    2000/03/23  K.Ikematsu	Added GetTheta method.
//*    2000/03/28  K.Ikematsu	Added Acol method.
//*    2001/02/16  K.Ikematsu	Added GetTrans and GetLong method.
//*
//* $Id$
//*************************************************************************
//
#include "ANL3DVector.h"
//_____________________________________________________________________
//  -----------------
//  Lockable 3DVector
//  -----------------
//

ClassImp(ANL3DVector)
