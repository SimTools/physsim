#ifndef EVENTSHAPE
#define EVENTSHAPE
//*************************************************************************
//* =======================
//*  ANLEventShape Classes 
//* =======================
//*
//* (Description)
//*    Jet finder classes for JLC analyses.
//* (Requires)
//* 	class TLorentzVector
//* 	class Lockable
//* 	class ANL4DVector
//* (Provides)
//* 	class ANLEventShape
//* (Usage)
//* (Update Recored)
//*    1999/08/11  K.Fujii	Original version derived from the LCD 
//*				version converted from java routines 
//*				written by Gary Bower.
//*				No essential changes except for the
//*				input event format.
//*    1999/09/05  K.Ikematsu	Replaced LockableLVector with ANL4DVector.
//*    1999/09/28  K.Fujii	Plugged memory leaks as recommended by
//*				M.Iwasaki.
//*
//*************************************************************************
//
#include <iostream.h>
#include "TClass.h"
#include "TMath.h"
#include "TMatrix.h"
#include "TRandom.h"
#include "TObjArray.h"
#include "ANL4DVector.h"
//_____________________________________________________________________
//  -------------------
//  ANLEventShape Class
//  -------------------
//

class ANLEventShape : public TObject{
public:
	ANLEventShape();
	~ANLEventShape();
	void setEvent(TObjArray* e);
   	void Initialize(const TObjArray& parts);
	void setThMomPower(Double_t tp);
	void setFast(Int_t nf);

	Double_t getThMomPower();
	Int_t    getFast();
	Double_t GetThrust() const;
	
	TVector3* thrustAxis();
	TVector3* majorAxis();
	TVector3* minorAxis();
	
	TVector3* thrust();
	Double_t  oblateness();

private:
	Double_t ulAngle(Double_t x, Double_t y);
	Double_t sign(Double_t a, Double_t b);

	void ludbrb(TMatrix *mom, 
				Double_t the, 
				Double_t phi, 
				Double_t bx, 
				Double_t by,
				Double_t bz);

	Int_t iPow(Int_t man, Int_t exp);

private:
	// PARU(41): Power of momentum dependence in sphericity finder.
	Double_t m_dSphMomPower; 

	// PARU(42): Power of momentum dependence in thrust finder.
	Double_t m_dDeltaThPower;
	
	// MSTU(44): # of initial fastest particles choosen to start search.
	Int_t m_iFast; 

	// PARU(48): Convergence criteria for axis maximization.
	Double_t m_dConv;

	// MSTU(45): # different starting configurations that must
	// converge before axis is accepted as correct.	
	Int_t m_iGood;

	// data: results
	// m_dAxes[1] is the Thrust axis.
	// m_dAxes[2] is the Major axis.
	// m_dAxes[3] is the Minor axis.
	TMatrix m_dAxes;

	TRandom m_random;

	Double_t m_dThrust[4];
	Double_t m_dOblateness;
	TVector3 m_EigenVector1;
	TVector3 m_EigenVector2;
	TVector3 m_EigenVector3;
	Double_t m_dEigenValue1;
	Double_t m_dEigenValue2;
	Double_t m_dEigenValue3;
	
	static Int_t m_maxpart;

	TVector3* m_ThrustAxis;
	TVector3* m_MajorAxis;
	TVector3* m_MinorAxis;
	TVector3* m_Thrust;

	ClassDef(ANLEventShape,0)	// ANLEventShape class
};

#endif
