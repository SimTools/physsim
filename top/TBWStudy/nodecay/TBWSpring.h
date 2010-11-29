#ifndef TBWSPRING_H
#define TBWSPRING_H

//////////////////////////////////////////////////////////////////////////
//
// TBWSpring
//
// e+e- -> tbwbar
//
//////////////////////////////////////////////////////////////////////////

#include "JSFConfig.h"

#include "TNamed.h"
#include "TMath.h"
#include "TDatime.h"

#include "JSFModule.h"
#include "JSFBases.h"
#include "JSFSpring.h"


// =====================================================================
//  COMMONs for shuffling random numbers.
// =====================================================================
typedef struct {
  Int_t nzz;
  Int_t ishufl[50];
} COMMON_BSHUFL;             //  Common /BSHUFL/

extern COMMON_BSHUFL bshufl_;

// =====================================================================
//  COMMONs for SM parameters, see USMPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  alfi, alfs, amsw, amsz, amsh, amst;
} COMMON_USMPRM;             //  Common /USMPRM/

extern COMMON_USMPRM usmprm_;

// =====================================================================
//  COMMONs for TBW calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts, polebm; // Beam parameters.
  Int_t    isrb;	// ISR and BM flag.
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class TBWBases  : public JSFBases {
protected:
  Double_t fRoots;	// CM Energy(GeV)
  Double_t fPolElectron;// e- polarization
  Int_t    fISRBM;	// Flag for ISR & Beamstrahlung
  Double_t fAlphai;	// 1/alpha(m_Z)
  Double_t fAlphas;	// alpha_s(m_Z)
  Double_t fMassW;	// m_W
  Double_t fMassZ;	// m_Z
  Double_t fMassTop;	// m_t
  Double_t fMassHiggs;	// m_h
  Int_t   fISHUFL[50];  // random number shuffler

public:
  TBWBases(const char *name="TBWBases", 
	     const char *title="TBWbar  Bases");
  virtual ~TBWBases(){} 

  void SetRoots(Double_t rs){ fRoots=rs; }
  Double_t GetRoots(){ return fRoots;}
  
  void SetPolElectron(Double_t polebm){ fPolElectron=polebm; }
  Double_t GetPolElectron(){ return fPolElectron;}
  
  void SetISRBM(Int_t isrbm){ fISRBM=isrbm; }
  Int_t GetISRBM(){ return fISRBM;}
   
  void SetAlphai(Double_t alphai){ fAlphai=alphai; }
  Double_t GetAlphai(){ return fAlphai;}
   
  void SetAlphas(Double_t alphas){ fAlphas=alphas; }
  Double_t GetAlphas(){ return fAlphas;}
   
  void SetMassW(Double_t mass){ fMassW=mass; }
  Double_t GetMassW(){ return fMassW;}
   
  void SetMassZ(Double_t mass){ fMassZ=mass; }
  Double_t GetMassZ(){ return fMassZ;}
   
  void SetMassHiggs(Double_t mass){ fMassHiggs=mass; }
  Double_t GetMassHiggs(){ return fMassHiggs;}
   
  void SetMassTop(Double_t mass){ fMassTop=mass; }
  Double_t GetMassTop(){ return fMassTop;}

  void Userin();   // Bases user initialization
  void Userout();  // Bases user output 
  Double_t Func(); // New style not yet
  Double_t Func(Double_t x[]); // Bases integration function.
  void PrintParameters(); // Print parameters

  ClassDef(TBWBases,1)  // Bases for e+e- -> tbwbar process
};

class TBWSpring;

class TBWSpringBuf : public JSFSpringBuf {
public:
//  <<+TBWbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-TBWbar>>
public:
  TBWSpringBuf(const char *name="TBWSpringBuf", 
	     const char *title="TBWSpring test event buffer",
	     TBWSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~TBWSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(TBWSpringBuf,1)  // TBWSpring event buffer
};

class TBWSpring : public JSFSpring {
public:
   TBWSpring(const char *name="TBWSpring", 
	      const char *title="TBWSpring test",
             TBWBases *bases=NULL);
   virtual ~TBWSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(TBWSpring,1)  // TBWSpring class
};
#endif
