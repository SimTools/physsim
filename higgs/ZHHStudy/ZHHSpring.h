#ifndef ZZSPRING_H
#define ZZSPRING_H

//////////////////////////////////////////////////////////////////////////
//
// ZHHSpring
//
// e+e- -> ZHH
//
//////////////////////////////////////////////////////////////////////////

#include "JSFConfig.h"

#include "TNamed.h"
#include "TMath.h"
#include "TDatime.h"

#ifndef __JSFModule__
#include "JSFModule.h"
#endif
#ifndef __JSFBases__
#include "JSFBases.h"
#endif
#ifndef __JSFSpring__
#include "JSFSpring.h"
#endif


// =====================================================================
//  COMMONs for shuffling random numbers.
// =====================================================================
typedef struct {
  Int_t nZHH;
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
//  COMMONs for ZHH calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts;	// sqrt(s).
  Float_t  polebm;	// electron beam polarization.
  Float_t  sgmebm;	// beam energy spread (fraction).
  Int_t    isrb;	// ISR and BM flag.
  Int_t    imd1lo;	// 1st Z0 decay mode ID.
  Int_t    imd1hi;	// lst Z0 decay mode ID.
  Int_t    imd2lo;	// 1st H0 decay mode ID (ffbar only).
  Int_t    imd2hi;	// lst H0 decay mode ID (ffbar only).
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class ZHHBases  : public JSFBases {
protected:
  Double_t fRoots;	// CM Energy - 2m_t(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;	// beam energy spread
  Int_t    fISRBM;	// Flag for ISR & Beamstrahlung
  Int_t    fZ0ModesLo;	// 1st Z0 decay mode ID.
  Int_t    fZ0ModesHi;	// lst Z0 decay mode ID.
  Double_t fAlphai;	// 1/alpha(m_Z)
  Double_t fAlphas;	// alpha_s(m_Z)
  Double_t fMassW;	// m_W
  Double_t fMassZ;	// m_Z
  Double_t fMassHiggs;	// m_H
  Double_t fMassTop;	// m_t
  Int_t   fISHUFL[50];  // random number shuffler
public:
  ZHHBases(const char *name="ZHHBases", 
	     const char *title="ZHHbar  Bases");
  virtual ~ZHHBases(){} 

  void SetRoots(Double_t deltrs){ fRoots=deltrs; }
  Double_t GetRoots(){ return fRoots;}
  
  void SetPolElectron(Double_t polebm){ fPolElectron=polebm; }
  Double_t GetPolElectron(){ return fPolElectron;}
  
  void SetSigmaEbeam(Double_t sgmebm){ fSigmaEbeam=sgmebm; }
  Double_t GetSigmaEbeam(){ return fSigmaEbeam;}
  
  void SetISRBM(Int_t isrbm){ fISRBM=isrbm; }
  Int_t GetISRBM(){ return fISRBM;}
   
  void SetZ0ModesLo(Int_t i){ fZ0ModesLo=i; }
  Int_t GetZ0ModesLo(){ return fZ0ModesLo;}
  
  void SetZ0ModesHi(Int_t i){ fZ0ModesHi=i; }
  Int_t GetZ0ModesHi(){ return fZ0ModesHi;}
   
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
  Double_t Func(Double_t x[]); // Bases integration function.
  Double_t Func(); // Bases integration function.
  void PrintParameters(); // Print parameters

  ClassDef(ZHHBases,1)  // Bases for e+e- -> ZHH process

};

class ZHHSpring;

class ZHHSpringBuf : public JSFSpringBuf {
public:
//  <<+ZHHbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-ZHHbar>>
public:
  ZHHSpringBuf(const char *name="ZHHSpringBuf", 
	     const char *title="ZHHSpring test event buffer",
	     ZHHSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ZHHSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(ZHHSpringBuf,1)  // ZHHSpring event buffer
};

class ZHHSpring : public JSFSpring {
public:
   ZHHSpring(const char *name="ZHHSpring", 
	      const char *title="ZHHSpring test",
             ZHHBases *bases=NULL);
   virtual ~ZHHSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(ZHHSpring,1)  // ZHHSpring class
};
#endif
