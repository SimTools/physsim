#ifndef __EEZZSpring__
#define __EEZZSpring__

//////////////////////////////////////////////////////////////////////////
//
// EEZZSpring
//
// e+ e- -> e+ e- Z Z
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
//  COMMONs for EEZZ calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts;	// sqrt(s).
  Float_t  polebm;	// electron beam polarization.
  Float_t  sgmebm;	// beam energy spread (fraction).
  Int_t    isrb;	// ISR and BM flag.
  Int_t    imd1lo;	// 1st Z1 decay mode ID.
  Int_t    imd1hi;	// lst Z1 decay mode ID.
  Int_t    imd2lo;	// 1st Z2 decay mode ID.
  Int_t    imd2hi;	// lst Z2 decay mode ID.
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class EEZZBases  : public JSFBases {
protected:
  Double_t fRoots;	// CM Energy - 2m_t(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;	// beam energy spread
  Int_t    fISRBM;	// Flag for ISR & Beamstrahlung
  Int_t    fZ1ModesLo;	// 1st Z1 decay mode ID.
  Int_t    fZ1ModesHi;	// lst Z1 decay mode ID.
  Int_t    fZ2ModesLo;	// 1st Z2 decay mode ID.
  Int_t    fZ2ModesHi;	// lst Z2 decay mode ID.
  Double_t fAlphai;	// 1/alpha(m_Z)
  Double_t fAlphas;	// alpha_s(m_Z)
  Double_t fMassW;	// m_W
  Double_t fMassZ;	// m_Z
  Double_t fMassHiggs;	// m_H
  Double_t fMassTop;	// m_t
  Int_t   fISHUFL[50];  // random number shuffler
public:
  EEZZBases(const char *name="EEZZBases", 
	     const char *title="EEZZbar  Bases");
  virtual ~EEZZBases(){} 

  void SetRoots(Double_t deltrs){ fRoots=deltrs; }
  Double_t GetRoots(){ return fRoots;}
  
  void SetPolElectron(Double_t polebm){ fPolElectron=polebm; }
  Double_t GetPolElectron(){ return fPolElectron;}
  
  void SetSigmaEbeam(Double_t sgmebm){ fSigmaEbeam=sgmebm; }
  Double_t GetSigmaEbeam(){ return fSigmaEbeam;}
  
  void SetISRBM(Int_t isrbm){ fISRBM=isrbm; }
  Int_t GetISRBM(){ return fISRBM;}
   
  void SetZ1ModesLo(Int_t i){ fZ1ModesLo=i; }
  Int_t GetZ1ModesLo(){ return fZ1ModesLo;}
   
  void SetZ1ModesHi(Int_t i){ fZ1ModesHi=i; }
  Int_t GetZ1ModesHi(){ return fZ1ModesHi;}
   
  void SetZ2ModesLo(Int_t i){ fZ2ModesLo=i; }
  Int_t GetZ2ModesLo(){ return fZ2ModesLo;}
   
  void SetZ2ModesHi(Int_t i){ fZ2ModesHi=i; }
  Int_t GetZ2ModesHi(){ return fZ2ModesHi;}
   
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

  ClassDef(EEZZBases,1)  // Bases for e+e- -> EEZZ process

};

class EEZZSpring;

class EEZZSpringBuf : public JSFSpringBuf {
public:
//  <<+EEZZbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-EEZZbar>>
public:
  EEZZSpringBuf(const char *name="EEZZSpringBuf", 
	     const char *title="EEZZSpring test event buffer",
	     EEZZSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~EEZZSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(EEZZSpringBuf,1)  // EEZZSpring event buffer
};

class EEZZSpring : public JSFSpring {
public:
   EEZZSpring(const char *name="EEZZSpring", 
	      const char *title="EEZZSpring test",
             EEZZBases *bases=NULL);
   virtual ~EEZZSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(EEZZSpring,1)  // EEZZSpring class
};


#endif









