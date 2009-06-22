#ifndef __ZZZSpring__
#define __ZZZSpring__

//////////////////////////////////////////////////////////////////////////
//
// ZZZSpring
//
// e+e- -> ZZZ
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
//  COMMONs for ZZZ calculation, see USRPRM.inc for more details.
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
  Int_t    imd3lo;	// 1st Z0 decay mode ID.
  Int_t    imd3hi;	// lst Z0 decay mode ID.
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class ZZZBases  : public JSFBases {
protected:
  Double_t fRoots;	// CM Energy - 2m_t(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;	// beam energy spread
  Int_t    fISRBM;	// Flag for ISR & Beamstrahlung
  Int_t    fWmModesLo;	// 1st Z1 decay mode ID.
  Int_t    fWmModesHi;	// 1st Z1 decay mode ID.
  Int_t    fWpModesLo;	// 2nd Z2 decay mode ID.
  Int_t    fWpModesHi;	// 2nd Z2 decay mode ID.
  Int_t    fZModesLo;	// 3rd Z0 decay mode ID.
  Int_t    fZModesHi;	// 3rd Z0 decay mode ID.
  Double_t fAlphai;	// 1/alpha(m_Z)
  Double_t fAlphas;	// alpha_s(m_Z)
  Double_t fMassW;	// m_W
  Double_t fMassZ;	// m_Z
  Double_t fMassHiggs;	// m_H
  Double_t fMassTop;	// m_t
  Int_t   fISHUFL[50];  // random number shuffler
public:
  ZZZBases(const char *name="ZZZBases", 
	     const char *title="ZZZbar  Bases");
  virtual ~ZZZBases(){} 

  void SetRoots(Double_t deltrs){ fRoots=deltrs; }
  Double_t GetRoots(){ return fRoots;}
  
  void SetPolElectron(Double_t polebm){ fPolElectron=polebm; }
  Double_t GetPolElectron(){ return fPolElectron;}
  
  void SetSigmaEbeam(Double_t sgmebm){ fSigmaEbeam=sgmebm; }
  Double_t GetSigmaEbeam(){ return fSigmaEbeam;}
  
  void SetISRBM(Int_t isrbm){ fISRBM=isrbm; }
  Int_t GetISRBM(){ return fISRBM;}
   
  void SetWmModesLo(Int_t i){ fWmModesLo=i; }
  Int_t GetWmModesLo(){ return fWmModesLo;}
   
  void SetWmModesHi(Int_t i){ fWmModesHi=i; }
  Int_t GetWmModesHi(){ return fWmModesHi;}
   
  void SetWpModesLo(Int_t i){ fWpModesLo=i; }
  Int_t GetWpModesLo(){ return fWpModesLo;}
   
  void SetWpModesHi(Int_t i){ fWpModesHi=i; }
  Int_t GetWpModesHi(){ return fWpModesHi;}
   
  void SetZModesLo(Int_t i){ fZModesLo=i; }
  Int_t GetZModesLo(){ return fZModesLo;}
   
  void SetZModesHi(Int_t i){ fZModesHi=i; }
  Int_t GetZModesHi(){ return fZModesHi;}
   
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

  ClassDef(ZZZBases,1)  // Bases for e+e- -> ZZZ process

};

class ZZZSpring;

class ZZZSpringBuf : public JSFSpringBuf {
public:
//  <<+ZZZbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-ZZZbar>>
public:
  ZZZSpringBuf(const char *name="ZZZSpringBuf", 
	     const char *title="ZZZSpring test event buffer",
	     ZZZSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ZZZSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(ZZZSpringBuf,1)  // ZZZSpring event buffer
};

class ZZZSpring : public JSFSpring {
public:
   ZZZSpring(const char *name="ZZZSpring", 
	      const char *title="ZZZSpring test",
             ZZZBases *bases=NULL);
   virtual ~ZZZSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(ZZZSpring,1)  // ZZZSpring class
};
#endif
