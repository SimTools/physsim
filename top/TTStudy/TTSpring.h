#ifndef __TTSpring__
#define __TTSpring__

//////////////////////////////////////////////////////////////////////////
//
// TTSpring
//
// e+e- -> ttbar
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
//  COMMONs for TT calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  deltrs, polebm, sgmebm;
  Int_t    isrb;	// ISR and BM flag.
  Int_t    imd1lo;	// 1st W- decay mode ID.
  Int_t    imd1hi;	// lst W- decay mode ID.
  Int_t    imd2lo;	// 1st W+ decay mode ID.
  Int_t    imd2hi;	// lst W+ decay mode ID.
  Float_t  vkmt, beth;
  Int_t    nqcd;
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class TTBases  : public JSFBases {
protected:
  Double_t fDeltaRoots;	// CM Energy - 2m_t(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;	// beam energy spread
  Int_t    fISRBM;	// Flag for ISR & Beamstrahlung
  Int_t    fWmModesLo;	// 1st W- decay mode ID.
  Int_t    fWmModesHi;	// lst W- decay mode ID.
  Int_t    fWpModesLo;	// 1st W+ decay mode ID.
  Int_t    fWpModesHi;	// lst W+ decay mode ID.
  Double_t fVkmTB;	// Normalized V_km[tb]
  Double_t fBetaH;	// Normalized Top's Yukawa
  Int_t    fNRQCD;	// NRQCD switch
  Double_t fAlphai;	// 1/alpha(m_Z)
  Double_t fAlphas;	// alpha_s(m_Z)
  Double_t fMassW;	// m_W
  Double_t fMassZ;	// m_Z
  Double_t fMassHiggs;	// m_H
  Double_t fMassTop;	// m_t
  Int_t   fISHUFL[50];  // random number shuffler
public:
  TTBases(const char *name="TTBases", 
	     const char *title="TTbar  Bases");
  virtual ~TTBases(){} 

  void SetDeltaRoots(Double_t deltrs){ fDeltaRoots=deltrs; }
  Double_t GetDeltaRoots(){ return fDeltaRoots;}
  
  void SetRoots(Double_t rs){ fDeltaRoots=rs-2*fMassTop; }
  Double_t GetRoots(){ return (fDeltaRoots+2*fMassTop);}
  
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
   
  void SetVkmTB(Double_t vkmtb){ fVkmTB=vkmtb; }
  Double_t GetVkmTB(){ return fVkmTB;}
   
  void SetBetaH(Double_t betah){ fBetaH=betah; }
  Double_t GetBetaH(){ return fBetaH;}
   
  void SetNRQCD(Int_t nrqcd){ fNRQCD=nrqcd; }
  Int_t GetNRQCD(){ return fNRQCD;}
   
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

  ClassDef(TTBases,1)  // Bases for e+e- -> ttbar process

};

class TTSpring;

class TTSpringBuf : public JSFSpringBuf {
public:
//  <<+TTbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-TTbar>>
public:
  TTSpringBuf(const char *name="TTSpringBuf", 
	     const char *title="TTSpring test event buffer",
	     TTSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~TTSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(TTSpringBuf,1)  // TTSpring event buffer
};

class TTSpring : public JSFSpring {
public:
   TTSpring(const char *name="TTSpring", 
	      const char *title="TTSpring test",
             TTBases *bases=NULL);
   virtual ~TTSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(TTSpring,1)  // TTSpring class
};


#endif







