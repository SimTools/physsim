#ifndef __ZZSpring__
#define __ZZSpring__

//////////////////////////////////////////////////////////////////////////
//
// ZHSpring
//
// e+e- -> ZH
//
//////////////////////////////////////////////////////////////////////////

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
  Int_t nZH;
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
//  COMMONs for ZH calculation, see USRPRM.inc for more details.
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

class ZHBases  : public JSFBases {
protected:
  Double_t fRoots;	// CM Energy - 2m_t(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;	// beam energy spread
  Int_t    fISRBM;	// Flag for ISR & Beamstrahlung
  Int_t    fZ0ModesLo;	// 1st Z0 decay mode ID.
  Int_t    fZ0ModesHi;	// lst Z0 decay mode ID.
  Int_t    fH0ModesLo;	// 1st H0 decay mode ID (ffbar only).
  Int_t    fH0ModesHi;	// lst H0 decay mode ID (ffbar only).
  Double_t fAlphai;	// 1/alpha(m_Z)
  Double_t fAlphas;	// alpha_s(m_Z)
  Double_t fMassW;	// m_W
  Double_t fMassZ;	// m_Z
  Double_t fMassHiggs;	// m_H
  Double_t fMassTop;	// m_t
  Int_t   fISHUFL[50];  // random number shuffler
public:
  ZHBases(const char *name="ZHBases", 
	     const char *title="ZHbar  Bases");
  virtual ~ZHBases(){} 

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
   
  void SetH0ModesLo(Int_t i){ fH0ModesLo=i; }
  Int_t GetH0ModesLo(){ return fH0ModesLo;}
   
  void SetH0ModesHi(Int_t i){ fH0ModesHi=i; }
  Int_t GetH0ModesHi(){ return fH0ModesHi;}
   
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
  void PrintParameters(); // Print parameters

  ClassDef(ZHBases,1)  // Bases for e+e- -> ZH process

};

class ZHSpring;

class ZHSpringBuf : public JSFSpringBuf {
public:
//  <<+ZHbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-ZHbar>>
public:
  ZHSpringBuf(const char *name="ZHSpringBuf", 
	     const char *title="ZHSpring test event buffer",
	     ZHSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~ZHSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(ZHSpringBuf,1)  // ZHSpring event buffer
};

class ZHSpring : public JSFSpring {
public:
   ZHSpring(const char *name="ZHSpring", 
	      const char *title="ZHSpring test",
             ZHBases *bases=NULL);
   virtual ~ZHSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(ZHSpring,1)  // ZHSpring class
};


#endif

