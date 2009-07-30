#ifndef __TTBBSpring__
#define __TTBBSpring__

//////////////////////////////////////////////////////////////////////////
//
// TTBBSpring
//
// e+e- -> ttbar gluon -> ttbar bbar
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
//  COMMONs for TTBB calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts, polebm, sgmebm;
  Int_t    isrb;
  Int_t    imd1lo;	// W- decay mode ID.
  Int_t    imd1hi;	// W- decay mode ID.
  Int_t    imd2lo;	// W+ decay mode ID.
  Int_t    imd2hi;	// W+ decay mode ID.
  Float_t  polpbm;      // positron polarization.
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class TTBBBases  : public JSFBases {
protected:
  Double_t fRoots;	// CM Energy (GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fPolPositron;	// e- polarization
  Double_t fSigmaEbeam;	// beam energy spread
  Int_t    fISRBM;	// Flag for ISR & Beamstrahlung
  Int_t    fWmModesLo;	// W- decay mode ID.
  Int_t    fWmModesHi;	// W- decay mode ID.
  Int_t    fWpModesLo;	// W+ decay mode ID.
  Int_t    fWpModesHi;	// W+ decay mode ID.
  Double_t fAlphai;	// 1/alpha(m_Z)
  Double_t fAlphas;	// alpha_s(m_Z)
  Double_t fMassW;	// m_W
  Double_t fMassZ;	// m_Z
  Double_t fMassHiggs;	// m_H
  Double_t fMassTop;	// m_t
  Int_t   fISHUFL[50];  // random number shuffler
public:
  TTBBBases(const char *name="TTBBBases", 
	     const char *title="TTBBbar  Bases");
  virtual ~TTBBBases(){} 

  void SetRoots(Double_t sqrts){ fRoots=sqrts; }
  Double_t GetRoots(){ return fRoots;}
  
  void SetPolElectron(Double_t polebm){ fPolElectron=polebm; }
  Double_t GetPolElectron(){ return fPolElectron;}

  void SetPolPositron(Double_t polpbm){ fPolPositron=polpbm; }
  Double_t GetPolPositron(){ return fPolPositron;}
  
  void SetSigmaEbeam(Double_t sgmebm){ fSigmaEbeam=sgmebm; }
  Double_t GetSigmaEbeam(){ return fSigmaEbeam;}
  
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

  ClassDef(TTBBBases,1)  // Bases for e+e- -> ttbar process

};

class TTBBSpring;

class TTBBSpringBuf : public JSFSpringBuf {
public:
//  <<+TTBBbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-TTBBbar>>
public:
  TTBBSpringBuf(const char *name="TTBBSpringBuf", 
	     const char *title="TTBBSpring test event buffer",
	     TTBBSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~TTBBSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(TTBBSpringBuf,1)  // TTBBSpring event buffer
};

class TTBBSpring : public JSFSpring {
public:
   TTBBSpring(const char *name="TTBBSpring", 
	      const char *title="TTBBSpring test",
             TTBBBases *bases=NULL);
   virtual ~TTBBSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(TTBBSpring,1)  // TTBBSpring class
};


#endif







