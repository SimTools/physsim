#ifndef __XN2XN2Spring__
#define __XN2XN2Spring__

//////////////////////////////////////////////////////////////////////////
//
// XN2XN2Spring
//
// e+e- -> XN2XN2
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
//  COMMONs for SS parameters, see USSPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  am0, amu, am2, tanb, ama, asft[3];
} COMMON_USSPRM;             //  Common /USSPRM/

extern COMMON_USSPRM ussprm_;

// =====================================================================
//  COMMONs for XN2XN2 calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts;
  Float_t  polebm;
  Float_t  sgmebm;
  Float_t  gamsz2;
  Int_t    isrb;
  Int_t    imd1lo;      // 1st Z1 decay mode ID.
  Int_t    imd1hi;      // lst Z1 decay mode ID.
  Int_t    imd2lo;      // 1st Z2 decay mode ID.
  Int_t    imd2hi;      // lst Z2 decay mode ID.
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class XN2XN2Bases  : public JSFBases {
protected:
  Double_t fRoots;		// CM Energy(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;		// beam energy spread
  Int_t    fISRBM;		// Flag for ISR & Beamstrahlung
  Int_t    fZ1ModesLo;  	// 1st Z1 decay mode ID.
  Int_t    fZ1ModesHi;  	// lst Z1 decay mode ID.
  Int_t    fZ2ModesLo;  	// 1st Z2 decay mode ID.
  Int_t    fZ2ModesHi;  	// lst Z2 decay mode ID.
  Double_t fAlphai;		// 1/alpha(m_Z)
  Double_t fAlphas;		// alpha_s(m_Z)
  Double_t fMassW;		// m_W
  Double_t fMassZ;		// m_Z
  Double_t fMassHiggs;		// m_H
  Double_t fMassTop;		// m_t
  Double_t fm0;			// m0 (GeV)
  Double_t fmu;			// mu (GeV)
  Double_t fM2;			// M2 (GeV)
  Double_t ftanb;		// tan(beta)
  Double_t fmA;			// mA (GeV)
  Double_t fAsoft[3];           // (A_tau, A_t, A_b) (GeV)
  Double_t fWidthChin2;		// XN2 width (GeV)
  Int_t    fISHUFL[50];  	// random number shuffler

public:
  XN2XN2Bases(const char *name="XN2XN2Bases", 
	     const char *title="XN2XN2bar  Bases");
  virtual ~XN2XN2Bases(){} 

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
   
  void SetM0(Double_t m){ fm0 = m; }
  Double_t GetM0(){ return fm0;}
   
  void SetMu(Double_t m){ fmu = m; }
  Double_t GetMu(){ return fmu;}
    
  void SetM2(Double_t m){ fM2 = m; }
  Double_t GetM2(){ return fM2;}
  
  void SetTanb(Double_t m){ ftanb = m; }
  Double_t GetTanb(){ return ftanb;}
  
  void SetMA(Double_t m){ fmA = m; }
  Double_t GetMA(){ return fmA;}

  void SetAsoft(Int_t i, Double_t m){ fAsoft[i] = m; }
  Double_t GetAsoft(Int_t i) { return fAsoft[i];}
  
  void SetWidthChin2(Double_t g){ fWidthChin2 = g; }
  Double_t GetWidthChin2(){ return fWidthChin2;}

  void Userin();   // Bases user initialization
  void Userout();  // Bases user output 
  Double_t Func(Double_t x[]); // Bases integration function.
  Double_t Func(); // Bases integration function.
  void PrintParameters(); // Print parameters

  ClassDef(XN2XN2Bases,1)  // Bases for e+e- -> XN2XN2bar process

};

class XN2XN2Spring;

class XN2XN2SpringBuf : public JSFSpringBuf {
public:
//  <<+XN2XN2bar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-XN2XN2bar>>
public:
  XN2XN2SpringBuf(const char *name="XN2XN2SpringBuf", 
	     const char *title="XN2XN2Spring test event buffer",
	     XN2XN2Spring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~XN2XN2SpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(XN2XN2SpringBuf,1)  // XN2XN2Spring event buffer
};

class XN2XN2Spring : public JSFSpring {
public:
   XN2XN2Spring(const char *name="XN2XN2Spring", 
	      const char *title="XN2XN2Spring test",
             XN2XN2Bases *bases=NULL);
   virtual ~XN2XN2Spring();

   virtual Bool_t Initialize(); // 

   ClassDef(XN2XN2Spring,1)  // XN2XN2Spring class
};


#endif
