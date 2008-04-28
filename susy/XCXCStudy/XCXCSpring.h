#ifndef __XCXCSpring__
#define __XCXCSpring__

//////////////////////////////////////////////////////////////////////////
//
// XCXCSpring
//
// e+e- -> X+X-bar
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
//  COMMONs for SS parameters, see USSPRMP.inc for more details.
// =====================================================================
typedef struct {
  Float_t  am0, amu, am2, tanb, ama, asft[3];
  Float_t  am1, am3;
  Int_t    isgut;
} COMMON_USSPRMP;             //  Common /USSPRMP/

extern COMMON_USSPRMP ussprmp_;

// =====================================================================
//  COMMONs for XCXC calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts;
  Float_t  polebm;
  Float_t  sgmebm;
  Float_t  gamsw1;
  Int_t    isrb;
  Int_t    imd1lo;      // 1st W- decay mode ID.
  Int_t    imd1hi;      // lst W- decay mode ID.
  Int_t    imd2lo;      // 1st W+ decay mode ID.
  Int_t    imd2hi;      // lst W+ decay mode ID.
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class XCXCBases  : public JSFBases {
protected:
  Double_t fRoots;		// CM Energy(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;		// beam energy spread
  Int_t    fISRBM;		// Flag for ISR & Beamstrahlung
  Int_t    fWmModesLo;  	// 1st W- decay mode ID.
  Int_t    fWmModesHi;  	// lst W- decay mode ID.
  Int_t    fWpModesLo;  	// 1st W+ decay mode ID.
  Int_t    fWpModesHi;  	// lst W+ decay mode ID.
  Int_t    fGUT;                // (0,1)=(No GUT, GUT)
  Double_t fAlphai;		// 1/alpha(m_Z)
  Double_t fAlphas;		// alpha_s(m_Z)
  Double_t fMassW;		// m_W
  Double_t fMassZ;		// m_Z
  Double_t fMassHiggs;		// m_H
  Double_t fMassTop;		// m_t
  Double_t fm0;			// m0 (GeV)
  Double_t fmu;			// mu (GeV)
  Double_t fM1;			// M1 (GeV)
  Double_t fM2;			// M2 (GeV)
  Double_t ftanb;		// tan(beta)
  Double_t fmA;			// mA (GeV)
  Double_t fAsoft[3];           // (A_tau, A_t, A_b) (GeV)
  Double_t fWidthChic1;		// XC1 width (GeV)
  Int_t    fISHUFL[50];  	// random number shuffler

public:
  XCXCBases(const char *name="XCXCBases", 
	     const char *title="XCXCbar  Bases");
  virtual ~XCXCBases(){} 

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
    
  void SetGUT(Double_t m){ fGUT = m; }
  Int_t GetGUT(){ return fGUT;}
    
  void SetM1(Double_t m){ fM1 = m; }
  Double_t GetM1(){ return fM1;}
    
  void SetM2(Double_t m){ fM2 = m; }
  Double_t GetM2(){ return fM2;}
  
  void SetTanb(Double_t m){ ftanb = m; }
  Double_t GetTanb(){ return ftanb;}
  
  void SetMA(Double_t m){ fmA = m; }
  Double_t GetMA(){ return fmA;}

  void SetAsoft(Int_t i, Double_t m){ fAsoft[i] = m; }
  Double_t GetAsoft(Int_t i) { return fAsoft[i];}
  
  void SetWidthChic1(Double_t g){ fWidthChic1 = g; }
  Double_t GetWidthChic1(){ return fWidthChic1;}

  void Userin();   // Bases user initialization
  void Userout();  // Bases user output 
  Double_t Func(Double_t x[]); // Bases integration function.
  Double_t Func(); // Bases integration function.
  void PrintParameters(); // Print parameters

  ClassDef(XCXCBases,1)  // Bases for e+e- -> XCXCbar process

};

class XCXCSpring;

class XCXCSpringBuf : public JSFSpringBuf {
public:
//  <<+XCXCbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-XCXCbar>>
public:
  XCXCSpringBuf(const char *name="XCXCSpringBuf", 
	     const char *title="XCXCSpring test event buffer",
	     XCXCSpring *spring=NULL)
	     : JSFSpringBuf(name,title,(JSFSpring*)spring) {} 
  virtual ~XCXCSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(XCXCSpringBuf,1)  // XCXCSpring event buffer
};

class XCXCSpring : public JSFSpring {
public:
   XCXCSpring(const char *name="XCXCSpring", 
	      const char *title="XCXCSpring test",
             XCXCBases *bases=NULL);
   virtual ~XCXCSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(XCXCSpring,1)  // XCXCSpring class
};


#endif









