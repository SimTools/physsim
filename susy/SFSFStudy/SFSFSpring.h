#ifndef __SFSFSpring__
#define __SFSFSpring__

//////////////////////////////////////////////////////////////////////////
//
// SFSFSpring
//
// e+e- -> SFSFbar
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
//  COMMONs for SFSF calculation, see USRPRM.inc for more details.
// =====================================================================
typedef struct {
  Float_t  sqrts;
  Float_t  polebm;
  Float_t  sgmebm;
  Int_t	   ignsf;
  Int_t    ihndm;
  Int_t    ihndp;
  Int_t    idotu;
  Float_t  htum;
  Float_t  htup;
  Int_t    isrb;
} COMMON_USRPRM;             //  Common /USRPRM/

extern COMMON_USRPRM usrprm_;

class SFSFBases  : public JSFBases {
protected:
  Double_t fRoots;		// CM Energy(GeV)
  Double_t fPolElectron;	// e- polarization
  Double_t fSigmaEbeam;		// beam energy spread
  Int_t    fISRBM;		// Flag for ISR & Beamstrahlung
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
  Double_t fAsoft[3];		// (A_tau, A_t, A_b) (GeV)
  Int_t    fSfGeneration;	// sfermion generation
  Int_t    fHandMinus;		// sfermion- handedness (1,2) for (L,R) 
  Int_t    fHandPlus;		// sfermion+ handedness (1,2) for (L,R) 
  Int_t    fIDoTau;		// do special treatment for stau
  Double_t fHelTauMinus;	// tau- polarization
  Double_t fHelTauPlus;		// tau+ polarization
  Int_t    fISHUFL[50];  	// random number shuffler

public:
  SFSFBases(const char *name="SFSFBases", 
	     const char *title="SFSFbar  Bases");
  virtual ~SFSFBases(){} 

  void SetRoots(Double_t deltrs){ fRoots=deltrs; }
  Double_t GetRoots(){ return fRoots;}
  
  void SetPolElectron(Double_t polebm){ fPolElectron=polebm; }
  Double_t GetPolElectron(){ return fPolElectron;}
  
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
  
  void SetSfGeneration(Int_t i){ fSfGeneration=i; }
  Int_t GetSfGeneration(){ return fSfGeneration;}
  
  void SetHandMinus(Int_t i){ fHandMinus=i; }
  Int_t GetHandMinus(){ return fHandMinus;}
  
  void SetHandPlus(Int_t i){ fHandPlus=i; }
  Int_t GetHandPlus(){ return fHandPlus;}
  
  void SetIDoTau(Int_t i){ fIDoTau=i; }
  Int_t GetIDoTau(){ return fIDoTau;}
  
  void SetHelTauMinus(Double_t m){ fHelTauMinus = m; }
  Double_t GetHelTauMinus(){ return fHelTauMinus;}
  
  void SetHelTauPlus(Double_t m){ fHelTauPlus = m; }
  Double_t GetHelTauPlus(){ return fHelTauPlus;}

  void Userin();   // Bases user initialization
  void Userout();  // Bases user output 
  Double_t Func(Double_t x[]); // Bases integration function.
  Double_t Func(); // Bases integration function.
  void PrintParameters(); // Print parameters

  ClassDef(SFSFBases,1)  // Bases for e+e- -> SFSFbar process

};

class SFSFSpring : public JSFSpring {
public:
   SFSFSpring(const char *name="SFSFSpring", 
	      const char *title="SFSFSpring test",
             SFSFBases *bases=NULL);
   virtual ~SFSFSpring();

   virtual Bool_t Initialize(); // 

   ClassDef(SFSFSpring,1)  // SFSFSpring class
};


class SFSFSpringBuf : public JSFSpringBuf {
public:
//  <<+SFSFbar>>
   Double_t fX[2];  // Two bases parameters, costh and phi.
//  <<-SFSFbar>>
public:
  SFSFSpringBuf(const char *name="SFSFSpringBuf", 
	     const char *title="SFSFSpring test event buffer",
	     SFSFSpring *spring=NULL)
	     : JSFSpringBuf(name,title,spring) {} 
  virtual ~SFSFSpringBuf(){}
  virtual void Spevnt(Int_t &iret);

  ClassDef(SFSFSpringBuf,1)  // SFSFSpring event buffer
};


#endif








