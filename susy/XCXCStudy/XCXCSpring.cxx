//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  XCXCSpring
//  
//  e+e- -> XCXC
//  
//  In this program, meanings of integration variables are as follows.
//
//  Definition of vairables
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : bremsstrahlng
//     ( 4) : X+ decay mode
//     ( 5) : X- decay mode
//     ( 6) : e- helicity
//     ( 7) : m(+ab)**2
//     ( 8) : m(+ac)**2
//     ( 9) : m(-ab)**2
//     (10) : m(-ac)**2
//     (11) : m(X+)**2
//     (12) : m(X-)**2
//     (13) : cos(theta_X-)
//     (14) : phi_X-
//     (15) : cos(theta_a)     in X+ rest frame
//     (16) : phi_a            in X+ rest frame
//     (17) : phi_b            in X+ rest frame
//     (18) : cos(theta_a)     in X- rest frame
//     (19) : phi_a            in X- rest frame
//     (20) : phi_b            in X- rest frame
//     (21) : E_beam spread of E-
//     (22) : E_beam spread of E+
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "XCXCSpring.h"

ClassImp(XCXCSpring)
ClassImp(XCXCSpringBuf)
ClassImp(XCXCBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
};

//_____________________________________________________________________________
XCXCSpring::XCXCSpring(const char *name, const char *title,
			 XCXCBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new XCXCSpringBuf("XCXCSpringBuf", 
	 "XCXCSpring event buffer", this);
  if( !bases ) { 
    XCXCBases *bs=new XCXCBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
XCXCSpring::~XCXCSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t XCXCSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void XCXCSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________

XCXCBases::XCXCBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("XCXCBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("XCXCBases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("XCXCBases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("XCXCBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("XCXCBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("XCXCBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 17;
  	fNWILD = 7;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 18;
  	fNWILD = 8;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 22;
  	fNWILD = 10;
  	fIOFF  = 0;
  } else {
  	printf("XCXCBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get XCXCbar specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"XCXCBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("XCXCBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("XCXCBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("XCXCBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("XCXCBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("XCXCBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("XCXCBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("XCXCBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("XCXCBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("XCXCBases.MassTop","170."),"%lg",&fMassTop);  
  sscanf(gJSF->Env()->GetValue("XCXCBases.m0","70."),"%lg",&fm0);
  sscanf(gJSF->Env()->GetValue("XCXCBases.mu","400."),"%lg",&fmu);
  sscanf(gJSF->Env()->GetValue("XCXCBases.M2","250."),"%lg",&fM2);
  sscanf(gJSF->Env()->GetValue("XCXCBases.tanb","+2."),"%lg",&ftanb);
  sscanf(gJSF->Env()->GetValue("XCXCBases.mA","-9999."),"%lg",&fmA);
  sscanf(gJSF->Env()->GetValue("XCXCBases.WidthChic1","-1."),"%lg",&fWidthChic1);

  fPrintInfo = gJSF->Env()->GetValue("XCXCBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("XCXCBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void XCXCBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->XCXC generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  XC_1 width = %g\n",fWidthChic1);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  return ;

}

//_____________________________________________________________________________
Double_t XCXCBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void XCXCBases::Userin()
{
//
//   Initialize User parameters for Bases
//
  JSFBases::Userin();  // Call JSFBases::Userin() for standard setup.

  // Copy class data member into common /usmprm/
  usmprm_.alfi = fAlphai;
  usmprm_.alfs = fAlphas;
  usmprm_.amsw = fMassW;
  usmprm_.amsz = fMassZ;
  usmprm_.amsh = fMassHiggs;
  usmprm_.amst = fMassTop;

  // Copy class data member into common /usmprm/
  ussprm_.am0  = fm0;
  ussprm_.amu  = fmu;
  ussprm_.am2  = fM2;
  ussprm_.tanb = ftanb;
  ussprm_.ama  = fmA;

  // Copy class data member into common /usrprm/
  usrprm_.sqrts  = fRoots;
  usrprm_.polebm = fPolElectron;
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;
  usrprm_.gamsw1 = fWidthChic1;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms
  Double_t ptmx = fRoots/2.;
  Double_t qmx  = fRoots/2;
  Xhinit( 1, -1.0, 1.0, 50,"cos(theta_X-)        ");
  Xhinit( 2,  0.0,360., 50,"phi_X-               ");
  Xhinit( 3, -1.0, 1.0, 50,"cos(theta_a+)        ");
  Xhinit( 4,  0.0,360., 50,"phi_b+               ");
  Xhinit( 5, -1.0, 1.0, 50,"cos(theta_a-)        ");
  Xhinit( 6,  0.0,360., 50,"phi_a-               ");
  Xhinit( 7,  0.0, qmx, 50,"M_ab+                ");
  Xhinit( 8,  0.0, qmx, 50,"M_ac+                ");
  Xhinit( 9,  0.0,ptmx, 50,"Missing PT           ");
  Xhinit(10,  0.0,180., 50,"Acop                 ");
  Xhinit(11,  0.0, qmx, 50,"M_ab-                ");
  Xhinit(12,  0.0, qmx, 50,"M_ac-                ");
  Xhinit(13,  1.0, 33., 32,"Hel. comb.           ");
  Xhinit(14,  0.0, qmx, 50,"E_e/mu               ");
  Xhinit(15, -1.0, 1.0, 50,"cos(theta_l-)        ");
  Xhinit(16,  0.0, qmx, 50,"E_qqbar              ");
  Xhinit(17, -1.0, 1.0, 50,"cos(theta_J-)        ");
  Xhinit(18,  1.0,  4.,  3,"Topology (LL,LJ,JJ)  ");
  Xhinit(19,  1.0, 13., 12,"Decay mode( X+ )     ");
  Xhinit(20,  1.0, 13., 12,"Decay mode( X- )     ");
  Xhinit(21,   .0,  1., 50,"sqrt(rs)/sqrt(root)  ");
  Xhinit(22, -1.0, 1.0, 50,"cos(theta_fdbar)     ");
  Xhinit(23, -1.0, 1.0, 50,"cos(theta_fd)        ");
  Xhinit(24,  0.0, qmx, 50,"M_X+                 ");
  Xhinit(25,  0.0, qmx, 50,"M_X-                 ");
  Xhinit(26, -1.0, 1.0, 50,"cos(theta_l+)_35     ");
  Xhinit(27, -1.0, 1.0, 50,"cos(theta_l-)_68     ");
  return ;
}

//_____________________________________________________________________________
void XCXCBases::Userout()
{
  usrout_();
}






















