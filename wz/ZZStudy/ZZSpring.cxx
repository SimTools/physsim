//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  ZZSpring
//  
//  e+e- -> ZZ
//  
//  Integration variables.
//    Z( 1) : e- beam
//     ( 2) : e- beam gaussian spread
//     ( 3) : e+ beam
//     ( 4) : e+ beam gaussian spread
//     ( 5) : bremsstrahlng
//     ( 6) : e- helicity
//     ( 7) : helicity combination for final states.
//     ( 8) : m(Z)**2
//     ( 9) : m(Z)**2
//     (10) : cos(theta_Z)
//     (11) : phi_Z
//     (12) : cos(theta_fb)     in Z rest frame
//     (13) : phi_fb            in Z rest frame
//     (14) : cos(theta_fb)     in Z rest frame
//     (15) : phi_fb            in Z rest frame
//     (16) : final state combination
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "ZZSpring.h"

ClassImp(ZZSpring)
ClassImp(ZZSpringBuf)
ClassImp(ZZBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
};

//_____________________________________________________________________________
ZZSpring::ZZSpring(const char *name, const char *title,
			 ZZBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new ZZSpringBuf("ZZSpringBuf", 
	 "ZZSpring event buffer", this);
  if( !bases ) { 
    ZZBases *bs=new ZZBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
ZZSpring::~ZZSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t ZZSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void ZZSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
ZZBases::ZZBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("ZZBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("ZZBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("ZZBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("ZZBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("ZZBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("ZZBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 11;
  	fNWILD = 4;
  	fIOFF  = 5;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 12;
  	fNWILD = 5;
  	fIOFF  = 4;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 16;
  	fNWILD = 9;
  	fIOFF  = 0;
  } else {
  	printf("ZZBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get ZZ specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"ZZBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("ZZBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("ZZBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("ZZBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("ZZBases.Z1ModesLo","1"),"%d",&fZ1ModesLo);
  sscanf(gJSF->Env()->GetValue("ZZBases.Z1ModesHi","12"),"%d",&fZ1ModesHi);
  sscanf(gJSF->Env()->GetValue("ZZBases.Z2ModesLo","1"),"%d",&fZ2ModesLo);
  sscanf(gJSF->Env()->GetValue("ZZBases.Z2ModesHi","12"),"%d",&fZ2ModesHi);
  sscanf(gJSF->Env()->GetValue("ZZBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("ZZBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("ZZBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("ZZBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("ZZBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("ZZBases.MassTop","170."),"%lg",&fMassTop);  

  fPrintInfo = gJSF->Env()->GetValue("ZZBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("ZZBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void ZZBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->ZZ generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  W- Decey Mode Lo =%d\n",fZ1ModesLo);
  printf("                Hi =%d\n",fZ1ModesHi);
  printf("  W+ Decey Mode Lo =%d\n",fZ2ModesLo);
  printf("                Hi =%d\n",fZ2ModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  return ;

}

//_____________________________________________________________________________
Double_t ZZBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void ZZBases::Userin()
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

  // Copy class data member into common /usrprm/
  usrprm_.sqrts = fRoots;
  usrprm_.polebm = fPolElectron;
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;
  usrprm_.imd1lo = fZ1ModesLo;
  usrprm_.imd1hi = fZ1ModesHi;
  usrprm_.imd2lo = fZ2ModesLo;
  usrprm_.imd2hi = fZ2ModesHi;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

  Xhinit( 1,  0.0,  1., 50,"rs/roots      ");
  Xhinit( 2,  60.,110., 50,"m(Z1)         ");
  Xhinit( 3,  60.,110., 50,"m(Z2)         ");
  Xhinit( 4,  0.0,  1., 50,"miss/roots    ");
  Xhinit( 5, -1.0, 1.0, 50,"cos(theta_Z)  ");
  Xhinit( 6,  0.0,360., 50,"phi_Z         ");
  Xhinit( 7, -1.0, 1.0, 50,"cos(theta_f1) ");
  Xhinit( 8,  0.0,360., 50,"phi_f1        ");
  Xhinit( 9, -1.0, 1.0, 50,"cos(theta_f2) ");
  Xhinit(10,  0.0,360., 50,"phi_f2        ");
  Xhinit(11,  1.0,  9.,  8,"helicity      ");
  Xhinit(12,  1.0, 25., 24,"Z decay mode  ");
  return ;
}

//_____________________________________________________________________________
void ZZBases::Userout()
{
  usrout_();
}



















