//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  ENWSpring
//  
//  e+e- -> e- nubar W+
//  
//  Integration variables.
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : bremsstrahlng
//     ( 4) : helicity
//     ( 5) : xi   (E_e)
//     ( 6) : m(W)**2
//     ( 7) : eta  (cos_e)
//     ( 8) : zeta (cos_n) in (nu W) frame
//     ( 9) : phi_e
//     (10) : phi_n   	   in (nu W) frame
//     (11) : cos_fd       in W+     frame
//     (12) : phi_fd       in W+     frame
//     (13) : final state combination.
//     (14) : E_beam spread of E-
//     (15) : E_beam spread of E+
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "ENWSpring.h"

ClassImp(ENWSpring)
ClassImp(ENWSpringBuf)
ClassImp(ENWBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
};

//_____________________________________________________________________________
ENWSpring::ENWSpring(const char *name, const char *title,
			 ENWBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new ENWSpringBuf("ENWSpringBuf", 
	 "ENWSpring event buffer", this);
  if( !bases ) { 
    ENWBases *bs=new ENWBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
ENWSpring::~ENWSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t ENWSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void ENWSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
ENWBases::ENWBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("ENWBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("ENWBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("ENWBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("ENWBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("ENWBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("ENWBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 10;
  	fNWILD = 5;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 11;
  	fNWILD = 6;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 15;
  	fNWILD = 8;
  	fIOFF  = 0;
  } else {
  	printf("ENWBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get ENW specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"ENWBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("ENWBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("ENWBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("ENWBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("ENWBases.WpModesLo","1"),"%d",&fWpModesLo);
  sscanf(gJSF->Env()->GetValue("ENWBases.WpModesHi","12"),"%d",&fWpModesHi);
  sscanf(gJSF->Env()->GetValue("ENWBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("ENWBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("ENWBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("ENWBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("ENWBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("ENWBases.MassTop","170."),"%lg",&fMassTop);  

  fPrintInfo = gJSF->Env()->GetValue("ENWBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("ENWBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void ENWBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->ENW generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  W+ Decey Mode Lo =%d\n",fWpModesLo);
  printf("                Hi =%d\n",fWpModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  return ;

}

//_____________________________________________________________________________
Double_t ENWBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void ENWBases::Userin()
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
  usrprm_.imdmin = fWpModesLo;
  usrprm_.imdmax = fWpModesHi;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

  Xhinit( 1, -1.0, 1.0, 50,"cos(e)        ");
  Xhinit( 2,  0.0,360., 50,"phi_e         ");
  Xhinit( 3,   0., 1.0, 50,"E_e/E_bm      ");
  Xhinit( 4,  60.,110., 50,"m(W)          ");
  Xhinit( 5, -1.0, 1.0, 50,"cos(nu)       ");
  Xhinit( 6,  0.0,360., 50,"phi_nu        ");
  Xhinit( 7, -1.0, 1.0, 50,"cos(theta_fd) ");
  Xhinit( 8,  0.0,360., 50,"phi_fd        ");
  Xhinit( 9,  0.0, 1.0, 50,"RS/ROOTS      ");
  Xhinit(10,  1.0,13.0, 12,"W decay mode  ");
  Xhinit(11,  1.0,17.0, 16,"CP/helicity   ");
  Xhinit(12, -10.,20.0, 40,"eta_e         ");
  Xhinit(13,   0., 1.0, 50,"PT_W/E_bm     ");
  Xhinit(14,   0., 1.0, 50,"E_e/E_bm      ");
  Xhinit(15,   0., 1.0, 50,"E_nu/E_bm     ");
  Xhinit(16,   0., 1.0, 50,"E_W/E_bm      ");
  return ;
}

//_____________________________________________________________________________
void ENWBases::Userout()
{
  usrout_();
}



















