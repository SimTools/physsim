//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  NNZSpring
//  
//  e+e- -> nu nu Z
//  
//  Integration variables.
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : bremsstrahlng
//     ( 4) : helicity
//     ( 5) : xi   (E_nu)
//     ( 6) : m(Z)**2
//     ( 7) : eta  (cos_nu)
//     ( 8) : zeta (cos_nb) in (nb Z) frame
//     ( 9) : phi_nu
//     (10) : phi_nb  	   in (nb Z) frame
//     (11) : cos_fb       in Z      frame
//     (12) : phi_f        in Z      frame
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
#include "NNZSpring.h"

ClassImp(NNZSpring)
ClassImp(NNZSpringBuf)
ClassImp(NNZBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
};

//_____________________________________________________________________________
NNZSpring::NNZSpring(const char *name, const char *title,
			 NNZBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new NNZSpringBuf("NNZSpringBuf", 
	 "NNZSpring event buffer", this);
  if( !bases ) { 
    NNZBases *bs=new NNZBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
NNZSpring::~NNZSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t NNZSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void NNZSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }


//_____________________________________________________________________________
NNZBases::NNZBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("NNZBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("NNZBases.ACC1","0.2"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("NNZBases.ACC2","0.1"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("NNZBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("NNZBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("NNZBases.NCALL","80000"),"%d",&fNCALL);

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
  	printf("NNZBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }
#ifdef ZEROWIDTH
	fNWILD -= 3;
#endif

//
// Get NNZ specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"NNZBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("NNZBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("NNZBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("NNZBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("NNZBases.ZModesLo","1"),"%d",&fZModesLo);
  sscanf(gJSF->Env()->GetValue("NNZBases.ZModesHi","12"),"%d",&fZModesHi);
  sscanf(gJSF->Env()->GetValue("NNZBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("NNZBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("NNZBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("NNZBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("NNZBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("NNZBases.MassTop","170."),"%lg",&fMassTop);  

  fPrintInfo = gJSF->Env()->GetValue("NNZBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("NNZBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void NNZBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->NNZ generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  Z0 Decey Mode Lo =%d\n",fZModesLo);
  printf("                Hi =%d\n",fZModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  return ;

}

//_____________________________________________________________________________
Double_t NNZBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void NNZBases::Userin()
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
  usrprm_.sqrts  = fRoots;
  usrprm_.polebm = fPolElectron;
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;
  usrprm_.imdmin = fZModesLo;
  usrprm_.imdmax = fZModesHi;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms
  Xhinit( 1, -1.0, 1.0, 50,"cos(n)        ");
  Xhinit( 2,  0.0,360., 50,"phi_n         ");
  Xhinit( 3,   0., 1.0, 50,"E_n/E_bm      ");
  Xhinit( 4,  60.,110., 50,"m(Z)          ");
  Xhinit( 5, -1.0, 1.0, 50,"cos(nb)       ");
  Xhinit( 6,  0.0,360., 50,"phi_nb        ");
  Xhinit( 7, -1.0, 1.0, 50,"cos(theta_fb) ");
  Xhinit( 8,  0.0,360., 50,"phi_fb        ");
  Xhinit( 9,  0.0, 1.0, 50,"RS/ROOTS      ");
  Xhinit(10,  1.0,13.0, 12,"Z decay mode  ");
  Xhinit(11,  1.0,17.0, 16,"helicity      ");
  Xhinit(12, -10.,20.0, 40,"eta_n         ");
  Xhinit(13,   0., 1.0, 50,"PT_Z/E_bm     ");
  Xhinit(14,   0., 1.0, 50,"E_n/E_bm      ");
  Xhinit(15,   0., 1.0, 50,"E_nb/E_bm     ");
  Xhinit(16,   0., 1.0, 50,"E_Z/E_bm      ");
  return ;
}

//_____________________________________________________________________________
void NNZBases::Userout()
{
  usrout_();
}






















