//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  SFSFSpring
//  
//  e+e- -> SFSFbar
//  
//  In this program, meanings of integration variables are as follows.
//
//  Definition of vairables
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : bremsstrahlng
//     ( 4) : e- helicity
//     ( 5) : m(X-)**2
//     ( 6) : m(X+)**2
//     ( 7) : cos(theta_X-)
//     ( 8) : phi_X-
//     ( 9) : cos(theta_l-)    in X- rest frame
//     (10) : phi_l-           in X- rest frame
//     (11) : cos(theta_l+)    in X+ rest frame
//     (12) : phi_l+           in X+ rest frame
//     (13) : e- beam energy spread
//     (14) : e+ beam energy spread
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "SFSFSpring.h"

ClassImp(SFSFSpring)
ClassImp(SFSFSpringBuf)
ClassImp(SFSFBases)

extern "C" {
extern void usrout_();
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
};

//_____________________________________________________________________________
SFSFSpring::SFSFSpring(const char *name, const char *title,
			 SFSFBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new SFSFSpringBuf("SFSFSpringBuf", 
	 "SFSFSpring event buffer", this);
  if( !bases ) { 
    SFSFBases *bs=new SFSFBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
SFSFSpring::~SFSFSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t SFSFSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void SFSFSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________

SFSFBases::SFSFBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//
// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("SFSFBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("SFSFBases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("SFSFBases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("SFSFBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("SFSFBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("SFSFBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 9;
  	fNWILD = 3;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 10;
  	fNWILD = 4;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 14;
  	fNWILD = 6;
  	fIOFF  = 0;
  } else {
  	printf("SFSFBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get SFSFbar specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"SFSFBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("SFSFBases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("SFSFBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("SFSFBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("SFSFBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("SFSFBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("SFSFBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("SFSFBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("SFSFBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("SFSFBases.MassTop","170."),"%lg",&fMassTop);  
  sscanf(gJSF->Env()->GetValue("SFSFBases.m0","70."),"%lg",&fm0);
  sscanf(gJSF->Env()->GetValue("SFSFBases.mu","400."),"%lg",&fmu);
  sscanf(gJSF->Env()->GetValue("SFSFBases.M2","250."),"%lg",&fM2);
  sscanf(gJSF->Env()->GetValue("SFSFBases.tanb","+2."),"%lg",&ftanb);
  sscanf(gJSF->Env()->GetValue("SFSFBases.mA","-9999."),"%lg",&fmA);
  sscanf(gJSF->Env()->GetValue("SFSFBases.SfGeneration","2"),"%d",&fSfGeneration);
  sscanf(gJSF->Env()->GetValue("SFSFBases.HandMinus","2"),"%d",&fHandMinus);
  sscanf(gJSF->Env()->GetValue("SFSFBases.HandPlus","2"),"%d",&fHandPlus);
  sscanf(gJSF->Env()->GetValue("SFSFBases.IDoTau","0"),"%d",&fIDoTau);
  sscanf(gJSF->Env()->GetValue("SFSFBases.HelTauMinus","0."),"%lg",&fHelTauMinus);
  sscanf(gJSF->Env()->GetValue("SFSFBases.HelTauPlus","0."),"%lg",&fHelTauPlus);

  fPrintInfo = gJSF->Env()->GetValue("SFSFBases.PrintInfo",kTRUE);
  fPrintHist = gJSF->Env()->GetValue("SFSFBases.PrintHist",kTRUE);

}


//_____________________________________________________________________________
void SFSFBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->SFSF generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  Generation = %d\n",fSfGeneration);
  printf("  Handedness for sf- = %d\n",fHandMinus);
  printf("  Handedness for sf+ = %d\n",fHandPlus);
  printf("  Flag for special stau treatment =%d\n",fIDoTau);
  printf("       = 0 ; Off\n");
  printf("       = 1 ; On\n");
  printf("  tau- Pol = %g\n",fHelTauMinus);
  printf("  tau+ Pol = %g\n",fHelTauPlus);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  return ;

}

//_____________________________________________________________________________
Double_t SFSFBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void SFSFBases::Userin()
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
  usrprm_.ignsf  = fSfGeneration;
  usrprm_.ihndm  = fHandMinus;
  usrprm_.ihndp  = fHandPlus;
  usrprm_.idotu  = fIDoTau;
  usrprm_.htum   = fHelTauMinus;
  usrprm_.htup   = fHelTauPlus;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

      Double_t qmx  = fRoots/2.;
      Double_t ptmx = 150.;
      Double_t elmx = fRoots/2.;
      Xhinit( 1, -1.0, 1.0, 50,"cos(theta_X-)      ");
      Xhinit( 2,  0.0,360., 50,"phi_X-             ");
      Xhinit( 3, -1.0, 1.0, 50,"cos(theta_l-)      ");
      Xhinit( 4,  0.0,360., 50,"phi_l-             ");
      Xhinit( 5, -1.0, 1.0, 50,"cos(theta_l+)      ");
      Xhinit( 6,  0.0,360., 50,"phi_l+             ");
      Xhinit( 7,  0.0, qmx, 50,"M_X-               ");
      Xhinit( 8,  0.0, qmx, 50,"M_X+               ");
      Xhinit( 9,  0.0,ptmx, 50,"Missing PT         ");
      Xhinit(10,  0.0,180., 50,"Acop               ");
      Xhinit(11,  0.0,elmx, 50,"E_l-               ");
      Xhinit(12,  0.0,elmx, 50,"E_l+               ");
      Xhinit(13, -1.0, 1.0, 50,"cos(theta_l-)_lab  ");
      Xhinit(14, -1.0, 1.0, 50,"-cos(theta_l+)_lab ");
      Xhinit(15, -1.0, 1.0, 50,"cos(theta_l)_lab   ");
  return ;
}

//_____________________________________________________________________________
void SFSFBases::Userout()
{
  usrout_();
}

















