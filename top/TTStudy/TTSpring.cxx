//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  TTSpring
//  
//  e+e- -> ttbar
//  
//  In this program, meanings of integration variables are as follows.
//
//  Definition of vairables
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : bremsstrahlng
//     ( 4) : e- helicity
//     ( 5) : m(t_bar)**2
//     ( 6) : m(t)**2
//     ( 7) : m(W-)**2
//     ( 8) : m(W+)**2
//     ( 9) : cos(theta_t)
//     (10) : phi_t
//     (11) : cos(theta_q_bar) in t_bar rest frame
//     (12) : phi_q_bar        in t_bar rest frame
//     (13) : cos(theta_f)     in W- rest frame
//     (14) : phi_f            in W- rest frame
//     (15) : cos(theta_q)     in t rest frame
//     (16) : phi_q            in t rest frame
//     (17) : cos(theta_f_bar) in W+ rest frame
//     (18) : phi_f_bar        in W+ rest frame
//     (19) : final state combination.
//     (20) : e- beam gaussian spread
//     (21) : e+ beam gaussian spread
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "TTSpring.h"

ClassImp(TTSpring)
ClassImp(TTSpringBuf)
ClassImp(TTBases)

// Variable conversion.

#define fNDIM  ndim 
#define fNWILD nwild
#define fIG    ig
#define fXL    xl
#define fXU    xu
#define fNCALL ncall
#define fACC1  acc1
#define fACC2  acc2
#define fITMX1 itmx1
#define fITMX2 itmx2

#define Xhinit(id,xlo,xhi,n,title) H1Init(id,title,n,xlo,xhi)
#define Dhinit(id,xlo,xhi,nx,ylo,yhi,ny,title) H2Init(id,title,nx,xlo,xhi,ny,ylo,yhi)

extern "C" {
extern void userin_();
extern Double_t func_(double x[]);
extern void spevnt_(Int_t *nret);
extern void exit(int);
JSFBases *bso;			// need to map xhfill to h1fill
void xhfill_(char *t, double *x, double *w, int len)
{
	char tmp[1024];
	int i;
	for (i=0; i<len; i++) tmp[i] = t[i];
	tmp[len] = '\0';
	bso->H1Fill(tmp,*x,*w);
}
void dhfill_(char *t, double *x, double *y, double *w, int len)
{
	char tmp[1024];
	int i;
	for (i=0; i<len; i++) tmp[i] = t[i];
	tmp[len] = '\0';
	bso->H2Fill(tmp,*x,*y,*w);
}
};

//_____________________________________________________________________________
TTSpring::TTSpring(const char *name, const char *title,
			 TTBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new TTSpringBuf("TTSpringBuf", 
	 "TTSpring event buffer", this);
  if( !bases ) { 
    TTBases *bs=new TTBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
TTSpring::~TTSpring()
{
  if( !fEventBuf ) { delete fEventBuf; fEventBuf = 0; }
}


//_____________________________________________________________________________
Bool_t TTSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void TTSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
TTBases::TTBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("TTBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("TTBases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("TTBases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("TTBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("TTBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("TTBases.NCALL","40000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 16;
  	fNWILD = 5;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 17;
  	fNWILD = 6;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 21;
  	fNWILD = 8;
  	fIOFF  = 0;
  } else {
  	printf("TTBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get ttbar specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"TTBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("TTBases.DeltaRoots","-0.595"),"%lg",&fDeltaRoots);
  sscanf(gJSF->Env()->GetValue("TTBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("TTBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("TTBases.WmModesLo","1"),"%d",&fWmModesLo);
  sscanf(gJSF->Env()->GetValue("TTBases.WmModesHi","12"),"%d",&fWmModesHi);
  sscanf(gJSF->Env()->GetValue("TTBases.WpModesLo","1"),"%d",&fWpModesLo);
  sscanf(gJSF->Env()->GetValue("TTBases.WpModesHi","12"),"%d",&fWpModesHi);
  sscanf(gJSF->Env()->GetValue("TTBases.VkmTB","1."),"%lg",&fVkmTB);
  sscanf(gJSF->Env()->GetValue("TTBases.BetaH","1."),"%lg",&fBetaH);
  sscanf(gJSF->Env()->GetValue("TTBases.NRQCD","1"),"%d",&fNRQCD);
  sscanf(gJSF->Env()->GetValue("TTBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("TTBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("TTBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("TTBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("TTBases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("TTBases.MassTop","170."),"%lg",&fMassTop);  

  Userin();
}


//_____________________________________________________________________________
void TTBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->TT generator\n");
  printf("  DeltaRoots  = %g (GeV)\n",fDeltaRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  VkmTB       = %g\n",fVkmTB);
  printf("  BetaH       = %g\n",fBetaH);
  printf("  Flag for NRQCD =%d\n",fNRQCD);
  printf("       = 0 ; Off\n");
  printf("       = 1 ; On\n");
  printf("  W- Decey Mode Lo =%d\n",fWmModesLo);
  printf("                Hi =%d\n",fWmModesHi);
  printf("  W+ Decey Mode Lo =%d\n",fWpModesLo);
  printf("                Hi =%d\n",fWpModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);


}

//_____________________________________________________________________________
Double_t TTBases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  TTBases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t TTBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void TTBases::Userin()
{
//
//   Initialize User parameters for Bases
//

  // Copy class data member into common /usmprm/
  usmprm_.alfi = fAlphai;
  usmprm_.alfs = fAlphas;
  usmprm_.amsw = fMassW;
  usmprm_.amsz = fMassZ;
  usmprm_.amsh = fMassHiggs;
  usmprm_.amst = fMassTop;

  // Copy class data member into common /usrprm/
  usrprm_.deltrs = fDeltaRoots;
  usrprm_.polebm = fPolElectron;
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;
  usrprm_.imd1lo = fWmModesLo;
  usrprm_.imd1hi = fWmModesHi;
  usrprm_.imd2lo = fWpModesLo;
  usrprm_.imd2hi = fWpModesHi;
  usrprm_.vkmt   = fVkmTB;
  usrprm_.beth   = fBetaH;
  usrprm_.nqcd   = fNRQCD;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

  Xhinit("h01", -1.0, 1.0, 50,"cos(theta_t) ");
  Xhinit("h02",  0.0,360., 50,"phi_t        ");
  Xhinit("h03", -1.0, 1.0, 50,"cos(theta_b_bar) ");
  Xhinit("h04",  0.0,360., 50,"phi_b_bar        ");
  Xhinit("h05", -1.0, 1.0, 50,"cos(theta_b) ");
  Xhinit("h06",  0.0,360., 50,"phi_b        ");
  Xhinit("h07",  0.0,  1., 50,"RS/ROOTS     ");
  Xhinit("h08", -10., 20., 50,"ROOTS-2*AMT  ");
  Xhinit("h09", 16900.,32400., 50,"S_1          ");
  Xhinit("h10", 16900.,32400., 50,"S_2          ");
  Xhinit("h11",  0.0, 50., 50,"|p|_cm       ");
  Xhinit("h12",  0.0, 50., 50,"|p|_lab      ");
  Xhinit("h13",  1.0, 13., 12,"W- decay mode ");
  Xhinit("h14",  1.0, 13., 12,"W+ decay mode ");
}

//_____________________________________________________________________________
void TTBases::Userout()
{
  Double_t fRoots = fMassTop*2 + fDeltaRoots;
  printf("End of TTBases\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iterations = %d\n",GetNoOfIterate());  
}

