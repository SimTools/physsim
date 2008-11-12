//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  TTZSpring
//  
//  e+e- -> ttbar Z
//  
//  In this program, meanings of integration variables are as follows.
//
//  Definition of vairables
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : ISR
//     ( 4) : e- helicity
//     ( 5) : helicity combination for final states.
//     ( 6) : m(t-bar)**2
//     ( 7) : m(t)**2
//     ( 8) : m(J)**2
//     ( 9) : m(t-t_bar)**2
//     (10) : m(W-)**2
//     (11) : m(W+)**2
//     (12) : cos(theta_H)
//     (13) : phi_H
//     (14) : cos(theta_t)     in t-t_bar rest frame
//     (15) : phi_t            in t-t_bar rest frame
//     (16) : cos(theta_q_bar) in t_bar rest frame
//     (17) : phi_q_bar        in t_bar rest frame
//     (18) : cos(theta_f)     in W- rest frame
//     (19) : phi_f            in W- rest frame
//     (20) : cos(theta_q)     in t rest frame
//     (21) : phi_q            in t rest frame
//     (22) : cos(theta_f_bar) in W+ rest frame
//     (23) : phi_f_bar        in W+ rest frame
//     (24) : cos(theta_f_bar) in Z rest frame
//     (25) : phi_f_bar        in Z rest frame
//     (26) : final state combination.
//     (27) : e- beam gaussian spread
//     (28) : e+ beam gaussian spread
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "TTZSpring.h"

ClassImp(TTZSpring)
ClassImp(TTZSpringBuf)
ClassImp(TTZBases)

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
TTZSpring::TTZSpring(const char *name, const char *title,
			 TTZBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new TTZSpringBuf("TTZSpringBuf", 
	 "TTZSpring event buffer", this);
  if( !bases ) { 
    TTZBases *bs=new TTZBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
TTZSpring::~TTZSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t TTZSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);

  if (fFile->IsWritable()) {
    TTZBases *bs = (TTZBases *)GetBases();
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd();
    cerr << ">>>>>> TTZBases written to file" << endl;
  }

  return kTRUE;
}


//_____________________________________________________________________________
void TTZSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
TTZBases::TTZBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("TTZBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("TTZBases.ZModesLo","1"),"%d",&fZModesLo);
  sscanf(gJSF->Env()->GetValue("TTZBases.ZModesHi","12"),"%d",&fZModesHi);
  sscanf(gJSF->Env()->GetValue("TTZBases.WmModesLo","1"),"%d",&fWmModesLo);
  sscanf(gJSF->Env()->GetValue("TTZBases.WmModesHi","12"),"%d",&fWmModesHi);
  sscanf(gJSF->Env()->GetValue("TTZBases.WpModesLo","1"),"%d",&fWpModesLo);
  sscanf(gJSF->Env()->GetValue("TTZBases.WpModesHi","12"),"%d",&fWpModesHi);
  sscanf(gJSF->Env()->GetValue("TTZBases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("TTZBases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("TTZBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("TTZBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("TTZBases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 23;
  	fNWILD = 7;
  	fIOFF  = 3;
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 24;
  	fNWILD = 8;
  	fIOFF  = 2;
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 28;
  	fNWILD = 11;
  	fIOFF  = 0;
  } else {
  	printf("TTZBases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }
#ifdef ZEROWIDTH
        fNWILD -= 6;
#endif

//
// Get ttbar specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"TTZBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("TTZBases.Roots","700."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("TTZBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("TTZBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("TTZBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("TTZBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("TTZBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("TTZBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("TTZBases.MassHiggs","120."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("TTZBases.MassTop","170."),"%lg",&fMassTop);  
}


//_____________________________________________________________________________
void TTZBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->TTZ generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  Z- Decey Mode Lo =%d\n",fZModesLo);
  printf("                Hi =%d\n",fZModesHi);
  printf("  W- Decey Mode Lo =%d\n",fWmModesLo);
  printf("                Hi =%d\n",fWmModesHi);
  printf("  W+ Decey Mode Lo =%d\n",fWpModesLo);
  printf("                Hi =%d\n",fWpModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);


}

//_____________________________________________________________________________
Double_t TTZBases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  TTZBases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t TTZBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void TTZBases::Userin()
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
  usrprm_.sqrts  = fRoots;
  usrprm_.polebm = fPolElectron;
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;
  usrprm_.imdzlo = fZModesLo;
  usrprm_.imdzhi = fZModesHi;
  usrprm_.imd1lo = fWmModesLo;
  usrprm_.imd1hi = fWmModesHi;
  usrprm_.imd2lo = fWpModesLo;
  usrprm_.imd2hi = fWpModesHi;

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

  Xhinit("h01", -1.0, 1.0, 50,"cos(theta_Z) ");
  Xhinit("h02",  0.0,360., 50,"phi_Z        ");
  Xhinit("h03",  0.0, 1.0, 50,"m(t-tb)/roots");
  Xhinit("h04", -1.0, 1.0, 50,"cos(theta_t) ");
  Xhinit("h05",  0.0,360., 50,"phi_t        ");
  Xhinit("h06", -1.0, 1.0, 50,"cos(theta_t)_lab    ");
  Xhinit("h07",  1.0,  5.,  4,"helicity combination");
  Dhinit("hd09",0.,1.,50,-1.,1.,50,"E_Z/E_bm-cos(th_Z)");
  Dhinit("hd10",0.,1.,50, 0.,1.,50,"E_t/E_bm-E_tb/E_bm");
  Xhinit("h12",  1.0, 13., 12,"Z  decay mode ");
  Xhinit("h13",  1.0, 13., 12,"W- decay mode ");
  Xhinit("h14",  1.0, 13., 12,"W+ decay mode ");
}

//_____________________________________________________________________________
void TTZBases::Userout()
{
  printf("End of TTZBases\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iterations = %d\n",GetNoOfIterate());  
}

