//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  TTHSpring
//  
//  e+e- -> ttbar H
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
#include "TTHSpring.h"

ClassImp(TTHSpring)
ClassImp(TTHSpringBuf)
ClassImp(TTHBases)

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
TTHSpring::TTHSpring(const char *name, const char *title,
			 TTHBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new TTHSpringBuf("TTHSpringBuf", 
	 "TTHSpring event buffer", this);
  if( !bases ) { 
    TTHBases *bs=new TTHBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
TTHSpring::~TTHSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t TTHSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);

  if (fFile->IsWritable()) {
    TTHBases *bs = (TTHBases *)GetBases();
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> TTHBases written to file" << endl;
  }

  return kTRUE;
}


//_____________________________________________________________________________
void TTHSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
TTHBases::TTHBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("TTHBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("TTHBases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("TTHBases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("TTHBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("TTHBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("TTHBases.NCALL","80000"),"%d",&fNCALL);

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
  	printf("TTHBases: Invalid ISRBM = %d\n",fISRBM);
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
    sprintf(pname,"TTHBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("TTHBases.Roots","700."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("TTHBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("TTHBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("TTHBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("TTHBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("TTHBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("TTHBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("TTHBases.MassHiggs","120."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("TTHBases.MassTop","170."),"%lg",&fMassTop);  
}


//_____________________________________________________________________________
void TTHBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->TTH generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);


}

//_____________________________________________________________________________
Double_t TTHBases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  TTHBases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t TTHBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void TTHBases::Userin()
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

  // Copy class data member into common /bshufl/
  bshufl_.nzz = fNDIM;
  for (Int_t i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms

  Xhinit("h01", -1.0, 1.0, 50,"cos(theta_H) ");
  Xhinit("h02",  0.0,360., 50,"phi_H        ");
  Xhinit("h03",  0.0, 1.0, 50,"m(t-tb)/roots");
  Xhinit("h04", -1.0, 1.0, 50,"cos(theta_t) ");
  Xhinit("h05",  0.0,360., 50,"phi_t        ");
  Xhinit("h06", -1.0, 1.0, 50,"cos(theta_t)_lab    ");
  Xhinit("h07",  1.0,  5.,  4,"helicity combination");
  Xhinit("h08", 90.0,200., 50,"m_tbar       ");
  Xhinit("h09", 90.0,200., 50,"m_t          ");
  Xhinit("h10", 50.0,150., 50,"m_H          ");
  Xhinit("h11", 40.0,120., 50,"m_W-         ");
  Xhinit("h12", 40.0,120., 50,"m_W+         ");
  Dhinit("hd21",0.,1.,50,-1.,1.,50,"E_H/E_bm-cos(th_H)");
  Dhinit("hd22",0.,1.,50, 0.,1.,50,"E_t/E_bm-E_tb/E_bm");
}

//_____________________________________________________________________________
void TTHBases::Userout()
{
  printf("End of TTHBases\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iterations = %d\n",GetNoOfIterate());  
}





















