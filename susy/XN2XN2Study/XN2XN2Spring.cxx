//*LastUpdate:  v.01.01 19-December-1998 Keisuke Fujii
//*-- Author :  	19-December-1998 Keisuke Fujii

//////////////////////////////////////////////////////////////////
//
//  XN2XN2Spring
//  
//  e+e- -> XN2XN2
//  
//  In this program, meanings of integration variables are as follows.
//
//  Definition of vairables
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : bremsstrahlng
//     ( 4) : Z2 decay mode
//     ( 5) : e- helicity
//     ( 6) : m(ab)**2 for X0i
//     ( 7) : m(ac)**2
//     ( 8) : m(ab)**2 for X0j
//     ( 9) : m(ac)**2
//     (10) : m(X0i)**2
//     (11) : m(X0j)**2
//     (12) : cos(theta_X0i)
//     (13) : phi_X0i
//     (14) : cos(theta_a)     in X0i rest frame
//     (15) : phi_a            in X0i rest frame
//     (16) : phi_b            in X0i rest frame
//     (17) : cos(theta_a)     in X0j rest frame
//     (18) : phi_a            in X0j rest frame
//     (19) : phi_b            in X0j rest frame
//     (20) : E_beam spread of E-
//     (21) : E_beam spread of E+
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "XN2XN2Spring.h"

ClassImp(XN2XN2Spring)
ClassImp(XN2XN2SpringBuf)
ClassImp(XN2XN2Bases)

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
};

//_____________________________________________________________________________
XN2XN2Spring::XN2XN2Spring(const char *name, const char *title,
			 XN2XN2Bases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new XN2XN2SpringBuf("XN2XN2SpringBuf", 
	 "XN2XN2Spring event buffer", this);
  if( !bases ) { 
    XN2XN2Bases *bs=new XN2XN2Bases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
XN2XN2Spring::~XN2XN2Spring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t XN2XN2Spring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void XN2XN2SpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________

XN2XN2Bases::XN2XN2Bases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;				// set this for xhfill_

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 16;
  	fNWILD = 6;
  	fIOFF  = 3;
        for(Int_t i=0;i<fNDIM;i++){
           fISHUFL[i] = i + fIOFF + 1;
        }
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 17;
  	fNWILD = 7;
  	fIOFF  = 2;
        fISHUFL[0] = 4;
        fISHUFL[1] = 5;
        fISHUFL[2] = 6;
        fISHUFL[3] = 7;
        fISHUFL[4] = 8;
        fISHUFL[5] = 9;
        fISHUFL[6] = 3;
        for(Int_t i=7;i<fNDIM;i++){
           fISHUFL[i] = i + fIOFF + 1;
        }
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 21;
  	fNWILD = 7;
  	fIOFF  = 0;
        fISHUFL[0] = 4;
        fISHUFL[1] = 5;
        fISHUFL[2] = 6;
        fISHUFL[3] = 7;
        fISHUFL[4] = 8;
        fISHUFL[5] = 9;
        fISHUFL[6] = 3;
        fISHUFL[7] = 1;
        fISHUFL[8] = 2;
        for(Int_t i=9;i<fNDIM;i++){
           fISHUFL[i] = i + fIOFF + 1;
        }
  } else {
  	printf("XN2XN2Bases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get XN2XN2bar specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"XN2XN2Bases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
  }

  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.Z1ModesLo","1"),"%d",&fZ1ModesLo);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.Z1ModesHi","12"),"%d",&fZ1ModesHi);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.Z2ModesLo","1"),"%d",&fZ2ModesLo);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.Z2ModesHi","12"),"%d",&fZ2ModesHi);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.MassTop","170."),"%lg",&fMassTop);  
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.m0","70."),"%lg",&fm0);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.mu","400."),"%lg",&fmu);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.GUT","1"),"%d",&fGUT);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.M2","250."),"%lg",&fM2);
  if (fGUT != 1) {
     sscanf(gJSF->Env()->GetValue("XN2XN2Bases.M1","125."),"%lg",&fM1);
  }
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.tanb","+2."),"%lg",&ftanb);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.mA","-9999."),"%lg",&fmA);
  sscanf(gJSF->Env()->GetValue("SFSFBases.Atau","0."),"%lg",&fAsoft[0]);
  sscanf(gJSF->Env()->GetValue("SFSFBases.At","0."),"%lg",&fAsoft[1]);
  sscanf(gJSF->Env()->GetValue("SFSFBases.Ab","0."),"%lg",&fAsoft[2]);
  sscanf(gJSF->Env()->GetValue("XN2XN2Bases.WidthChin2","-1."),"%lg",&fWidthChin2);
  if (fmA == -9999.) fmA = TMath::Sqrt(fm0*fm0+fmu*fmu);
}


//_____________________________________________________________________________
void XN2XN2Bases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->XN2XN2 generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  XN_2 width = %g\n",fWidthChin2);
  printf("  Z1 Decey Mode Lo =%d\n",fZ1ModesLo);
  printf("                Hi =%d\n",fZ1ModesHi);
  printf("  Z2 Decey Mode Lo =%d\n",fZ2ModesLo);
  printf("                Hi =%d\n",fZ2ModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);
}
//_____________________________________________________________________________
Double_t XN2XN2Bases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  XN2XN2Bases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t XN2XN2Bases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void XN2XN2Bases::Userin()
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

  // Copy class data member into common /usmprm/
  ussprmp_.am0  = fm0;
  ussprmp_.amu  = fmu;
  ussprmp_.am2  = fM2;
  ussprmp_.isgut = fGUT;
  if (fGUT != 1) {
    ussprmp_.am1  = fM1;
    Double_t s2w = (1-fMassW/fMassZ)*(1+fMassW/fMassZ);
    ussprmp_.am3  = fAlphas*s2w*fAlphai*fM2;
  }
  ussprmp_.tanb = ftanb;
  ussprmp_.ama  = fmA;
  ussprmp_.asft[0] = fAsoft[0];
  ussprmp_.asft[1] = fAsoft[1];
  ussprmp_.asft[2] = fAsoft[2];

  // Copy class data member into common /usrprm/
  usrprm_.sqrts  = fRoots;
  usrprm_.polebm = fPolElectron;
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;
  usrprm_.imd1lo = fZ1ModesLo;
  usrprm_.imd1hi = fZ1ModesHi;
  usrprm_.imd2lo = fZ2ModesLo;
  usrprm_.imd2hi = fZ2ModesHi;
  usrprm_.gamsz2 = fWidthChin2;

  // Copy class data member into common /bshufl/
  Int_t i;
  bshufl_.nzz = fNDIM;
  for (i=0; i<fNDIM; i++) bshufl_.ishufl[i] = fISHUFL[i];
  
  // Initialize physical constants, etc.
  userin_();

  // Print parameters.
  PrintParameters();

  // Define histograms
  Double_t ptmx = fRoots/2.;
  Double_t qmx  = fRoots/2;
  Xhinit("h01", -1.0, 1.0, 50,"cos(theta_X0i)       ");
  Xhinit("h02",  0.0,360., 50,"phi_X0i              ");
  Xhinit("h03", -1.0, 1.0, 50,"cos(theta_a) X0i     ");
  Xhinit("h04",  0.0,360., 50,"phi_b                ");
  Xhinit("h05",  0.0, qmx, 50,"M_ab                 ");
  Xhinit("h06",  0.0, qmx, 50,"M_ac                 ");
  Xhinit("h07", -1.0, 1.0, 50,"cos(theta_a) X0j     ");
  Xhinit("h08",  0.0,360., 50,"phi_b                ");
  Xhinit("h09",  0.0, qmx, 50,"M_ab                 ");
  Xhinit("h10",  0.0, qmx, 50,"M_ac                 ");
  Xhinit("h11",  0.0,ptmx, 50,"Missing PT           ");
  Xhinit("h12",  1.0, 33., 32,"Hel. comb.           ");
  Xhinit("h13",  0.0, qmx, 50,"E_f    X0i           ");
  Xhinit("h14",  0.0, qmx, 50,"E_fbar X0i           ");
  Xhinit("h15",  0.0, qmx, 50,"E_f    X0j           ");
  Xhinit("h16",  0.0, qmx, 50,"E_fbar X0j           ");
  Xhinit("h17", -1.0, 1.0, 50,"cos(theta_f)    X0i  ");
  Xhinit("h18", -1.0, 1.0, 50,"cos(theta_fbar) X0i  ");
  Xhinit("h19", -1.0, 1.0, 50,"cos(theta_f)    X0j  ");
  Xhinit("h20", -1.0, 1.0, 50,"cos(theta_fbar) X0j  ");
  Xhinit("h21",   .0,  1., 50,"sqrt(rs)/sqrt(root)  ");
  Xhinit("h22",  0.0, qmx, 50,"M_X0i                ");
  Xhinit("h23",  0.0, qmx, 50,"M_X0j                ");
  Xhinit("h24",  1.0, 25., 24,"Decay mode           ");
}

//_____________________________________________________________________________
void XN2XN2Bases::Userout()
{ 
  printf("End of Bases of ee --> XN2 XN2 process\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iteration  = %d\n",GetNoOfIterate());  
}
