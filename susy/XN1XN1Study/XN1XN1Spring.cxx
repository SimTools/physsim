//////////////////////////////////////////////////////////////////
//
//  XN1XN1Spring
//  
//  e+e- -> XN1XN1
//  
//  In this program, meanings of integration variables are as follows.
//
//  Definition of vairables
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : bremsstrahlng
//     ( 4) : XN1 decay mode
//     ( 5) : e- helicity
//     ( 6) : m(X0i)**2
//     ( 7) : m(X0j)**2
//     ( 8) : cos(theta_X0i)
//     ( 9) : phi_X0i
//     (10) : cos(theta_a)     in X0i rest frame
//     (11) : phi_a            in X0i rest frame
//     (12) : cos(theta_a)     in X0j rest frame
//     (13) : phi_a            in X0j rest frame
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
#include "XN1XN1Spring.h"

ClassImp(XN1XN1Spring)
ClassImp(XN1XN1SpringBuf)
ClassImp(XN1XN1Bases)

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
XN1XN1Spring::XN1XN1Spring(const char *name, const char *title,
			 XN1XN1Bases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new XN1XN1SpringBuf("XN1XN1SpringBuf", 
	 "XN1XN1Spring event buffer", this);
  if( !bases ) { 
    XN1XN1Bases *bs=new XN1XN1Bases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
XN1XN1Spring::~XN1XN1Spring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t XN1XN1Spring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  if (fFile->IsWritable()) {
    XN1XN1Bases *bs = static_cast<XN1XN1Bases *>(GetBases());
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    cerr << ">>>>>> Ecm = " << bs->GetRoots() << " [GeV] " << endl;
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> XN1XN1Bases written to file" << endl;
  }

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);
  return kTRUE;
}


//_____________________________________________________________________________
void XN1XN1SpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________

XN1XN1Bases::XN1XN1Bases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;				// set this for xhfill_

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.NCALL","80000"),"%d",&fNCALL);

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 10;
  	fNWILD = 5;
  	fIOFF  = 3;
        for(Int_t i=0;i<fNDIM;i++){
           fISHUFL[i] = i + fIOFF + 1;
        }
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 11;
  	fNWILD = 6;
  	fIOFF  = 2;
        fISHUFL[0] = 4;
        fISHUFL[1] = 5;
        fISHUFL[2] = 6;
        fISHUFL[3] = 7;
        fISHUFL[4] = 8;
        fISHUFL[5] = 3;
        for(Int_t i=6;i<fNDIM;i++){
           fISHUFL[i] = i + fIOFF + 1;
        }
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 15;
  	fNWILD = 6;
  	fIOFF  = 0;
        fISHUFL[0] = 4;
        fISHUFL[1] = 5;
        fISHUFL[2] = 6;
        fISHUFL[3] = 7;
        fISHUFL[4] = 8;
        fISHUFL[5] = 3;
        fISHUFL[6] = 1;
        fISHUFL[7] = 2;
        for(Int_t i=8;i<fNDIM;i++){
           fISHUFL[i] = i + fIOFF + 1;
        }
  } else {
  	printf("XN1XN1Bases: Invalid ISRBM = %d\n",fISRBM);
  	printf("    Will STOP immediately\n");
  	exit(1);
  }

//
// Get XN1XN1bar specific parameters.
//
  Char_t pname[40];
  for(Int_t i=0;i<fNDIM;i++){
    sprintf(pname,"XN1XN1Bases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
  }

  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.Roots","500."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.X1ModesLo","1"),"%d",&fX1ModesLo);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.X1ModesHi","2"),"%d",&fX1ModesHi);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.X2ModesLo","1"),"%d",&fX2ModesLo);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.X2ModesHi","2"),"%d",&fX2ModesHi);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.MassHiggs","9999."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.MassTop","170."),"%lg",&fMassTop);  
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.m0","70."),"%lg",&fm0);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.mu","400."),"%lg",&fmu);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.GUT","1"),"%d",&fGUT);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.M2","250."),"%lg",&fM2);
  if (fGUT != 1) {
     sscanf(gJSF->Env()->GetValue("XN1XN1Bases.M1","125."),"%lg",&fM1);
  }
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.tanb","+2."),"%lg",&ftanb);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.mA","-9999."),"%lg",&fmA);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.Atau","0."),"%lg",&fAsoft[0]);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.At","0."),"%lg",&fAsoft[1]);
  sscanf(gJSF->Env()->GetValue("XN1XN1Bases.Ab","0."),"%lg",&fAsoft[2]);
  if (fmA == -9999.) fmA = TMath::Sqrt(fm0*fm0+fmu*fmu);
}


//_____________________________________________________________________________
void XN1XN1Bases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->XN1XN1 generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  // printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  XN width = %g\n",fWidthChin2);
  printf("  1st XN Decey Mode Lo =%d\n",fX1ModesLo);
  printf("                    Hi =%d\n",fX1ModesHi);
  printf("  2nd XN Decey Mode Lo =%d\n",fX2ModesLo);
  printf("                    Hi =%d\n",fX2ModesHi);

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);
}
//_____________________________________________________________________________
Double_t XN1XN1Bases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  XN1XN1Bases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t XN1XN1Bases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void XN1XN1Bases::Userin()
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
  usrprm_.imd1lo = fX1ModesLo;
  usrprm_.imd1hi = fX1ModesHi;
  usrprm_.imd2lo = fX2ModesLo;
  usrprm_.imd2hi = fX2ModesHi;
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
  Xhinit("h05",  0.0, qmx, 50,"M_X0i                ");
  Xhinit("h06", -1.0, 1.0, 50,"cos(theta_a) X0j     ");
  Xhinit("h07",  0.0,360., 50,"phi_b                ");
  Xhinit("h08",  0.0, qmx, 50,"M_X0j                ");
  Xhinit("h09",  0.0,  1., 50,"RS/ROOTS             ");
  Xhinit("h10",  1.0,  9.,  8,"Hel. comb.           ");
  Xhinit("h11",  1.0,  5.,  4,"Decay mode           ");
}

//_____________________________________________________________________________
void XN1XN1Bases::Userout()
{ 
  printf("End of Bases of ee --> XN1 XN1 process\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iteration  = %d\n",GetNoOfIterate());  
}
