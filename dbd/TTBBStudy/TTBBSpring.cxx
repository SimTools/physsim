//*(Update Record)
//*  2009/07/30  T.Tanabe	Initial version based on TTHStudy by K.Fujii

//////////////////////////////////////////////////////////////////
//
//  TTBBSpring
//  
//  e+e- -> ttbar gluon -> ttbar bbar
//  
//  In this program, meanings of integration variables are as follows.
//
//  Definition of vairables
//    Z( 1) : e- beam
//     ( 2) : e+ beam
//     ( 3) : ISR for e-
//     ( 4) : ISR for e+ ( or e-/e+ selection when ISR_LLR_Order=0)
//     ( 5) : e- helicity
//     ( 6) : helicity combination for final states.
//     ( 7) : m(t-bar)**2
//     ( 8) : m(t)**2
//     ( 9) : m(J)**2
//     (10) : m(t-t_bar)**2
//     (11) : m(W-)**2
//     (12) : m(W+)**2
//     (13) : cos(theta_g)
//     (14) : phi_g
//     (15) : cos(theta_t)     in t-t_bar rest frame
//     (16) : phi_t            in t-t_bar rest frame
//     (17) : cos(theta_q_bar) in t_bar rest frame
//     (18) : phi_q_bar        in t_bar rest frame
//     (19) : cos(theta_f)     in W- rest frame
//     (20) : phi_f            in W- rest frame
//     (21) : cos(theta_q)     in t rest frame
//     (22) : phi_q            in t rest frame
//     (23) : cos(theta_f_bar) in W+ rest frame
//     (24) : phi_f_bar        in W+ rest frame
//     (25) : cos(theta_f_bar) in Z rest frame
//     (26) : phi_f_bar        in Z rest frame
//     (27) : final state combination.
//     (28) : e- beam gaussian spread
//     (29) : e+ beam gaussian spread
//     (30) : random variable for Pt of ISR photon from e- beam
//     (31) : random variable for Pt of ISR photon from e+ beam
//     (32) : random variable for azimuthal angle of ISR photon from e- beam
//     (33) : random variable for azimuthal angle of ISR photon from e+ beam
//
//////////////////////////////////////////////////////////////////

#include "JSFModule.h"
#include "JSFLCFULL.h"
#include "JSFSpring.h"
#include "JSFSteer.h"
#include "JSFBases.h"
#include "JSFHadronizer.h"
#include "TTBBSpring.h"
#include "TSystem.h"

ClassImp(TTBBSpring)
ClassImp(TTBBSpringBuf)
ClassImp(TTBBBases)

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
TTBBSpring::TTBBSpring(const char *name, const char *title,
			 TTBBBases *bases)
  : JSFSpring(name, title, bases)
{
  fEventBuf = new TTBBSpringBuf("TTBBSpringBuf", 
	 "TTBBSpring event buffer", this);
  if( !bases ) { 
    TTBBBases *bs=new TTBBBases();
    SetBases(bs);
  }
}


//_____________________________________________________________________________
TTBBSpring::~TTBBSpring()
{
  if( !fEventBuf ) delete fEventBuf;
}


//_____________________________________________________________________________
Bool_t TTBBSpring::Initialize()
{
  // Make sure to set JSFHadronize::fCopySpringClassDataToBank = kFALSE

  JSFSpring::Initialize();

  JSFHadronizer *had=(JSFHadronizer*)gJSF->FindModule("JSFHadronizer");
  if(had) had->SetCopySpringClassDataToBank(kFALSE);

  if (fFile->IsWritable()) {
    TTBBBases *bs = (TTBBBases *)GetBases();
    TDirectory *last = gDirectory;
    fFile->cd("/conf");
    TList *dlist = gDirectory->GetListOfKeys();
    if (!dlist->FindObject("init")) gDirectory->mkdir("init");
    fFile->cd("/conf/init");
    bs->Write();
    last->cd(); 
    cerr << ">>>>>> TTBBBases written to file" << endl;
  }

  return kTRUE;
}


//_____________________________________________________________________________
void TTBBSpringBuf::Spevnt(Int_t &iret) { spevnt_(&iret); }

//_____________________________________________________________________________
TTBBBases::TTBBBases(const char *name, const char *title)
           : JSFBases(name, title)
{
//  Constructor of bases.  Default parameter should be initialized here
//

  bso = this;

// Get parameters from jsf.conf, if specified.

  sscanf(gJSF->Env()->GetValue("TTBBBases.ISRBM","3"),"%d",&fISRBM);
  sscanf(gJSF->Env()->GetValue("TTBBBases.WmModesLo","1"),"%d",&fWmModesLo);
  sscanf(gJSF->Env()->GetValue("TTBBBases.WmModesHi","12"),"%d",&fWmModesHi);
  sscanf(gJSF->Env()->GetValue("TTBBBases.WpModesLo","1"),"%d",&fWpModesLo);
  sscanf(gJSF->Env()->GetValue("TTBBBases.WpModesHi","12"),"%d",&fWpModesHi);
  sscanf(gJSF->Env()->GetValue("TTBBBases.ACC1","0.4"),"%lg",&fACC1);
  sscanf(gJSF->Env()->GetValue("TTBBBases.ACC2","0.2"),"%lg",&fACC2);
  sscanf(gJSF->Env()->GetValue("TTBBBases.ITMX1","5"),"%d",&fITMX1);
  sscanf(gJSF->Env()->GetValue("TTBBBases.ITMX2","5"),"%d",&fITMX2);
  sscanf(gJSF->Env()->GetValue("TTBBBases.NCALL","80000"),"%d",&fNCALL);
  sscanf(gJSF->Env()->GetValue("TTBBBases.AltColorSinglets","0"),"%d",&fAltColorSinglets);

#ifdef WITH_DBD_STANDARD
 sscanf(gJSF->Env()->GetValue("TTBBBases.ISR_LLA_Order","0"),"%d",&fISR_LLA_Order);
// if fISR_LLA_Order is gt 0, LLA form with order, LLA_Order, is used 
// to generate ISR spectrum
  sscanf(gJSF->Env()->GetValue("TTBBBases.Store_Remnants","1"),"%d",&fStoreRemnants);
// (0, 1)=(not save, save) ISR photon info.
  sscanf(gJSF->Env()->GetValue("TTBBBases.Store_Beams","0"),"%d",&fStoreBeams);
// (0,1)=(not save, save) initial e+/e- beam after beamstrahlung
  if( fISRBM < 100 &&  ( fStoreRemnants != 0 || fISR_LLA_Order != 0 )) {
    std::cout << "Fatal error in TTBBBases .. " << std::endl;
    std::cout << "  Input parameter .. fISRBM=" << fISRBM << " fStoreRemnants=" << fStoreRemnants
              << " fISR_LLA_Order=" << fISR_LLA_Order << std::endl;
    std::cout << "  But fISRBM should be > 100 if fStoreRemnants or fISR_LLA_Order .ne. 0 " << std::endl;
    exit(-1);
  }
  fLumiFileDirectory=gJSF->Env()->GetValue("TTBBBases.LumiFileDirectory","");
#endif

  Int_t fIOFF;
  if ( fISRBM == 1 ) {
  	fNDIM  = 23;
  	fNWILD = 7;
  	fIOFF  = 3;
#ifdef WITH_DBD_STANDARD
        fNDIM  = 23;
        fNWILD = 8;
        fIOFF  = 4;
#endif
  } else if ( fISRBM == 2 ) {
   	fNDIM  = 24;
  	fNWILD = 8;
  	fIOFF  = 2;
#ifdef WITH_DBD_STANDARD
        fNDIM  = 25;
        fNWILD = 8;
        fIOFF  = 2;
#endif
  } else if ( fISRBM == 3 ) {
  	fNDIM  = 28;
  	fNWILD = 11;
  	fIOFF  = 0;
#ifdef WITH_DBD_STANDARD
  } else if ( fISRBM >100 ) {
        fNDIM  = 33;
        fNWILD = 12;
        fIOFF  = 0;
#endif
  } else {
  	printf("TTBBBases: Invalid ISRBM = %d\n",fISRBM);
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
    sprintf(pname,"TTBBBases.X%i2.2Range",i+1);
    sscanf(gJSF->Env()->GetValue(pname,"0.0 1.0"),"%lg%lg",&fXL[i],&fXU[i]);
    if ( i >= fNWILD ) { 
       fIG[i] = 0;
    } else {
       fIG[i] = 1;
    }
    fISHUFL[i] = i + fIOFF + 1;
  }

  sscanf(gJSF->Env()->GetValue("TTBBBases.Roots","700."),"%lg",&fRoots);
  sscanf(gJSF->Env()->GetValue("TTBBBases.PolElectron","0."),"%lg",&fPolElectron);
  sscanf(gJSF->Env()->GetValue("TTBBBases.PolPositron","0."),"%lg",&fPolPositron);
  sscanf(gJSF->Env()->GetValue("TTBBBases.SigmaEbeam","0.005"),"%lg",&fSigmaEbeam);
  sscanf(gJSF->Env()->GetValue("TTBBBases.Alphai","128."),"%lg",&fAlphai);
  sscanf(gJSF->Env()->GetValue("TTBBBases.Alphas","0.120"),"%lg",&fAlphas);
  sscanf(gJSF->Env()->GetValue("TTBBBases.MassW","80.0"),"%lg",&fMassW);
  sscanf(gJSF->Env()->GetValue("TTBBBases.MassZ","91.18"),"%lg",&fMassZ);
  sscanf(gJSF->Env()->GetValue("TTBBBases.MassHiggs","120."),"%lg",&fMassHiggs);
  sscanf(gJSF->Env()->GetValue("TTBBBases.MassTop","170."),"%lg",&fMassTop);  

#ifdef WITH_DBD_STANDARD
   fFinalStatesMix = gJSF->Env()->GetValue("TTBBBases.FinalStatesMix",0);

  if( fISRBM > 100 ) {
    if( fLumiFileDirectory.size() < 1 ) { 
      std::cout << "Error!! TTBBBases.LumiFile_Directory is not set, though ISRBM> 100" << std::endl;    
      std::cout << "  ISRBM=" << fISRBM << std::endl;
      exit(-1);
    }
    gSystem->Setenv("LUMI_LINKER",(fLumiFileDirectory+std::string("/lumi_linker_000")).c_str()); 
    gSystem->Setenv("PHOTONS_B1",(fLumiFileDirectory+std::string("/photons_beam1_linker_000")).c_str()); 
    gSystem->Setenv("PHOTONS_B2",(fLumiFileDirectory+std::string("/photons_beam2_linker_000")).c_str()); 
    gSystem->Setenv("EBEAM",(fLumiFileDirectory+std::string("/ebeam_in_linker_000")).c_str()); 
    gSystem->Setenv("PBEAM",(fLumiFileDirectory+std::string("/pbeam_in_linker_000")).c_str()); 
  }
#endif

}


//_____________________________________________________________________________
void TTBBBases::PrintParameters()
{
//  Print parameters
//
  
  printf("Parameters for ee->TTBB generator\n");
  printf("  Roots  = %g (GeV)\n",fRoots);
  printf("  PolElectron = %g\n",fPolElectron);
  printf("  PolPositron = %g\n",fPolPositron);
  printf("  SigmaEbeam  = %g\n",fSigmaEbeam);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("       >100 ; DBD standard lumi spectrum \n");
  printf("  W- Decey Mode Lo =%d\n",fWmModesLo);
  printf("                Hi =%d\n",fWmModesHi);
  printf("  W+ Decey Mode Lo =%d\n",fWpModesLo);
  printf("                Hi =%d\n",fWpModesHi);
  printf("    W decay mode: (1,2,3)=(n1e1, n2e2, n3e3) (4,5,6)=(ud, cd, td)\n");
  printf("                  (7,8,9)=(us, cs, ts)    (10,11,12)=(ub, cb, tb)\n");
  printf("  Alternative Color Singlets =%d\n",fAltColorSinglets);

#ifdef WITH_DBD_STANDARD
  printf("  Final state mode mix    = %d\n",fFinalStatesMix);
  printf("     If not 0, W- decay mode setting applies to one W and W+ \n");
  printf("     decay mode applies to another W.  Set non 0 value to generate \n");
  printf("     W-W+ --> l-nu qq and l+nu qq simultaneously.\n");
  printf("  LLA Order for ISR =%d\n",fISR_LLA_Order);
  printf("  Store ISR remnants (0,1)=(no, yes): %d\n",fStoreRemnants);
  printf("  Store e+e- beams after beamstrahlung (0,1)=(no,yes) : %d\n",fStoreBeams);
  if ( fISRBM > 100 ) {
    printf("  Lumi File Directory : %s\n",fLumiFileDirectory.c_str());
  }
#endif

  printf("  Bases integration parameters..\n");
  printf("  ITMX1=%d  ITMX2=%d  NCALL=%d\n",fITMX1, fITMX2, fNCALL);
  printf("  ACC1 =%g  ACC2 =%g\n",fACC1,fACC2);

  printf(" ========= End of TTBBBases parameters =============================================\n");
}

//_____________________________________________________________________________
Double_t TTBBBases::Func()		// new style not yet implemented
{
  cerr << ":::::: ERROR "
       << "  TTBBBases::Func() not implemented !!!"
       << "  Will STOP immediately." << endl;
       exit(1);
  return 0.;
}

//_____________________________________________________________________________
Double_t TTBBBases::Func(Double_t x[])
{
//  Bases Integrand
//
  double val=func_(x);
  return val;

}

//_____________________________________________________________________________
void TTBBBases::Userin()
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
  usrprm_.polpbm = fPolPositron;
  usrprm_.sgmebm = fSigmaEbeam;
  usrprm_.isrb   = fISRBM;
  usrprm_.imd1lo = fWmModesLo;
  usrprm_.imd1hi = fWmModesHi;
  usrprm_.imd2lo = fWpModesLo;
  usrprm_.imd2hi = fWpModesHi;
  usrprm_.altcs  = fAltColorSinglets;

#ifdef WITH_DBD_STANDARD
  usrprm_.ndcysl = fFinalStatesMix;
  usrprm_.isr_lla_order = fISR_LLA_Order;
  usrprm_.ns_remnants  = fStoreRemnants;
  usrprm_.ns_beams   = fStoreBeams;
#endif

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
  Xhinit("h03", 300.,500., 50,"m(t-tb)      ");
  Xhinit("h04", -1.0, 1.0, 50,"cos(theta_t) ");
  Xhinit("h05",  0.0,360., 50,"phi_t        ");
  Xhinit("h06", -1.0, 1.0, 50,"cos(theta_t)_lab    ");
  Xhinit("h07",  1.0,  5.,  4,"helicity combination");
  Xhinit("h08", 90.0,200., 50,"m_tbar       ");
  Xhinit("h09", 90.0,200., 50,"m_t          ");
  Xhinit("h10",  0.0,150., 50,"m_g          ");
  Xhinit("h11", 40.0,120., 50,"m_W-         ");
  Xhinit("h12", 40.0,120., 50,"m_W+         ");
  Xhinit("h13",  0.0, 1.0, 50,"rsh/roots    ");
  Xhinit("h15",  1.0, 13., 12,"W- decay mode ");
  Xhinit("h16",  1.0, 13., 12,"W+ decay mode ");
  Dhinit("hd21",0.,1.,50,-1.,1.,50,"E_H/E_bm-cos(th_H)");
  Dhinit("hd22",0.,1.,50, 0.,1.,50,"E_t/E_bm-E_tb/E_bm");
  Dhinit("dalitz",0.,90000.,50, 0.,90000.,50,"m^2(tg) vs. m^2(t~g)");
}

//_____________________________________________________________________________
void TTBBBases::Userout()
{
  printf("End of TTBBBases\n");
  printf("ISRBM = %d\n",fISRBM);
  printf("  Flag for ISR/BM Effects(ISRBM) =%d\n",fISRBM);
  printf("       = 1 ; None\n");
  printf("       = 2 ; ISR only\n");
  printf("       = 3 ; ISR + BM\n");
  printf("  W- Decey Mode Lo =%d\n",fWmModesLo);
  printf("                Hi =%d\n",fWmModesHi);
  printf("  W+ Decey Mode Lo =%d\n",fWpModesLo);
  printf("                Hi =%d\n",fWpModesHi);
  printf("Ecm                  = %g (GeV)\n",fRoots);
  printf("Total Cross section  = %g +- %g (fb)\n",GetEstimate(),GetError());
  printf("Number of iterations = %d\n",GetNoOfIterate());  
}

// __________________________________________________________________________
void TTBBSpring::GetPy6frmProb(int nseq, double prob[7])
{
  prob[0]=1.0;  //p12
  prob[1]=0.0;  //p13
  prob[2]=0.0;  //p21
  prob[3]=0.0;  //p23
  prob[4]=0.0;  //p31
  prob[5]=0.0;  //p32
  prob[6]=1.0;  //ptop
  return ;
}

