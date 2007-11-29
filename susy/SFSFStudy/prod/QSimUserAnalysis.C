//*************************************************************************
//* ================
//*  UserAnalysis.C
//* ================
//*
//* (Description)
//*    A very primitive sample script to study e+e- -> SFSF.
//* (Usage)
//*    	 $ jsf gui.C
//*    gui.C then invokes this script from within before event loop
//*    starts.
//* (Update Recored)
//*    1999/05/26  K.Fujii	Original version.
//*    1999/08/10  K.Fujii	Revised to use Anlib classes.
//*
//*************************************************************************
//
#include <iostream.h>
//*------------------------*//
//* User Analysis          *//
//*------------------------*//

Bool_t gDEBUG = kFALSE;

// Canvas //

TCanvas    *cHist;

#define LEPTONID 13

// Hists //

TH1F *hStat;
TH1F *hNleptons;
TH1F *hElepton;
TH1F *hEvis;
TH1F *hPt;
TH1F *hCoslm;
TH1F *hCoslp;
TH1F *hMll;
TH1F *hAcop;

// Event Information //

Double_t fEcm=350.0;	// CM energy (GeV)
Int_t    fNleptons;	// Lepton multiplicity
Double_t fElepm;	// l^- energy
Double_t fElepp;	// l^+ energy
Double_t fEvis;		// Visible energy
Double_t fPt;		// Pt
Double_t fCoslm;	// cos(theta_lepton^-)
Double_t fCoslp;	// cos(theta_lepton^+)
Double_t fMll;		// m(l^+l^-)
Double_t fAcop;		// Acoplanarity

// Event selection condition //

static const Int_t maxcut = 50;
Char_t cutName[maxcut][100];	// Cut names

Int_t    xNleptons =  2   ; 	// Number of leptons
Double_t xElepton  =  2.00; 	// Lepton E minimum
Double_t xEvis     = 10.00; 	// Minimum visible energy
Double_t xPt       = 10.00; 	// Pt minimum
Double_t xCosl     =  0.90; 	// |cos(theta_l)| maximum
Double_t xCoslw    =  0.75; 	// -Q_l*cos(theta_l) maximum
Double_t xMll      = 10.00; 	// |m_ll-m_Z| minimum
Double_t xAcop     = 30.00; 	// Acoplanarity minimum

// Some constants //

static const Double_t kMassW   = 80.00; // W mass
static const Double_t kMassZ   = 91.19; // Z mass
static const Double_t kSigmaMw =   4.0; // W mass resolution
static const Double_t kSigmaMz =   4.0; // W mass resolution
static const Int_t    kZoneX   =     3;	// No. X Zones in the Canvas
static const Int_t    kZoneY   =     3;	// No. Y Zones in the Canvas

// Some global variables //

static Int_t Ngoods = 0;	// Number of good events

JSFQuickSim *sim=0;
JSFSIMDST *simdst=0;

//________________________________________________________
void UserModuleDefine()
{
  Char_t *outputFileName=jsf->Env()->GetValue("JSFGUI.OutputFileName","jsf.root");
  Char_t *inputFileName=jsf->Env()->GetValue("JSFGUI.InputFileName","");

  ofile = new TFile(outputFileName,"RECREATE");
  file  = new TFile(inputFileName);
  jsf->SetIOFiles();
  jsf->SetOutput(*ofile);

  sim=new JSFQuickSim();
  simdst=new JSFSIMDST();
  simdst->SetQuickSimParam(sim->Param());
}

//_________________________________________________________
void UserInitialize()
{
  //  This function is called at the begining of the job or when
  //  "reset hist" action is selected in the gui menu.
  //  This is used to define/reset histograms.

  hStat     = new TH1F("hStat","Cut Statistics"   ,  20,  0.0,  20.0);
  hNleptons = new TH1F("hNleptons","No. leptons"  ,  20,  0.0,  20.0);
  hElepton  = new TH1F("hElepton","Lepton Energy" ,  50,  0.0, 100.0);
  hEvis     = new TH1F("hEvis","Visible energy"   ,  30,  0.0, 150.0);
  hPt       = new TH1F("hPt","Missing Pt"         ,  50,  0.0, 100.0);
  hCoslm    = new TH1F("hCoslm","cos(theta_lm)"   ,  50, -1.0,  +1.0);
  hCoslp    = new TH1F("hCoslp","cos(theta_lp)"   ,  50, -1.0,  +1.0);
  hMll      = new TH1F("hMll","m_ll"              ,  50, 50.0, 150.0);
  hAcop     = new TH1F("hAcop","Acoplanarity"     ,  90,  0.0, 180.0);
}

//_________________________________________________________
void DrawHist()
{
  //  This function is called to draw histograms during the interactive 
  //  session.  Thus you can see the accumulation of the histogram
  //  interactively.  

  TDirectory *last = gDirectory;
  if( !cHist ) {
    cHist = new TCanvas("cHist","Canvas 1",50, 50, kZoneX*200, kZoneY*200);
    cHist->Divide(kZoneX,kZoneY);
  } 
  else {
    cHist->cd();
  }
  cHist->cd(1);	hStat->Draw();
  cHist->cd(2);	hNleptons->Draw();
  cHist->cd(3);	hElepton->Draw();
  cHist->cd(4);	hEvis->Draw();
  cHist->cd(5);	hPt->Draw();
  cHist->cd(6);	hCoslm->Draw();
  cHist->cd(7);	hCoslp->Draw();
  cHist->cd(8);	hMll->Draw();
  cHist->cd(9);	hAcop->Draw();

  cHist->Update();
  last->cd();
}

//_________________________________________________________
void UserSetOptions()
{
  // This function is called only once, soon after jsf is started.
  // This function can be used to define parameters which is not 
  // defined in jsf.conf file.

}

//_________________________________________________________
void UserAnalysis()
{
  // This function is called when the processing of one event is completed.
  // Any data processing of the event can be performed in this function.
  //

  Char_t msg[60];

  // Analysis starts here.
  
  Float_t selid = -0.5;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) strcpy(&cutName[(Int_t)selid][0],"No cut");

  // Get event buffer and make combined tracks accessible.

  JSFSIMDST    *sds     = (JSFSIMDST*)jsf->FindModule("JSFSIMDST");
  JSFSIMDSTBuf *evt     = (JSFSIMDSTBuf*)sds->EventBuf();
  Int_t         ntracks = evt->GetNLTKCLTracks(); 	// No. of tracks 
  TObjArray *tracks     = evt->GetLTKCLTracks(); 	// combined tracks

  // Cut on No. of tracks.
  
  if ( ntracks != 2 ) return;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"N_tracks = 2");
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Loop over combined tracks.
  
  ANL4DVector qsum;
  ANL4DVector qm;
  ANL4DVector qp;
   
  fNleptons = 0;
  for ( Int_t i = 0; i < ntracks; i++ ) {
    JSFLTKCLTrack *t = (JSFLTKCLTrack*)tracks->UncheckedAt(i);
    if ( t->GetType() == LEPTONID ) {
      fNleptons++;			// No. of leptons
      ANL4DVector qt(t->GetPV());	// track 4 momentum
      if ( t->GetCharge() < 0 ) {
        qm = qt;			// track 4 momentum
        fElepm = qm(0);			// l^- energy
        fCoslm = qm(3)/qm.GetMag();	// l^- cos(theta)
      } else {
        qp = qt;			// track 4 momentum
        fElepp = qp(0);			// l^+ energy
        fCoslp = qp(3)/qp.GetMag();	// l^+ cos(theta)
      }
      qsum += qt;
    }
  }
  fEvis = qsum(0);		// E_vis
  fPt   = qsum.GetPt();		// P_t
  fMll  = qsum.GetMass();	// m_ll
  fAcop = qm.Acop(qp);		// theta_Acop (degrees)
   
  if (gDEBUG) {
    cerr << "Evis = " << fEvis << " Pt = " << fPt 
    	 << " M_ll = " << fMll << " Acop = " << fAcop << endl;
    cerr << "l-: "; qm.DebugPrint();
    cerr << "l+: "; qp.DebugPrint();
    cerr << "ll: "; qsum.DebugPrint();
  }

  // Cut on No. of leptons.
  
  hNleptons->Fill(fNleptons);
  if ( fNleptons != 2 ) return;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"N_leptons = 2");
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on lepton energies.

  hElepton->Fill(fElepm);
  hElepton->Fill(fElepp);
  if ( fElepm < xElepton || fElepp < xElepton ) return;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"E_lepton > %g",xElepton);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Evis.

  hEvis->Fill(fEvis);
  if ( fEvis < xEvis ) return;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"E_vis > %g",xEvis);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
 
  // Cut on Pt.

  hPt->Fill(fPt);
  if ( fPt < xPt ) return;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Pt > %g",xPt);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
   
  // Cut on |cos(theta_l)|.

  hCoslm->Fill(fCoslm);
  hCoslp->Fill(fCoslp);
  if ( TMath::Abs(fCoslm) > xCosl || TMath::Abs(fCoslp) > xCosl ) return;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_l)| < %g",xCosl);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
   
  // Cut on -Q_l*cos(theta_l).

  if ( fCoslm > xCoslw && -fCoslp > xCoslw ) return;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"-Q*cos(theta_l) < %g",xCoslw);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
   
  // Cut on m_ll.

  hMll->Fill(fMll);
  if ( (( kMassZ - xMll ) < fMll) && (fMll < ( kMassZ + xMll )) ) return;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"M_ll > %g",xMll);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  // Cut on Acop.

  hAcop->Fill(fAcop);
  if ( fAcop < xAcop ) return;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Acop > %g",xAcop);
    strcpy(&cutName[(Int_t)selid][0],msg);
    ++selid;
    strcpy(&cutName[(Int_t)selid][0],"END");
  }

 Ngoods++;
 
  cerr << "------------------------------------------" << endl
       << "Event " << jsf->GetEventNumber()  << endl
       << "------------------------------------------" << endl;
  //  Uncomment the following to print generator particle list.

  /*
  cerr << " # Generator Particles = " << sds->GetNGeneratorParticles() << endl;
  TClonesArray *gen = sds->GetGeneratorParticles();
  for ( Int_t i = 0; i<sds->GetNGeneratorParticles(); i++ ) {
    JSFGeneratorParticle *g = gen->UncheckedAt(i);
    Int_t ndau = g->GetNDaughter();
    if ( ndau != 0 ) continue;
    // cerr << " ndau = " << ndau << endl;
    g->ls();
  }
  */

 }

//_________________________________________________________
void UserTerminate()
{
  // This function is called at the end of job.
  cout << endl;
  cout << "  =============" << endl;
  cout << "   Cut Summary " << endl;
  cout << "  =============" << endl;
  cout << endl;
  cout << "  -----------------------------------------------------------" << endl;
  cout << "   ID   No.Events    Cut Description" << endl;
  cout << "  -----------------------------------------------------------" << endl;
  Int_t i;
  for ( i = 0; strncmp(&cutName[i][0],"END",4) && i < maxcut ; i++ ) {
    printf("  %3d  %10d  : %s\n",i,(int)hStat->GetBinContent(i+1),&cutName[i][0]);
  } 
  cout << "  -----------------------------------------------------------" << endl;
}







