//* $Id$
//*************************************************************************
//* ========================
//*  SFSF2LAnalysis Classes
//* ========================
//*
//* (Description)
//*   A sample user analysis classes for JLC analyses.
//*   This reads and analyzes sfsf data to select 2-lepton events.
//* (Requires)
//* 	library Anlib
//* 	library SFSFStudy
//* (Provides)
//* 	class SFSF2LAnalysis
//* 	class SFSF2LAnalysisBuf
//* (Usage)
//*   Take a look at anl2L2J.C.
//* (Update Recored)
//*
//*************************************************************************
//
#include "SFSF2LAnalysis.h"
#include "JSFSpringParton.h"
#include "JSFSpring.h"
#include "ANLTrack.h"

static const Double_t kMassW   = 80.00; 	// W mass
static const Double_t kMassZ   = 91.19; 	// Z mass
static const Double_t kMasst   = 170.0; 	// top mass
static const Double_t kSigmaMw =   4.0; 	// W mass resolution
static const Double_t kSigmaMz =   4.0; 	// Z mass resolution
static const Double_t kSigmaMt =  15.0; 	// top mass resolution
static const Int_t    kZoneX   =     3;		// No. X Zones in the Canvas
static const Int_t    kZoneY   =     3;		// No. Y Zones in the Canvas

Int_t SFSF2LAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

typedef enum { kElectron = 11, kMuon = 13 } EPID;

//_____________________________________________________________________
//  -----------------------
//  SFSF2LAnalysisBuf Class
//  -----------------------
//
//
ClassImp(SFSF2LAnalysisBuf)

//_________________________________________________________
SFSF2LAnalysisBuf::SFSF2LAnalysisBuf(const Char_t         *name, 
                                     const Char_t         *title,
                                           SFSF2LAnalysis *module) 
      : JSFEventBuf(name, title, (JSFModule*)module),
        fSpringName(0), fEcm(350.), fNtracks(-9999),
        fNleptons(-9999), fElepm(-9999.), fElepp(-9999.), fEvis(-9999.),
        fPt(-9999.), fPl(-9999.), fCoslm(-9999.), fCoslp(-9999.),
        fCosSFm1(-9999.), fCosSFm2(-9999.), fCosSFmG(-9999.),
        fMll(-9999.), fAcop(-9999.)
{
   SetSpringName("SFSFSpring");
}

//_________________________________________________________
SFSF2LAnalysisBuf::SFSF2LAnalysisBuf(      SFSF2LAnalysis *module, 
                                     const Char_t         *name,
                                     const Char_t         *title) 
      : JSFEventBuf(name, title, (JSFModule*)module),
        fSpringName(0), fEcm(350.), fNtracks(-9999),
        fNleptons(-9999), fElepm(-9999.), fElepp(-9999.), fEvis(-9999.),
        fPt(-9999.), fPl(-9999.), fCoslm(-9999.), fCoslp(-9999.),
        fCosSFm1(-9999.), fCosSFm2(-9999.), fCosSFmG(-9999.),
        fMll(-9999.), fAcop(-9999.)
{
   SetSpringName("SFSFSpring");
}

//_____________________________________________________________________
//  --------------------
//  SFSF2LAnalysis Class
//  --------------------
//
//

ClassImp(SFSF2LAnalysis)

SFSF2LAnalysis::SFSF2LAnalysis(const Char_t *name, const Char_t *title)
  : JSFModule(name, title), cHist(0)
{
  fEventBuf = new SFSF2LAnalysisBuf(this);
  SetBufferSize(2000);  // buffer size for event data.
  cout << "SFSF2LAnalysisBuf is created...fEventBuf is "
       << (Int_t)fEventBuf << endl;
}

//_____________________________________________________________________
SFSF2LAnalysis::~SFSF2LAnalysis()
{
  cout << "SFSF2LAnalysisBuf will be deleted...fEventBuf is "
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}

//_____________________________________________________________________
Bool_t SFSF2LAnalysis::Initialize()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat     = new TH1F("hStat","Cut Statistics"   ,  20,  0.0,  20.0);
  hNtracks  = new TH1F("hNtracks","No. tracks"    ,  20,  0.0,  20.0);
  hNleptons = new TH1F("hNleptons","No. leptons"  ,  20,  0.0,  20.0);
  hElepton  = new TH1F("hElepton","Lepton Energy" , 200,  0.0, 100.0);
  hEvis     = new TH1F("hEvis","Visible energy"   ,  30,  0.0, 150.0);
  hPt       = new TH1F("hPt","Missing Pt"         ,  50,  0.0, 100.0);
  hCoslm    = new TH1F("hCoslm","cos(theta_lm)"   ,  50, -1.0,  +1.0);
  hCoslp    = new TH1F("hCoslp","cos(theta_lp)"   ,  50, -1.0,  +1.0);
  hMll      = new TH1F("hMll","m_ll"              ,  50, 50.0, 150.0);
  hAcop     = new TH1F("hAcop","Acoplanarity"     ,  90,  0.0, 180.0);

  hMllFinal      = new TH1F("hMllFinal","m_ll final"      ,  50, 50.0, 150.0);
  hElmFinal      = new TH1F("hElmFinal","l- Energy final" , 200,  0.0, 100.0);
  hElpFinal      = new TH1F("hElpFinal","l+ Energy final" , 200,  0.0, 100.0);
  hElFinal       = new TH1F("hElFinal","l Energy final"   , 200,  0.0,  50.0);
  hCoslmFinal    = new TH1F("hCoslmFinal","cos_l^- final" ,  50, -1.0,  +1.0);
  hCoslpFinal    = new TH1F("hCoslpFinal","cos_l^+ final" ,  50, -1.0,  +1.0);
  hPtFinal       = new TH1F("hPtFinal","Missing Pt final" ,  50,  0.0, 100.0);
  hAcopFinal     = new TH1F("hAcopFinal","Acop final"     ,  90,  0.0, 180.0);
  hCosSFmFinal   = new TH1F("hCosSFmFinal","cos_SF^- final", 50, -1.0,  +1.0);
  hEvisPlFinal   = new TH2F("hEvisPlFinal","(Evis,Pl) final", 30, 0.0, 150.,
                                                             100,-200.,200.);
  hCosSFm1G      = new TH2F("hCosSFm1G","(cos(sf-)1,cos(sf-)G)",50,-1.,1.,
                                                                50,-1.,1.);
  hCosSFm2G      = new TH2F("hCosSFm2G","(cos(sf-)2,cos(sf-)G)",50,-1.,1.,
                                                                50,-1.,1.);
  hCosSFm12      = new TH2F("hCosSFm12","(cos(sf-)1,cos(sf-)2)",50,-1.,1.,
                                                                50,-1.,1.);
  hCosSFmGB      = new TH2F("hCosSFmGB","(cos(sf-)2,cos(sf-)G)",50,-1.,1.,
                                                                50,-1.,1.);

  xNleptons =   2   ;     // Number of leptons
  xElLo     =   5.00;     // Lepton E minimum
  xElHi     = 125.00;     // Lepton E minimum
  xEvisLo   =  20.00;     // Minimum visible energy
  xEvisHi   = 250.00;     // Maximum visible energy
  xPt       =  15.00;     // Pt minimum
  xCosl     =   0.90;     // |cos(theta_l)| maximum
  xCoslw    =   0.75;     // -Q_l*cos(theta_l) maximum
  xMll      =  10.00;     // |m_ll-m_Z| minimum
  xAcop     =  30.00;     // Acoplanarity minimum

  for (Int_t i = 0; i < MAXCUT; i++) {
    strcpy(&cutName[i][0],"     ");
  }

  last->cd();
  return 0;
}

//_________________________________________________________
void SFSF2LAnalysis::DrawHist()
{
  TDirectory *last = gDirectory;
  if (!cHist) {
    cHist = new TCanvas("cHist","Canvas 1",10, 10, kZoneX*200, kZoneY*200);
    cHist->Divide(kZoneX,kZoneY);
  } else {
    cHist->cd();
  }

  Int_t Ihist = 0;
  cHist->cd(Ihist++); hStat->Draw();
  cHist->cd(Ihist++); hNleptons->Draw();
  cHist->cd(Ihist++); hElepton->Draw();
  cHist->cd(Ihist++); hEvis->Draw();
  cHist->cd(Ihist++); hPt->Draw();
  cHist->cd(Ihist++); hCoslm->Draw();
  cHist->cd(Ihist++); hCoslp->Draw();
  cHist->cd(Ihist++); hMll->Draw();
  cHist->cd(Ihist++); hAcop->Draw();

  cHist->Update();

  last->cd();
}

//_________________________________________________________
Bool_t SFSF2LAnalysis::Process(Int_t ev)
{
  // Local copies of SFSF2LAnalysisBuf data members.

  Double_t fEcm;          // CM energy (GeV)
  Int_t    fNtracks;      // Multiplicity
  Int_t    fNleptons;     // Lepton multiplicity
  Double_t fElepm;        // l^- energy
  Double_t fElepp;        // l^+ energy
  Double_t fEvis;         // Visible energy
  Double_t fPt;           // Pt
  Double_t fPl;           // Pl
  Double_t fCoslm;        // cos(theta_lepton^-)
  Double_t fCoslp;        // cos(theta_lepton^+)
  Double_t fCosSFm1;      // cos(theta_SF^-)_1
  Double_t fCosSFm2;      // cos(theta_SF^-)_2
  Double_t fCosSFmG;      // cos(theta_SF^-)_G
  Double_t fMll;          // m(l^+l^-)
  Double_t fAcop;         // Acoplanarity
  Double_t fMsf;          // Msf
  Double_t fMlsp;         // Mlsp

  Char_t msg[256];

  // ---------------------
  // Analysis starts here.
  // ---------------------

  Int_t selid = -1;
  hStat->Fill(++selid);
  if (Ngoods == 0) strcpy(&cutName[selid][0],"No cut");

  // Get event buffer and make combined tracks accessible.

  JSFSIMDST         *sds = (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
  JSFSIMDSTBuf      *evt = (JSFSIMDSTBuf*)sds->EventBuf();
  SFSF2LAnalysisBuf *ua  = (SFSF2LAnalysisBuf *)fEventBuf;
  SFSF2LAnalysisBuf &a   = *ua;
  a.SetEcm(evt->GetEcm());
  fEcm     = a.GetEcm();
 
  fNtracks         = evt->GetNLTKCLTracks(); 	// No. of tracks 
  TObjArray  *trks = evt->GetLTKCLTracks(); 	// combined tracks

  // Cut on No. of tracks.

  hNtracks->Fill(fNtracks);
  if (fNtracks < xNtracks) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"N_tracks > %i",xNtracks);
    strcpy(&cutName[selid][0],msg);
  }

  // Select good tracks and store them in "TObjArray tracks".

  ANLTrack qsum;
  ANLTrack qm;
  ANLTrack qp;
  Int_t chargesum = 0;

  fNleptons = 0;
  for ( Int_t i = 0; i < fNtracks; i++ ) {
    JSFLTKCLTrack *t = (JSFLTKCLTrack*)trks->UncheckedAt(i);
    if (t->GetType() == kMuon) {
      fNleptons++;                      // No. of leptons
      ANLTrack qt(t);                   // track 4 momentum
      if (t->GetCharge() < 0) {
        qm = qt;                        // track 4 momentum
        fElepm = qm(0);                 // l^- energy
        fCoslm = qm(3)/qm.GetMag();     // l^- cos(theta)
      } else {
        qp = qt;                        // track 4 momentum
        fElepp = qp(0);                 // l^+ energy
        fCoslp = qp(3)/qp.GetMag();     // l^+ cos(theta)
      }
      qsum      += qt;
      chargesum += t->GetCharge();
    }
  }

  fEvis = qsum(0);		// E_vis
  fPt   = qsum.GetPt();		// P_t
  fPl   = qsum(3);		// P_l
  fMll  = qsum.GetMass();       // m_ll
  fAcop = qm.Acop(qp);          // theta_Acop (degrees)

  if (gDEBUG) {
    cerr << "Evis = " << fEvis << " Pt = " << fPt 
         << " M_ll = " << fMll << " Acop = " << fAcop << endl;
    cerr << "l-: "; qm.DebugPrint();
    cerr << "l+: "; qp.DebugPrint();
    cerr << "ll: "; qsum.DebugPrint();
  }

  // Cut on lepton energies.

  hNleptons->Fill(fNleptons);
  if (fNleptons != 2) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"N_leptons = 2");
    strcpy(&cutName[selid][0],msg);
  }

  // Cut on chage balance.

  if (chargesum != 0) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"charge balance");
    strcpy(&cutName[selid][0],msg);
  }

  // Cut on lepton energies.

  hElepton->Fill(fElepm);
  hElepton->Fill(fElepp);
  if (fElepm < xElLo || fElepm > xElHi ||
      fElepp < xElLo || fElepp > xElHi) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"%g < E_lepton < %g",xElLo,xElHi);
    strcpy(&cutName[selid][0],msg);
  }

  // Cut on Evis.

  hEvis->Fill(fEvis);
  if (fEvis < xEvisLo || fEvis > xEvisHi) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"%g < E_vis < %g",xEvisLo,xEvisHi);
    strcpy(&cutName[selid][0],msg);
  }
   
  // Cut on |cos(theta_l)|.

  hCoslm->Fill(fCoslm);
  hCoslp->Fill(fCoslp);
  if (TMath::Abs(fCoslm) > xCosl || TMath::Abs(fCoslp) > xCosl) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"|cos(theta_l)| < %g",xCosl);
    strcpy(&cutName[selid][0],msg);
  }
   
  // Cut on -Q_l*cos(theta_l).

  if (fCoslm > xCoslw && -fCoslp > xCoslw) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"-Q*cos(theta_l) < %g",xCoslw);
    strcpy(&cutName[selid][0],msg);
  }
   
  // Cut on m_ll.

  hMll->Fill(fMll);
  if (((kMassZ - xMll) < fMll) && (fMll < (kMassZ + xMll))) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"M_ll > %g",xMll);
    strcpy(&cutName[selid][0],msg);
  }
  
  // Cut on Acop.

  hAcop->Fill(fAcop);
  if (fAcop < xAcop) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"Acop > %g",xAcop);
    strcpy(&cutName[selid][0],msg);
  }
 
  // Cut on Pt.

  hPt->Fill(fPt);
  if (fPt < xPt) return kFALSE;
  hStat->Fill(++selid);
  if (Ngoods == 0) {
    sprintf(msg,"Pt > %g",xPt);
    strcpy(&cutName[selid][0],msg);
    ++selid;
    strcpy(&cutName[selid][0],"END");
  }

  // ----------------------
  // End of event selection
  // ----------------------

  if (Ngoods == 0) {
    selid++;
    sprintf(msg,"END");
    strcpy(&cutName[selid][0],msg);
  }
  Ngoods++;

  cerr << "------------------------------------------" << endl
       << "Event " << gJSF->GetEventNumber()           << endl
       << "------------------------------------------" << endl;

  // ----------------------
  // Hists and plots
  // ----------------------

  hMllFinal   ->Fill(fMll);
  hElmFinal   ->Fill(fElepm);
  hElpFinal   ->Fill(fElepm);
  hElFinal    ->Fill(fElepm);
  hElFinal    ->Fill(fElepp);
  hCoslmFinal ->Fill(fCoslm);
  hCoslpFinal ->Fill(fCoslp);
  hPtFinal    ->Fill(fPt);
  hAcopFinal  ->Fill(fAcop);
  hEvisPlFinal->Fill(fEvis,fPl,1.);

  if (!strcmp(a.fSpringName,"SFSFSpring")) {
     JSFSpring       *sp       = (JSFSpring*)gJSF->FindModule(a.fSpringName);
     JSFSpringBuf    *spb      = (JSFSpringBuf*)sp->EventBuf();
     TClonesArray    *partons  = spb->GetPartons();

     JSFSpringParton *sf       = ((JSFSpringParton *)partons->UncheckedAt(0));
     JSFSpringParton *lsp      = ((JSFSpringParton *)partons->UncheckedAt(2));

     fMsf     = sf->GetMass();
     fMlsp    = lsp->GetMass();
     fCosSFmG = sf->GetCosth();

     ANL3DVector esol[2];
     if (a.SolveKinematics(fEcm, fMsf, fMlsp, qm, qp, esol)) {
        if (qm(3) < 0.) {
           fCosSFm1 = TMath::Min(esol[0](3), esol[1](3));
           fCosSFm2 = TMath::Max(esol[0](3), esol[1](3));
        } else {
           fCosSFm1 = TMath::Max(esol[0](3), esol[1](3));
           fCosSFm2 = TMath::Min(esol[0](3), esol[1](3));
        }
     }

     hCosSFmFinal->Fill(fCosSFm1);
     hCosSFmFinal->Fill(fCosSFm2);
     hCosSFm1G   ->Fill(fCosSFm1,fCosSFmG,1.);
     hCosSFm2G   ->Fill(fCosSFm2,fCosSFmG,1.);
     hCosSFm12   ->Fill(fCosSFm1,fCosSFm2,1.);

     if (TMath::Abs(fCosSFm1-fCosSFmG) > TMath::Abs(fCosSFm2-fCosSFmG)) {
        hCosSFmGB->Fill(fCosSFm2,fCosSFm1,1.);
     } else {
        hCosSFmGB->Fill(fCosSFm1,fCosSFm2,1.);
     }
  }

  // ------------------------------------
  // Now store this in SFSF2LAnalysisBuf.
  // ------------------------------------

  a.fEcm      = fEcm;          // CM energy (GeV)
  a.fNtracks  = fNtracks;      // Multiplicity
  a.fNleptons = fNleptons;     // Lepton multiplicity
  a.fElepm    = fElepm;        // l^- energy
  a.fElepp    = fElepp;        // l^+ energy
  a.fEvis     = fEvis;         // Visible energy
  a.fPt       = fPt;           // Pt
  a.fPl       = fPl;           // Pl
  a.fCoslm    = fCoslm;        // cos(theta_lepton^-)
  a.fCoslp    = fCoslp;        // cos(theta_lepton^+)
  a.fCosSFm1  = fCosSFm1;      // cos(theta_SF^-)_1
  a.fCosSFm2  = fCosSFm2;      // cos(theta_SF^-)_2
  a.fCosSFmG  = fCosSFmG;      // cos(theta_SF^-)_G
  a.fMll      = fMll;          // m(l^+l^-)
  a.fAcop     = fAcop;         // Acoplanarity		
  a.fMsf      = fMsf;          // Msf
  a.fMlsp     = fMlsp;         // Mlsp

  return kTRUE;
}

//_________________________________________________________
Bool_t SFSF2LAnalysis::Terminate()
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
  for ( i = 0; strncmp(&cutName[i][0],"END",4) && i < MAXCUT ; i++ ) {
    printf("  %3d  %10d  : %s\n",i,(int)hStat->GetBinContent(i+1),&cutName[i][0]);
  }
  cout << "  -----------------------------------------------------------" << endl;
  return 0;
}
