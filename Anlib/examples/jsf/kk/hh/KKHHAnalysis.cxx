//***************************************************************************
//*  ====================
//*  KKHHAnalysis Classes
//*  ====================
//*
//*  (Description)
//*	A user analysis class for JLC analyses.
//*	This reads and analyzes MC e+e- -> hh data.
//*  (Requires)
//*	library Anlib (in physsim-99a-1 if K.Fujii)
//*  (Provides)
//*	class KKHHAnalysis
//*	class KKHHAnalysisBuf
//*  (Usage)
//*	...
//*  (Update Record)
//*	09 Dec 2003	Nicolas Delerue Adapted from ZH analysis
//***************************************************************************

#include "KKHHAnalysis.h"

static const Double_t kMassH   = 120.0; // H mass
//static const Double_t kMassZ   = 91.19;	// Z mass
static const Double_t kSigmaMh =   4.0;	// H mass resolution
//static const Double_t kSigmaMz =   4.0;	// Z mass resolution
static const Int_t    kZoneX   =     4;	// No. of X Zones in the Canvas
static const Int_t    kZoneY   =     5;	// No. of Y Zones in the Canvas

Int_t KKHHAnalysis::Ngoods = 0;
Bool_t gDEBUG = kFALSE;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ---------------------
//  KKHHAnalysisBuf Class
//  ---------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ClassImp(KKHHAnalysisBuf)

// Constructors
KKHHAnalysisBuf::KKHHAnalysisBuf(const Char_t *name, const Char_t *title,
   KKHHAnalysis *module) : JSFEventBuf(name, title, (JSFModule*)module) {}

KKHHAnalysisBuf::KKHHAnalysisBuf(KKHHAnalysis *module, const Char_t *name,
   const Char_t *title) : JSFEventBuf(name, title, (JSFModule*)module) {}

KKHHAnalysisBuf::~KKHHAnalysisBuf() {}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//  ------------------
//  KKHHAnalysis Class
//  ------------------
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

ClassImp(KKHHAnalysis)

// Constructor
KKHHAnalysis::KKHHAnalysis(const Char_t *name, const Char_t *title)
               : JSFModule(name, title)
{
  fEventBuf = new KKHHAnalysisBuf(this);
  SetBufferSize(2000);	// buffer size for event data
  cout << "KKHHAnalysis is created... fEventBuf is "
       << (Int_t)fEventBuf << endl;

  //We set the cross section to a non physical value
  crossSection=-1;
}

// Destructor
KKHHAnalysis::~KKHHAnalysis()
{
  cout << "KKHHAnalysisBuf will be deleted... fEventBuf is "
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}


//####### Functions ###########

// *****CleanUp()***** //
void KKHHAnalysis::CleanUp(TObjArray *objs)
{
  TIter next(objs);
  TObject *obj;
  while ( (obj = next()) ) {
    objs->Remove(obj);
    delete obj;
  }
}
// *****End of CleanUp()***** //


// ******Initialize()****** //
Bool_t KKHHAnalysis::Initialize()
{
  TDirectory *last = gDirectory;
  gFile->cd("/");

  hStat       = new TH1F("hStat","Cut Statistics"   ,  20,   0.0,   20.0);
  hNtracks    = new TH1F("hNtracks","No. of Tracks" ,  50,   0.0,  100.0);
  hEvis       = new TH1F("hEvis","Visible Energy"   ,  50,   0.0,  2000.0);
  hPt         = new TH1F("hPt","Missing Pt"         ,  50,   0.0,  250.0);
  hNjets      = new TH1F("hNjets","No. of Jets"     ,  20,   0.0,   20.0);
  hEjet       = new TH1F("hEjet","Jet Energy"       ,  50,   0.0,  1000.0);
  hCosjet     = new TH1F("hCosjet","cos(theta_j)"   ,  50,  -1.0,   +1.0);
  hNsols      = new TH1F("hNsols","No. of solutions",  20,   0.0,   20.0);
  hChi2       = new TH1F("hChi2","Chi2"             ,  50,   0.0,   50.0);
  hEvisPl     = new TH2F("hEvisPl","(Evis,Pl)"      ,
                                     60,  0.0, 600.0,  50,-100.0, +100.0);
  hEzEh       = new TH2F("hEzEh","E_Z vs. E_H"      ,
				     50,  0.0, 500.0,  50,   0.0,  500.0);
  hCoszCosh   = new TH2F("hCoszCosh","Cos(theta_Z) vs. Cos(theta_H)",
                                     50, -1.0,  +1.0,  50,  -1.0,   +1.0);
  hAcop       = new TH1F("hAcop","Acoplanarity"     ,  90,   0.0,  180.0);
  hAcopZoom       = new TH1F("hAcopZoom","Acoplanarity"     ,  90,   0.0,  18.0);
  hMassh1      = new TH1F("hMassh1","M_h Distribution",  60,   60.0,  240.0);
  hMassh2      = new TH1F("hMassh2","M_h Distribution",  60,   60.0,  240.0);
  hMasshTot      = new TH1F("hMasshTot","M_h Distribution (squared avg)",  60,   60.0,  240.0);
  hMasshh      = new TH2F("hMasshh","M_h (h1 vs h2)Distribution",  60,   60.0,  240.0,  60,   60.0,  240.0);



  hMassh1raw      = new TH1F("hMassh1raw","M_h Distribution (raw)",  60,   60.0,  240.0);
  hMassh2raw      = new TH1F("hMassh2raw","M_h Distribution (raw)",  60,   60.0,  240.0);
  hMasshhraw      = new TH2F("hMasshhraw","M_h (h1 vs h2)Distribution",  60,   60.0,  240.0,  60,   60.0,  240.0);

  
  hMassh1ac      = new TH1F("hMassh1ac","M_h Distribution (after mjj cut)",  40,   80.0,  160.0);
  hMassh2ac      = new TH1F("hMassh2ac","M_h Distribution (after mjj cut)",  40,   80.0,  160.0);
  hMasshhac      = new TH2F("hMasshhac","M_h (h1 vs h2)Distribution",  40,   80.0,  160.0,  40,   80.0,  160.0);
  

  hCosTheta     = new TH1F("hTheta","cos(theta_h)"   ,  50,  -1.0,   +1.0);


  xNtracks  =     25;	// No. of Tracks
  xEtrack   =   0.10;	// track energy
  xEvis     = 100.00;	// Minimum visible energy
  xPt       =  10.00;	// Pt minimum
  xPl       = 999.00;	// Pl maximum
  xYcut     =  0.004;	// y_cut to force the event to 4 jets
  xNjets    =      4;	// No. of Jets
  xEjet     =   5.00;	// E_jet minimum
  xCosjet   =   0.99;	// |cos(theta_j)| maximum
  xCoszh    =   0.99;	// |cos(theta_Z)| and |cos(theta_H)| maximum
  xM2j      =  18.00;	// |m_jj-m_S| maximum , S = Z,H
  xAcop     =  30.00;	// Acoplanarity maximum

  last->cd();
  return 0;
}
// *****End of Initialize()***** //


// *****DrawHist()***** //
void KKHHAnalysis::DrawHist()
{
  TDirectory *last = gDirectory;
  if (!cHist) {
    cHist = new TCanvas("cHist","Canvas 1", 10, 10, kZoneX*200, kZoneY*200);
    cHist->Divide(kZoneX,kZoneY);
  } 
  else {
    cHist->cd();
  }

  cHist->cd(1);		hStat->Draw();
  cHist->cd(2);		hNtracks->Draw();
  cHist->cd(3);		hEvis->Draw();
  cHist->cd(4);		hPt->Draw();
  cHist->cd(5);		hNjets->Draw();
  cHist->cd(6);		hEjet->Draw();
  cHist->cd(7);		hCosjet->Draw();
  cHist->cd(8);		hNsols->Draw();
  cHist->cd(9);		hChi2->Draw();
  cHist->cd(10);	hMassh1->Draw();
  cHist->cd(11);	hMassh2->Draw();
  cHist->cd(12);	hMasshh->Draw();
  cHist->cd(13);	hEvis->Draw();
  cHist->cd(14);	hAcop->Draw();
  cHist->cd(15);	hAcopZoom->Draw();
  cHist->cd(17);	hMassh1raw->Draw();
  cHist->cd(18);	hMassh2raw->Draw();
  cHist->cd(19);	hMasshhraw->Draw();

  cHist->Update();

  last->cd();
}
// *****End of DrawHist()***** //


// *****Process()***** //
Bool_t KKHHAnalysis::Process(Int_t ev)
{
  // Local copies of KKHHAnalysisBuf data members.

  Int_t		fNtracks;	// track multiplicity
  Double_t	fEvis;		// visible energy
  Double_t	fPt;		// Pt
  Double_t	fPl;		// Pl
  Double_t	fYcut;		// y_cut to force the event to 4 jets
  Int_t		fNjets;		// jet multiplicity

  // Remember the previous directory.
  
  TDirectory *last = gDirectory;
  gFile->cd("/");

  Char_t msg[60];

  // Analysis starts here.

  Float_t selid = -0.5;
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) strcpy(&cutName[(Int_t)selid][0],"No Cut");

  // Get event buffer and make combined tracks accessible.

  JSFSIMDST     *sds	= (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
  JSFSIMDSTBuf	*evt	= (JSFSIMDSTBuf*)sds->EventBuf();
  KKHHAnalysisBuf *ua	= (KKHHAnalysisBuf*)fEventBuf;
  KKHHAnalysisBuf &a	= *ua;

  Int_t		ntrks	= evt->GetNLTKCLTracks();	// No. of tracks
  TObjArray	*trks	= evt->GetLTKCLTracks();	// Combined tracks

  // Select good tracks
  ANL4DVector qsum;
  TObjArray tracks(1000);
  fNtracks = 0;
  for ( Int_t i = 0; i < ntrks; i++) {
    JSFLTKCLTrack *t = (JSFLTKCLTrack*)trks->UncheckedAt(i);
    if ( t->GetE() > xEtrack ) {
      ANLTrack *qt = new ANLTrack(t);
      tracks.Add(qt);		// track 4-momentum
      qsum += *qt;		// total 4-mometum
      fNtracks++;		// *qt stays.
    }
  }
  if (gDEBUG) cerr << "Ntracks = " << fNtracks << endl;

  // Cut on No. of tracks.
  
  hNtracks->Fill(fNtracks);
  if ( fNtracks < xNtracks ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"N_tracks > %i",xNtracks);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  
  fEvis = qsum(0);	// E_vis
  fPt	= qsum.GetPt(); // P_t
  fPl	= qsum(3);	// P_l

  if (gDEBUG) cerr << "Evis = " << fEvis << " Pt = "
	<< fPt << " Pl = " << fPl << endl;

  // Cut on Evis.

  hEvis->Fill(fEvis);
  if ( fEvis < xEvis ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Evis > %g",xEvis);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Pt.
  
  hPt->Fill(fPt);
  if ( fPt > xPt ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg," Pt <= %g",xPt);
    strcpy(&cutName[(Int_t)selid][0], msg);
  }

  // Cut on Pl.
  
  if ( TMath::Abs(fPl) > xPl ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|Pl| <= %g",xPl);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Find jets.

  fYcut = xYcut;
  ANLJadeEJetFinder jclust(fYcut);
  jclust.Initialize(tracks);
  jclust.FindJets();
  fYcut  = jclust.GetYcut();
  fNjets = jclust.GetNjets();
  
  if (gDEBUG) cerr << "Ycut = " << fYcut << " Njets = " << fNjets << endl;

  // Cut on No. of jets.

  hNjets->Fill(fNjets);
  if ( fNjets < xNjets ) { CleanUp(&tracks); return kFALSE;}
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Njets >= %i for Ycut = %g",xNjets,xYcut);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Now force the event to be xNjets.

  jclust.ForceNJets(xNjets);
  fNjets = jclust.GetNjets();
  fYcut  = jclust.GetYcut();

  if(gDEBUG) cerr << "Ycut = " << fYcut << " Njets = " << fNjets << endl;

  // Make sure the No. of Jets is xNjets.
  
  if ( fNjets != xNjets ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Njets = %i",xNjets);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  TObjArray &jets = jclust.GetJets();
  TIter nextjet(&jets);
  ANLJet *jetp;
  Double_t ejetmin = 999999.;
  Double_t cosjmax = 0.;
  while ( (jetp = (ANLJet*)nextjet()) ) {
    ANLJet &jet = *jetp;
    if (gDEBUG) jet.DebugPrint();
    Double_t ejet = jet()(0);
    if ( ejet < ejetmin ) ejetmin = ejet;
    hEjet->Fill(ejet);
    Double_t cosj = jet.CosTheta();
    if ( TMath::Abs(cosj) > TMath::Abs(cosjmax) ) cosjmax = cosj;
    hCosjet->Fill(cosj);
  }

  // Cut on Ejet_min.

  if ( ejetmin < xEjet ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Ejet > %g",xEjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on |cos(theta_j)|_max.

  if ( TMath::Abs(cosjmax) > xCosjet ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_j)| <= %g", xCosjet);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }




  // Find H candidates in given relaxed mass window.
  
  TObjArray solutions(10);
  ANLPairCombiner hcandidates2(jets,jets);
  ANLPair *hp1, *hp2;
  while ( (hp2 = (ANLPair*)hcandidates2()) ) {
    ANLPair &h2 = *hp2;
    Double_t hmass2 = h2().GetMass();
    if (TMath::Abs(hmass2 - kMassH) > (4*xM2j)) continue;
    h2.LockChildren();
    ANLPairCombiner hcandidates1(hcandidates2);
    while ( (hp1 = (ANLPair*)hcandidates1()) ) {
      ANLPair &h1 = *hp1;
      if (h1.IsLocked()) continue;
      Double_t hmass1 = h1().GetMass();
      if (TMath::Abs(hmass1 - kMassH) > (4*xM2j)) continue;
      if (gDEBUG) {
        cerr << " M_h2 = " << hmass2 << " M_h1 = " << hmass1 << endl;
        cerr << " h2  = " << (void*)hp2
             << " h1  = " << (void*)hp1 << endl;
        cerr << " h2[0] = " << (void*)h2[0]
             << " h2[1] = " << (void*)h2[1]
             << " h1[0] = " << (void*)h1[0]
             << " h1[1] = " << (void*)h1[1] << endl;
      }
      Double_t chi2 = TMath::Power((hmass2 - kMassH)/kSigmaMh,2.)
                    + TMath::Power((hmass1 - kMassH)/kSigmaMh,2.);
      hMassh1raw->Fill(hmass1);
      hMassh2raw->Fill(hmass2);
      if (hmass1>hmass2) {
	hMasshhraw->Fill(hmass1,hmass2);
      } else {
	hMasshhraw->Fill(hmass2,hmass1);
      }
      hMasshTot->Fill((hmass2+hmass1)/2);
      solutions.Add(new ANLPair(hp2,hp1,chi2));
    }
    h2.UnlockChildren();
  }

  // Cut on number of solutions.
 
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"2 jets pairs");
    strcpy(&cutName[(Int_t)selid][0],msg);
  }


  //Tighten the mass window
  TIter hmasssol(&solutions);
  ANLPair *sol;
  while ( (sol = (ANLPair*)hmasssol()) ) {
    ANLPair &h1 = *(ANLPair*)(*sol)[0];
    ANLPair &h2 = *(ANLPair*)(*sol)[1];

    Double_t hmass1 = h1().GetMass();
    if (TMath::Abs(hmass1 - kMassH) > xM2j) {
      solutions.Remove(sol);
      delete sol;
      continue;
    }
        
    Double_t hmass2 = h2().GetMass();
    if (TMath::Abs(hmass2 - kMassH) > xM2j) {
      solutions.Remove(sol);
      delete sol;
      continue;
    }
        
  } // for each solution


  // Cut on number of solutions.
 
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|m_jj - m_h| <= %g",xM2j);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  
  //Filling histos after cuts
  if (1==1) {
    TIter hmasssol(&solutions);
    ANLPair *sol;
    while ( (sol = (ANLPair*)hmasssol()) ) {
      ANLPair &h1 = *(ANLPair*)(*sol)[0];
      ANLPair &h2 = *(ANLPair*)(*sol)[1];
      Double_t hmass1 = h1().GetMass();
      Double_t hmass2 = h2().GetMass();
      hMassh1ac->Fill(hmass1);
      hMassh2ac->Fill(hmass2);
      if (hmass1>hmass2) {
	hMasshhac->Fill(hmass1,hmass2);
      } else {
	hMasshhac->Fill(hmass2,hmass1);
      }
    }  // for each solution
  } // 1 == 1

  //B Tagging
  TIter btagsol(&solutions);
  ANLVTXTagger btag(3.,3);
  while ( (sol = (ANLPair*)btagsol()) ) {
    ANLPair &h1 = *(ANLPair*)(*sol)[0];
    ANLPair &h2 = *(ANLPair*)(*sol)[1];
    if (((!btag(*(ANLJet *)h1[0]) && !btag(*(ANLJet *)h1[1])))&&((!btag(*(ANLJet *)h2[0]) && !btag(*(ANLJet *)h2[1])))){
      solutions.Remove(sol);
      delete sol;
      continue;
    }
    
    if ((!btag(*(ANLJet *)h1[0]) && !btag(*(ANLJet *)h1[1]))) {
      solutions.Remove(sol);
      delete sol;
      continue;
    }
    
    
    if ((!btag(*(ANLJet *)h2[0]) && !btag(*(ANLJet *)h2[1]))) {
      solutions.Remove(sol);
      delete sol;
      continue;
    }
    
  } // for each solution


  // Cut on number of solutions.
 
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"B tagging 3.,3");
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on cos(theta_ZH).

  TIter nextsol(&solutions);
  while ( (sol = (ANLPair*)nextsol()) ) {
    ANL4DVector &h2 = *(ANL4DVector*)(*sol)[0];
    ANL4DVector &h1 = *(ANL4DVector*)(*sol)[1];
    Double_t eh2 = h2(0);
    Double_t eh1 = h1(0);
    hEzEh->Fill(eh2,eh1,1.0);
    Double_t cosh2 = h2.CosTheta();
    Double_t cosh1 = h1.CosTheta();
    hCoszCosh->Fill(cosh2,cosh1,1.);
    if (TMath::Abs(cosh2) > xCoszh || TMath::Abs(cosh1) > xCoszh) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"|cos(theta_h)| and |cos(theta_h)| <= %g", xCoszh);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }

  // Cut on Acop.

  nextsol.Reset();
  while ( (sol = (ANLPair*)nextsol()) ) {
    ANL4DVector &h2 = *(ANL4DVector*)(*sol)[0];
    ANL4DVector &h1 = *(ANL4DVector*)(*sol)[1];
    Double_t acop = h2.Acop(h1);
    hAcop->Fill(acop);
    hAcopZoom->Fill(acop);
    if (acop > xAcop) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Acop <= %g",xAcop);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }


  //Cut on total distance from Higgs Mass
  nextsol.Reset();
  while ( (sol = (ANLPair*)nextsol()) ) {
    ANLPair &h2 = *(ANLPair*)(*sol)[0];
    ANLPair &h1 = *(ANLPair*)(*sol)[1];
    Double_t h2mass = h2.GetMass();
    Double_t h1mass = h1.GetMass();
    if (TMath::Abs(((h1mass+h2mass)/2)-kMassH)>xMassTotDist) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if ( !solutions.GetEntries() ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if ( Ngoods == 0 ) {
    sprintf(msg,"Total distance from Higgs Mass %g",xMassTotDist);
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
    
 




  // End of event selection
  
  if ( Ngoods == 0 ) {
    selid++;
    sprintf(msg,"END");
    strcpy(&cutName[(Int_t)selid][0],msg);
  }
  Ngoods++;
  
    if (solutions.GetEntries()>1) {
  cerr << "------------------------------------------" << endl
       << "Event " << gJSF->GetEventNumber()
       << ": Number of solutions = " << solutions.GetEntries() << endl
       << "------------------------------------------" << endl;
   }

  // Sort the solutions in the ascending order of chi2 values.

  solutions.Sort();

  // Now store this in KKHHAnalysisBuf.

  a.fNtracks	= fNtracks;
  a.fEvis	= fEvis;
  a.fPt		= fPt;
  a.fPl		= fPl;
  a.fYcut	= fYcut;
  a.fNjets	= fNjets;

  // Hists and plots for selected events.

  if (gDEBUG) {
    Int_t nj = 0;
    nextjet.Reset();
    while ( (jetp = (ANLJet*)nextjet()) ) {
      cerr << "------" << endl
           << "Jet " << ++nj << endl
           << "------" << endl;
      jetp->DebugPrint();
    }
  }

  hNsols->Fill(solutions.GetEntries());
  hEvisPl->Fill(fEvis,fPl,1.);

  nextsol.Reset();
  Int_t nsols = 0;
  while ( (sol = (ANLPair*)nextsol()) ) {
    if ( nsols++ ) break;
    ANLPair &h2 = *(ANLPair*)(*sol)[0];
    ANLPair &h1 = *(ANLPair*)(*sol)[1];
    Double_t chi2   = sol->GetQuality();
    Double_t h2mass = h2.GetMass();
    hMassh2->Fill(h2mass);
    Double_t h1mass = h1.GetMass();
    hMassh1->Fill(h1mass);
    hMasshh->Fill(h1mass,h2mass);
    hChi2->Fill(chi2);
    hCosTheta->Fill(h1.CosTheta());
    hCosTheta->Fill(h2.CosTheta());
  }



  // Clean up

  CleanUp(&solutions);
  CleanUp(&tracks);

  last->cd();
  return kTRUE;
}
// *****End of Process()***** //


// *****Terminate()***** //
Bool_t KKHHAnalysis::Terminate()
{
  // This function is called at the end of job.
  cout << endl;
  cout << "  =============" << endl;
  cout << "   Cut Summary " << endl;
  cout << "  =============" << endl;
  cout << endl;
  cout << "Cross section= " << crossSection << " lumi= " << (hStat->GetBinContent(1)/crossSection) << endl; 
  cout << endl;
  cout << "  ---------------------------------------------------------" <<endl;   
  cout << "  ID   No.Events     Scaled     Cut Description" << endl;
  cout << "  ---------------------------------------------------------" <<endl;
  Int_t i;
  for ( i = 0; strncmp(&cutName[i][0],"END",4) && i < MAXCUT ; i++ ) {
    printf("  %3d  %10d  %4f \t   : %s\n",i,(int)hStat->GetBinContent(i+1),((Double_t)((hStat->GetBinContent(i+1)*crossSection)/hStat->GetBinContent(1))),&cutName[i][0]);
  } 
  cout << "  ---------------------------------------------------------" << endl;
  
  return 0;
}
// *****End of Terminate()***** //
