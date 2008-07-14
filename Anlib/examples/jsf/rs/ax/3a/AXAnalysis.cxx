//***************************************************************************
//*  ====================
//*  AXAnalysis Classes
//*  ====================
//*
//*  (Description)
//*	A user analysis class for JLC analyses.
//*	This reads and analyzes MC e+e- -> AX data.
//*  (Requires)
//*	library Anlib (in physsim-99a-1 if K.Fujii)
//*	library AXStudy (also in physsim)
//*  (Provides)
//*	class AXAnalysis
//*  (Usage)
//*	...
//*  (Update Record)
//*     2007/01/27  K.Fujii	Original version.
//***************************************************************************

#include "AXAnalysis.h"
#include "JSFGeneratorParticle.h"
#include "TNtupleD.h"
#include "ANLTrack.h"

#include <sstream>
#include <iomanip>
//#define __USEGENERATORINFO__

static const Double_t kMassX   = 120.0; // X mass
static const Double_t kMassZ   = 91.19;	// Z mass
static const Double_t kEcm     = 500.;  // Ecm
static const Double_t kEa      = (kEcm/2.)*(1.-kMassX/kEcm)*(1.+kMassX/kEcm);
static const Double_t kSigmaMx =   2.0;	// X mass resolution
static const Double_t kSigmaEa =   3.3;	// gamma energy resolution

Bool_t gDEBUG = kFALSE;

      Int_t  gNgoods  = 0;
const Int_t  kMaxCuts = 100;
      stringstream gCutName[kMaxCuts];   
      TH1D  *hStat = 0;

// -------------------------------------------------------------------------
//  ------------------
//  AXAnalysis Class
//  ------------------

ClassImp(AXAnalysis)

// Constructor
AXAnalysis::AXAnalysis(const Char_t *name, const Char_t *title)
            : JSFModule(name, title)
{
}

// Destructor
AXAnalysis::~AXAnalysis()
{
}

// -------------------------------------------------------------------------
Bool_t AXAnalysis::Initialize()
{

  xEtrack   =   0.10;	// track energy
  xEvis     = 100.00;	// Minimum visible energy
  xPt       = 999.00;	// Pt maximum
  xPl       = 999.00;	// Pl maximum
  xYcut     =  0.004;	// y_cut to force the event to 4 jets
  xNjets    =      3;	// No. of Jets
  xEjet     =   1.00;	// E_jet minimum
  xCosjet   =   0.99;	// |cos(theta_j)| maximum
  xCosrsax  =   0.99;	// |cos(theta_Z)| and |cos(theta_H)| maximum
  xM2j      =   9.00;	// |m_jj-m_S| maximum , S = Z,H

  return 0;
}

// -------------------------------------------------------------------------
Bool_t AXAnalysis::Process(Int_t ev)
{
  // Local copies of AXAnalysisBuf data members.

  Int_t		fNtracks;	// track multiplicity
  Double_t	fEvis;		// visible energy
  Double_t	fPt;		// Pt
  Double_t	fPl;		// Pl
  Double_t	fYcut;		// y_cut to force the event to 4 jets
  Int_t		fNjets;		// jet multiplicity

  // Remember the previous directory.
  
  TDirectory *last = gDirectory;
  gFile->cd("/");

  // Analysis starts here.

  Float_t selid = -0.5;
  if (!hStat) hStat = new TH1D("hStat", "", 20, 0., 20.);
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "No Cuts" << ends;
  }

  // Get event buffer and make combined tracks accessible.

  JSFSIMDST     *sds	= static_cast<JSFSIMDST*>(gJSF->FindModule("JSFSIMDST"));
  JSFSIMDSTBuf	*evt	= static_cast<JSFSIMDSTBuf*>(sds->EventBuf());

  // Select good tracks

  ANL4DVector qsum;
  TObjArray tracks(1000);
  tracks.SetOwner();
  fNtracks = 0;
#ifndef __USEGENERATORINFO__
  Int_t		ntrks	= evt->GetNLTKCLTracks();	// No. of tracks
  TObjArray	*trks	= evt->GetLTKCLTracks();	// Combined tracks
  for (Int_t i = 0; i < ntrks; i++) {
    JSFLTKCLTrack *t = static_cast<JSFLTKCLTrack*>(trks->UncheckedAt(i));
    if (t->GetE() > xEtrack) {
      ANLTrack *qt = new ANLTrack(t);
      tracks.Add(qt);		// track 4-momentum
      qsum += *qt;		// total 4-mometum
      fNtracks++;		// *qt stays.
    }
  }
#else
  Int_t         ntrks = evt->GetNGeneratorParticles();
  TClonesArray *trks   = evt->GetGeneratorParticles();
  for (Int_t i = 0; i < ntrks; i++) {
    JSFGeneratorParticle *t = static_cast<JSFGeneratorParticle *>(trks->UncheckedAt(i));
    if (t->GetE() > xEtrack) {
      ANL4DVector *qt = new ANL4DVector(t->GetPV());
      tracks.Add(qt);		// track 4-momentum
      qsum += *qt;		// total 4-mometum
      fNtracks++;		// *qt stays.
    }
  }
#endif
  if (gDEBUG) cerr << "Ntracks = " << fNtracks << endl;

  // Cut on No. of tracks.
  
  if (fNtracks < xNtracks) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Ntracks >= " << xNtracks << ends;
  }
  
  fEvis = qsum(0);	// E_vis
  fPt	= qsum.GetPt(); // P_t
  fPl	= qsum(3);	// P_l

  if (gDEBUG) cerr << "Evis = " << fEvis << " Pt = "
	<< fPt << " Pl = " << fPl << endl;

  // Cut on Evis.

  if (fEvis < xEvis) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0)  {
    gCutName[(Int_t)selid] << "Evis >= " << xEvis << ends;
  }

  // Cut on Pt.
  
  if (fPt > xPt) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Pt <= " << xPt << ends;
  }

  // Cut on Pl.
  
  if (TMath::Abs(fPl) > xPl) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|Pl| <= " << xPl << ends;
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

  if (fNjets < xNjets) { return kFALSE;}
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Njets >= " << xNjets << " for Ycut = " << xYcut << ends;
  }

  // Now force the event to be xNjets.

  jclust.ForceNJets(xNjets);
  fNjets = jclust.GetNjets();
  fYcut  = jclust.GetYcut();

  if(gDEBUG) cerr << "Ycut = " << fYcut << " Njets = " << fNjets << endl;

  // Make sure the No. of Jets is xNjets.
  
  if (fNjets != xNjets) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Njets = " << xNjets << ends;
  }

  TObjArray &jets = jclust.GetJets();
  TIter nextjet(&jets);
  ANLJet *jetp;
  Double_t ejetmin = 999999.;
  Double_t cosjmax = 0.;
  while ((jetp = static_cast<ANLJet*>(nextjet()))) {
    ANLJet &jet = *jetp;
    if (gDEBUG) jet.DebugPrint();
    Double_t ejet = jet()(0);
    if (ejet < ejetmin) ejetmin = ejet;
    Double_t cosj = jet.CosTheta();
    if (TMath::Abs(cosj) > TMath::Abs(cosjmax)) cosjmax = cosj;
  }

  // Cut on Ejet_min.

  if (ejetmin < xEjet) {  return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Ejet > " << xEjet << ends;
  }

  // Cut on |cos(theta_j)|_max.

  if (TMath::Abs(cosjmax) > xCosjet ) {  return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_j)| <= " << xCosjet << ends;
  }

  // Find AX candidates in given mass window.
  
  TObjArray        solutions(10);
  solutions.SetOwner();
  ANLPairCombiner  xcandidates(jets,jets);
  ANLPair         *xp;
  while ((xp = static_cast<ANLPair *>(xcandidates()))) {
    ANLPair &x = *xp;
    Double_t xmass = x().GetMass();
    if (TMath::Abs(xmass - kMassX) > xM2j) continue;
    x.LockChildren();
    TIter acandidates(&jets);
    ANLJet *ap;
    while ((ap = static_cast<ANLJet*>(acandidates()))) {
      ANLJet &a = *ap;
      if (a.IsLocked()) continue;
      if (gDEBUG) {
        cerr << " M_X = " << xmass << endl;
      }
      Double_t chi2 = TMath::Power((xmass - kMassX)/kSigmaMx,2.)
                      + TMath::Power((a.E() - kEa)/kSigmaEa,2);
      solutions.Add(new ANLPair(xp,ap,chi2));
    }
    x.UnlockChildren();
  }

  // Cut on number of solutions.
 
  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|m_jj - m_X| <= " << xM2j << ends;
  }

  // Cut on cos(theta_AX).

  TIter nextsol(&solutions);
  ANLPair *sol;
  while ((sol = static_cast<ANLPair*>(nextsol()))) {
    ANLPair    &x = *(ANLPair *)(*sol)[0];
    ANLJet     &a = *(ANLJet  *)(*sol)[1];
    Double_t cosx = x.CosTheta();
    Double_t cosa = a.CosTheta();
    if (TMath::Abs(cosx) > xCosrsax || TMath::Abs(cosa) > xCosrsax) {
      solutions.Remove(sol);
      delete sol;
    }
  }
  if (!solutions.GetEntries()) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_x)| and |cos(theta_a)| <= " << xCosrsax << ends;
  }

  // End of event selection
  
  gNgoods++;

  cerr << "------------------------------------------" << endl
       << "Event " << gJSF->GetEventNumber()
       << ": Number of solutions = " << solutions.GetEntries() << endl
       << "------------------------------------------" << endl;
  
  // Sort the solutions in the ascending order of chi2 values.

  solutions.Sort();

  // Hists and plots for selected events.

  if (gDEBUG) {
    Int_t nj = 0;
    nextjet.Reset();
    while ((jetp = static_cast<ANLJet*>(nextjet()))) {
      cerr << "------" << endl
           << "Jet " << ++nj << endl
           << "------" << endl;
      jetp->DebugPrint();
    }
  }

  //--
  // Save the best solution.
  //--
  static TNtupleD *hEvt = 0;
  if (!hEvt) {
    stringstream sout;
    sout << "chi2:evis:pt:pl:cosx:phix:cosa:phia:mx:ex"            << ":"
         << "ea1:pa1x:pa1y:pa1z:ea2:pa2x:pa2y:pa2z:ea:pax:pay:paz" << ":"
         << "cosa1h:phia1h:cosa2h:phia2h"                          << ends;
    hEvt = new TNtupleD("hEvt","",sout.str().data());
  }
  nextsol.Reset();
  Int_t nsols = 0;
  while ((sol = static_cast<ANLPair*>(nextsol()))) {
    if (nsols++) break;
    ANLPair &x = *static_cast<ANLPair *>((*sol)[0]);
    ANLJet  &a = *static_cast<ANLJet  *>((*sol)[1]);
    Double_t chi2  = sol->GetQuality();
    //--
    // Calculate helicity angles.
    //--
    TVector3    ez    = TVector3(0., 0., 1.);
    TVector3    exz   = x.Vect().Unit();
    TVector3    exx   = exz.Cross(ez).Unit();
    TVector3    exy   = exz.Cross(exx).Unit();
  
    const ANL4DVector &gam1 = *static_cast<ANL4DVector *>(x[0]);
    TVector3    bstx  = TVector3(0., 0., x.Vect().Mag()/x.E());
    ANL4DVector k1h   = ANL4DVector(gam1.E(), gam1.Vect()*exx,
                                              gam1.Vect()*exy,
                                              gam1.Vect()*exz);
    k1h.Boost(-bstx);
    Double_t    csa1h = k1h.CosTheta();
    Double_t    fia1h = k1h.Phi();
  
    const ANL4DVector &gam2 = *static_cast<ANL4DVector *>(x[1]);
    ANL4DVector k2h   = ANL4DVector(gam2.E(), gam2.Vect()*exx,
                                              gam2.Vect()*exy,
                                              gam2.Vect()*exz);
    k2h.Boost(-bstx);
    Double_t    csa2h = k2h.CosTheta();
    Double_t    fia2h = k2h.Phi();
#if 0
    cerr << " k1h = (" << k1h.E () << ", "
                       << k1h.Px() << ", "
                       << k1h.Py() << ", "
                       << k1h.Pz() << ") " << endl;
    cerr << " k2h = (" << k2h.E () << ", "
                       << k2h.Px() << ", "
                       << k2h.Py() << ", "
                       << k2h.Pz() << ") " << endl;
#endif
    //--
    // Save them to the Ntuple.
    //--
    Double_t data[100];
    data[ 0] = chi2;
    data[ 1] = fEvis;
    data[ 2] = fPt;
    data[ 3] = fPl;
    data[ 4] = x.CosTheta();
    data[ 5] = x.Phi();
    data[ 6] = a.CosTheta();
    data[ 7] = a.Phi();
    data[ 8] = x.GetMass();
    data[ 9] = x.E();
    data[10] = gam1.E();
    data[11] = gam1.Px();
    data[12] = gam1.Py();
    data[13] = gam1.Pz();
    data[14] = gam2.E();
    data[15] = gam2.Px();
    data[16] = gam2.Py();
    data[17] = gam2.Pz();
    data[18] = a.E();
    data[19] = a.Px();
    data[20] = a.Py();
    data[21] = a.Pz();
    data[22] = csa1h;
    data[23] = fia1h;
    data[24] = csa2h;
    data[25] = fia2h;
    hEvt->Fill(data);
  }

  last->cd();
  return kTRUE;
}


// *****Terminate()***** //
Bool_t AXAnalysis::Terminate()
{
  // This function is called at the end of job.
  //--
  // Print out cut statistics
  //--
  cerr << endl
       << "  =============" << endl
       << "   Cut Summary " << endl
       << "  =============" << endl
       << endl
       << "  -----------------------------------------------------------" << endl
       << "   ID   No.Events    Cut Description                         " << endl
       << "  -----------------------------------------------------------" << endl;
  for (int id=0; id<kMaxCuts && gCutName[id].str().data()[0]; id++) {
     cerr << "  " << setw( 3) << id 
          << "  " << setw(10) << static_cast<int>(hStat->GetBinContent(id+1))
          << "  : " << gCutName[id].str().data() << endl;
  }
  cerr << "  -----------------------------------------------------------" << endl;
  //--
  // That's it!
  //--
  return 0;
}
// *****End of Terminate()***** //
