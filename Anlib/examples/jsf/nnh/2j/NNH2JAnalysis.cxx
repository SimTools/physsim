//***************************************************************************
//*  =====================
//*  NNH2JAnalysis Classes
//*  =====================
//*
//*  (Description)
//*	A user analysis class for JLC analyses.
//*	This reads and analyzes MC e+e- -> NNH data.
//*  (Requires)
//*	library Anlib in LEDA
//*	library NNHStudy++ in physsim
//*  (Provides)
//*	class NNH2JAnalysis
//*  (Usage)
//*	...
//*  (Update Record)
//*    2010/04/28 K.Fujii	Original version.
//***************************************************************************

#include "NNH2JAnalysis.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TH1D.h"
#include "JSFSteer.h"
#include "JSFSIMDST.h"
#include "Anlib.h"
#include "ANLTrack.h"
#include "NNHSpring.h"

#include <sstream>
#include <iomanip>

using namespace std;

static const Double_t kMassH   = 120.0; // H mass
static const Double_t kMassZ   = 91.19;	// Z mass
static const Double_t kSigmaMh =   4.0;	// H mass resolution
static const Double_t kSigmaMz =   4.0;	// Z mass resolution
static const Int_t    kZoneX   =     4;	// No. of X Zones in the Canvas
static const Int_t    kZoneY   =     4;	// No. of Y Zones in the Canvas

typedef enum { kElectron = 11, kMuon = 13 } EPID;

Bool_t gDEBUG = kFALSE;
Int_t  gNgoods  = 0;
const Int_t  kMaxCuts = 100;
      stringstream gCutName[kMaxCuts];
      TH1D  *hStat = 0;

//--------------------------------------------------------------------
// NNH2JAnalysis Class
//--------------------------------------------------------------------

ClassImp(NNH2JAnalysis)

// Constructor
NNH2JAnalysis::NNH2JAnalysis(const Char_t *name,
		             const Char_t *title)
             : JSFModule(name, title),
               fNtracksCut(    5),   // No. of tracks
               fEtrackCut (  0.1),   // track energy
               fEvisLoCut (  0.0),   // Minimum visible energy
               fEvisHiCut (999.0),   // Maximum visible energy
               fPtCut     (  0.0),   // Pt minimum
               fPlCut     (999.0),   // Pl maximum
               fEleptonCut(  0.1),   // track energy
               fCosConeCut(999.0),   // |cos(theta_cone)|
               fEconeCut  (999.0),   // Econe maximum
               fNleptonCut(    0),   // No. of tracks
               fYcutCut   (0.004),   // y_cut to force the event to 4 jets
               fNjetsCut  (    2),   // No. of jets
               fEjetCut   (  5.0),   // E_jet minimum
               fCosjetCut (999.0),   // |cos(theta_jet)| maximum
               fCos2jCut  (999.0),   // |cos(theta_qq)| maximum
               fM2jCut    ( 18.0)    // |m_qq-m_H| maximum
{
  //--
  //  Read in Generator info.
  //--
  gJSF->GetInput()->cd("/conf/init");
  NNHBases *bsp = static_cast<NNHBases *>(gROOT->FindObject("NNHBases"));
  cerr << "------------------------------------" << endl
       << " Ecm = " << bsp->GetEcmInit() << " GeV" << endl
       << "------------------------------------" << endl;
  SetEcm(bsp->GetEcmInit()); 
}

// Destructor
NNH2JAnalysis::~NNH2JAnalysis()
{
  delete fEventBuf;
  fEventBuf = 0;
}


//####### Functions ###########


// ******Initialize()****** //
Bool_t NNH2JAnalysis::Initialize()
{
  return kTRUE;
}
// *****End of Initialize()***** //

// *****Process()***** //
Bool_t NNH2JAnalysis::Process(Int_t ev)
{
  using namespace std;

  //--
  // Remember the previous directory.
  //--
  TDirectory *last = gDirectory;
  gFile->cd("/");
  if (!hStat) hStat = new TH1D("hStat","Cut Statistics", 20, 0., 20.);

  static TNtupleD *hEvt = 0;
  if (!hEvt) {
    stringstream tupstr;
    tupstr << "nev:ecm:ntracks:evis:pt:pl:mqq:ej1:ej2:cosj1:cosj2"      << ends;
    hEvt = new TNtupleD("hEvt", "", tupstr.str().data());
  }

  // ---------------------
  // Analysis starts here.
  // ---------------------
  Double_t selid = -0.5;
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "No Cuts" << ends;
  }

  //--
  // Get event buffer and make combined tracks accessible.
  //--
  JSFSIMDST     *sdsp    = static_cast<JSFSIMDST *>(gJSF->FindModule("JSFSIMDST"));
  JSFSIMDSTBuf  *evtp    = static_cast<JSFSIMDSTBuf *>(sdsp->EventBuf());
  Int_t          ntrks   = evtp->GetNLTKCLTracks(); 	// No. of tracks 
  TObjArray     *trks    = evtp->GetLTKCLTracks(); 	// combined tracks

  //--
  // Select good tracks.
  //--
  ANL4DVector qsum;
  TObjArray tracks(1000);
  tracks.SetOwner();
  Int_t ntracks = 0;
  for (Int_t i = 0; i < ntrks; i++) {
    JSFLTKCLTrack *tp = static_cast<JSFLTKCLTrack*>(trks->UncheckedAt(i));
    if (tp->GetE() > fEtrackCut) {
      ANLTrack *qtp = new ANLTrack(tp);
      tracks.Add(qtp); 		// track 4-momentum
      qsum += *qtp;		// total 4-momentum
      ntracks++;
    }				// *qt stays.
  }
  if (gDEBUG) cerr << "Ntracks = " << ntracks << endl;

  //--
  // Cut on No. of tracks.
  //--
  if (ntracks < fNtracksCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Ntracks >= " << fNtracksCut << ends;
  }

  Double_t evis = qsum(0);		// E_vis
  Double_t pt   = qsum.GetPt();		// P_t
  Double_t pl   = qsum(3);		// P_l

  if (gDEBUG) cerr << "Evis = " << evis << " Pt = " << pt << " Pl = " << pl << endl;

  //--
  // Cut on Evis.
  //--
  if (evis < fEvisLoCut || evis > fEvisHiCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << fEvisHiCut << " <= Evis <= " << fEvisLoCut << ends;
  }

  //--
  // Cut on Pt.
  //--

  if (pt > fPtCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Pt <= " << fPtCut << ends;
  }

  //--
  // Cut on Pl.
  //--
  if (TMath::Abs(pl) > fPlCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|Pl| <= " << fPlCut << ends;
  }
  
  //--
  // Find Isolated Lepton.
  //--
  TObjArray lptracks(20);
  Int_t nlptracks = 0;
  TIter nexttrk (&tracks);
  ANLTrack *trkp;
  while ((trkp = static_cast<ANLTrack *>(nexttrk()))) {
    ANLTrack &trk = *trkp;
    if (!trk.IsLepton()) continue;
    Double_t elepton = trk.E();
    if (elepton < fEleptonCut) continue;
    Double_t econe = trk.GetConeEnergy(fCosConeCut, &tracks);
    if (econe <= fEconeCut) {
      lptracks.Add(trkp);
      nlptracks++;
      trk.Lock(); // lock isolated leptons
    }
  }

  //--
  // Require no isolated lepton.
  //--
  if (nlptracks != 0) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Nlptracks = 0" << ends;
  }

  //--
  // Find jets.
  //--
  ANLJadeEJetFinder jclust(fYcutCut);
  jclust.Initialize(tracks);
  jclust.FindJets();
  Double_t ycut  = jclust.GetYcut();
  Int_t    njets = jclust.GetNjets();

  if (gDEBUG) cerr << "Ycut = " << ycut << " Njets = " << njets << endl;
  //--
  // Cut on No. of jets.
  //--

  if (njets < fNjetsCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Njets >= " << fNjetsCut << " for Ycut = " << fYcutCut << ends;
  }

  //--
  // Now force the event to be fNjetsCut.
  //--
  jclust.ForceNJets(fNjetsCut);
  njets = jclust.GetNjets();
  ycut  = jclust.GetYcut();

  if (gDEBUG) cerr << "Ycut = " << ycut << " Njets = " << njets << endl;

  //--
  // Make sure that No. of jets is fNjetsCut.
  //--
  if (njets != fNjetsCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Njets = " << fNjetsCut << ends;
  }

  //--
  // Loop over jets and decide Ejet_min and |cos(theta_j)|_max.
  //--
  TObjArray &jets = jclust.GetJets();
  ANLJet &j1     = *static_cast<ANLJet *>(jets.At(0));
  Double_t ej1   = j1.E();
  Double_t cosj1 = j1.CosTheta();
  ANLJet &j2     = *static_cast<ANLJet *>(jets.At(1));
  Double_t ej2   = j2.E();
  Double_t cosj2 = j2.CosTheta();

  //--
  // Cut on Ejet_min.
  //--
  if (ej1 < fEjetCut || ej2 < fEjetCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Ejet >= " << fEjetCut << ends;
  }

  //--
  // Cut on |cos(theta_j)|_max.
  //--
  if (TMath::Abs(cosj1) > fCosjetCut || TMath::Abs(cosj2) > fCosjetCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_j)| <= " << fCosjetCut << ends;
  }

  //--
  // Find H candidates in the given mass window.
  //--
  Double_t hmass = qsum.GetMass();
  if (TMath::Abs(hmass - kMassH) > fM2jCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|m_qq - m_h| <= " << fM2jCut << ends;
  }

  //--
  // Cut on cos(theta_h).
  //--
  Double_t cosqq = qsum.CosTheta();
  if (TMath::Abs(cosqq) > fCos2jCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_qq)| <= " << fCos2jCut << ends;
  }

  // ----------------------
  // End of event selection
  // ----------------------
  gNgoods++;
  cerr << "--------------------------------------------------" << endl
       << "Event "                   << gJSF->GetEventNumber() << endl
       << "--------------------------------------------------" << endl;

    Double_t data[100];
    data[ 0] = gJSF->GetEventNumber();
    data[ 1] = GetEcm();
    data[ 2] = ntracks;
    data[ 3] = evis;
    data[ 4] = pt;
    data[ 5] = pl;
    data[ 6] = hmass;
    data[ 7] = ej1;
    data[ 8] = ej2;
    data[ 9] = cosj1;
    data[10] = cosj2;

    hEvt->Fill(data);

  last->cd();

  // ------------------------
  //  That's it.
  // ------------------------
  return kTRUE;
}
// *****End of Process()***** //


// *****Terminate()***** //
Bool_t NNH2JAnalysis::Terminate()
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
