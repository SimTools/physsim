//***************************************************************************
//*  =====================
//*  NNH2LAnalysis Classes
//*  =====================
//*
//*  (Description)
//*	A user analysis class for JLC analyses.
//*	This reads and analyzes MC e+e- -> NNH data.
//*  (Requires)
//*	library Anlib in LEDA
//*	library NNHStudy++ in physsim
//*  (Provides)
//*	class NNH2LAnalysis
//*  (Usage)
//*	...
//*  (Update Record)
//*    2010/04/28 K.Fujii	Original version.
//***************************************************************************

#include "NNH2LAnalysis.h"
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
// NNH2LAnalysis Class
//--------------------------------------------------------------------

ClassImp(NNH2LAnalysis)

// Constructor
NNH2LAnalysis::NNH2LAnalysis(const Char_t *name,
		             const Char_t *title)
             : JSFModule(name, title),
               fNtracksCut(    2),   // No. of tracks
               fEtrackCut (  0.1),   // track energy
               fEvisLoCut (  0.0),   // Minimum visible energy
               fEvisHiCut (999.0),   // Maximum visible energy
               fPtCut     (  0.0),   // Pt minimum
               fPlCut     (999.0),   // Pl maximum
               fEleptonCut(  0.1),   // track energy
               fCosConeCut(999.0),   // |cos(theta_cone)|
               fEconeCut  (999.0),   // Econe maximum
               fCosLepCut (999.0),   // |cos(theta_lepton)| maximum
               fM2lCut    ( 18.0)    // |m_ll-m_H| maximum
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
NNH2LAnalysis::~NNH2LAnalysis()
{
  delete fEventBuf;
  fEventBuf = 0;
}


//####### Functions ###########


// ******Initialize()****** //
Bool_t NNH2LAnalysis::Initialize()
{
  return kTRUE;
}
// *****End of Initialize()***** //

// *****Process()***** //
Bool_t NNH2LAnalysis::Process(Int_t ev)
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
    tupstr << "nev:ecm:ntracks:evis:pt:pl:mll:elm:elp:coslm:coslp"      << ends;
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
  Double_t Qsum     = 0.;
  Double_t lpcharge = 0.;
  Double_t elm      = 0.;
  Double_t elp      = 0.;
  Double_t coslm    = 0.;
  Double_t coslp    = 0.;
  while ((trkp = static_cast<ANLTrack *>(nexttrk()))) {
    ANLTrack &trk = *trkp;
    if (!trk.IsLepton()) continue;
    Double_t elepton = trk.E();
    if (elepton < fEleptonCut) continue;
    Double_t econe = trk.GetConeEnergy(fCosConeCut, &tracks);
    if (econe <= fEconeCut) {
      lptracks.Add(trkp);
      nlptracks++;
      trk.Lock();
      lpcharge = trk.GetCharge();
      Qsum += lpcharge;
      if (lpcharge > 0) {
         elp   = elepton;
         coslp = trk.CosTheta();
      } else {
         elm   = elepton;
         coslm = trk.CosTheta();
      }
    }
  }

  //--
  // Require two and only two isolated leptons.
  //--
  if (nlptracks != 2 || Qsum != 0.) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Nlptracks = 2 && Qsum = 0." << ends;
  }

  //--
  // Find H candidates in the given mass window.
  //--
  Double_t hmass = qsum.GetMass();
  if (TMath::Abs(hmass - kMassH) > fM2lCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|m_ll - m_h| <= " << fM2lCut << ends;
  }

  //--
  // Cut on cos(theta_h).
  //--
  Double_t cosll = qsum.CosTheta();
  if (TMath::Abs(cosll) > fCosLepCut) { return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "|cos(theta_ll)| <= " << fCosLepCut << ends;
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
    data[ 7] = elm;
    data[ 8] = elp;
    data[ 9] = coslm;
    data[10] = coslp;

    hEvt->Fill(data);

  last->cd();

  // ------------------------
  //  That's it.
  // ------------------------
  return kTRUE;
}
// *****End of Process()***** //


// *****Terminate()***** //
Bool_t NNH2LAnalysis::Terminate()
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
