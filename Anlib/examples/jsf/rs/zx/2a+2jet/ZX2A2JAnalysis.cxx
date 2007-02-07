//***************************************************************************
//* ========================
//*  ZX2A2JAnalysis Classes
//* ========================
//* (Description)
//*   A user analysis class for JLC analyses.
//*   This reads and analyzes MC e+e- -> ZX data.
//* (Requires)
//*   library Anlib     (in LEDA and JSF)
//*   library RSZXStudy (in physsim)
//* (Provides)
//*   class ZX2A2JAnalysis
//*   class ZX2A2JAnalysisBuf
//* (Usage)
//*	...
//* (Update Record)
//*   2007/02/03 K.Fujii          Original version.
//***************************************************************************

#include "ZX2A2JAnalysis.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TH1D.h"
#include "JSFSIMDST.h"
#include "Anlib.h"
#include "ANLTrack.h"
#include <sstream>

static const Double_t kMassX   = 120.0; // H mass
static const Double_t kMassZ   = 91.19;	// Z mass
static const Double_t kSigmaMx =   4.0;	// H mass resolution
static const Double_t kSigmaMz =   4.0;	// Z mass resolution

#if 1
      Bool_t gDEBUG  = kFALSE;
#else
      Bool_t gDEBUG  = kTRUE;
#endif
      Int_t  gNgoods  = 0;
const Int_t  kMaxCuts = 100;
      Char_t gCutName[kMaxCuts][100];
      TH1D  *hStat = 0;

//___________________________________________________________________________
// ----------------------------
//  ZX2A2JAnalysis Class
// ----------------------------

ClassImp(ZX2A2JAnalysis)

// Constructor
ZX2A2JAnalysis::ZX2A2JAnalysis(const Char_t *name, 
                               const Char_t *title)
              : JSFModule(name, title),
                fEtrackCut( 0.10),
                fEgammaCut(10.00),
                fCosCone  ( 0.94),
                fEconeCut ( 5.00),
                fYcutCut  (0.004),
                fM2aCut   (9.000)
{
  SetBufferSize(2000);	// buffer size for event data
  cout << "ZX2A2JAnalysis is created... fEventBuf is "
       << (Int_t)fEventBuf << endl;
}

// Destructor
ZX2A2JAnalysis::~ZX2A2JAnalysis()
{
  cout << "ZX2A2JAnalysisBuf will be deleted... fEventBuf is "
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}

// Functions
void ZX2A2JAnalysis::CleanUp(TObjArray *objs)
{
  TIter next(objs);
  TObject *obj;
  while ( (obj = next()) ) {
    objs->Remove(obj);
    delete obj;
  }
}

Bool_t ZX2A2JAnalysis::Process(Int_t ev)
{
  using namespace std;
  //--
  // Remember the previous directory.
  //--
  TDirectory *last = gDirectory;
  gFile->cd("/");
  if (!hStat) hStat = new TH1D("hStat","Cut Statistics", 20, 0., 20.);
  Char_t msg[60];

  static TNtupleD *hEvt = 0;
  if (!hEvt) {
    stringstream tupstr;
    tupstr << "ntracks:evis:pt:pl:ngams:njets:ycut:chi2"                << ":"
           << "ea1:pa1x:pa1y:pa1z:ea2:pa2x:pa2y:pa2z"                   << ":"
           << "npj1:ej1:pj1x:pj1y:pj1z:npj2:ej2:pj2x:pj2y:pj2z"         << ":"
           << "zmass:xmass:acop"                                        << ":"
           << "cosa1h:phia1h:cosa2h:phia2h:cosj1h:phij1h:cosj2h:phij2h" << ":"
           << "econe1:econe2"                                           << ends;

    hEvt = new TNtupleD("hEvt", "", tupstr.str().data());
  }

  //--
  // Analysis starts here.
  //--
  Double_t selid = -0.5;
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    sprintf(msg,"No Cuts");
    strcpy(&gCutName[(Int_t)selid][0],msg);
  }

  //--
  // Get event buffer and make combined tracks accessible.
  //--
  JSFSIMDST     *sdsp = static_cast<JSFSIMDST*>(gJSF->FindModule("JSFSIMDST"));
  JSFSIMDSTBuf	*evtp = static_cast<JSFSIMDSTBuf*>(sdsp->EventBuf());

  Int_t	     ntrks = evtp->GetNLTKCLTracks();	// No. of tracks
  TObjArray *trks  = evtp->GetLTKCLTracks ();	// Combined tracks

  //--
  // Select good tracks
  //--
  Int_t    ntracks = 0;
  ANL4DVector qsum;
  TObjArray tracks(1000);
  for (Int_t i=0; i<ntrks; i++) {
    JSFLTKCLTrack *tp = (JSFLTKCLTrack*)trks->UncheckedAt(i);
    if (tp->GetE() > fEtrackCut) {
      ANLTrack *qtp = new ANLTrack(tp);
      tracks.Add(qtp);		// track 4-momentum
      qsum += *qtp;		// total 4-mometum
      ntracks++;		// *qt stays.
    }
  }
  if (gDEBUG) cerr << "Ntracks = " << ntracks << endl;

  //--
  // Cut on No. of tracks.
  //--
  if ( ntracks < 4 ) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    sprintf(msg,"Ntracks < 4");
    strcpy(&gCutName[(Int_t)selid][0],msg);
  }
  
  Double_t evis  = qsum(0);	 // E_vis
  Double_t pt    = qsum.GetPt(); // P_t
  Double_t pl    = qsum(3);	 // P_l

  if (gDEBUG) cerr << " Evis = " << evis 
                   << " Pt = " << pt << " pl = " << pl << endl;
  //--
  // Find isolated gammas.
  //--
  Double_t eacone[100];
  Int_t    ngams   = 0;
  TObjArray gammas(20);
  TIter nexttrk (&tracks);
  ANLTrack *trkp;
  while ((trkp = static_cast<ANLTrack *>(nexttrk()))) {
    ANLTrack &trk = *trkp;
    if (!trk.IsGamma()) continue;
    Double_t egamma = trk.E();
    if (egamma < fEgammaCut) continue;
    Double_t econe = trk.GetConeEnergy(fCosCone, &tracks);
    if (gDEBUG) cerr << " econe  = " << econe 
                     << " egamma = " << egamma << endl;
    if (econe <= fEconeCut) {
      gammas.Add(trkp);
      eacone[ngams] = econe;
      ngams++;
    }
  }

  if (gDEBUG) cerr << " ngams = " << ngams << endl;

  //--
  // Require more than two isolated gammas.
  //--
  if (ngams < 2) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    sprintf(msg,"Ngammas >= 2");
    strcpy(&gCutName[(Int_t)selid][0],msg);
  }

  //--
  // Lock the best gamma combination as an X candidate.
  //--
  TObjArray solutions(10);
  ANLPairCombiner xcandidates(gammas,gammas);
  ANLPair *xp;
  while ((xp = static_cast<ANLPair *>(xcandidates()))) {
    ANLPair &x = *xp;
    Double_t xmass = x().GetMass();
    if (TMath::Abs(xmass - kMassX) > fM2aCut) continue;
    Double_t chi2  = TMath::Power((xmass - kMassX)/kSigmaMx,2.);
    solutions.Add(new ANLPair(x[0],x[1],chi2));    
  }

  if (gDEBUG) cerr << " solutions = " << solutions.GetEntries() << endl;

  //--
  // Require at least one solution.
  //--
  if (solutions.GetEntries() < 1) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    sprintf(msg,"Solutions > 0");
    strcpy(&gCutName[(Int_t)selid][0],msg);
  }

  solutions.Sort();  
  ANLPair &x = *static_cast<ANLPair *>(solutions.At(0));
  x.LockChildren();
  ANLTrack &gam1 = *static_cast<ANLTrack *>(x[0]);
  ANLTrack &gam2 = *static_cast<ANLTrack *>(x[1]);

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
  if (njets < 2) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    sprintf(msg,"Njets >= 2 for Ycut = %g",fYcutCut);
    strcpy(&gCutName[(Int_t)selid][0],msg);
  }

  //--
  // Now force the event to be 2-jet.
  //--
  jclust.ForceNJets(2);
  njets = jclust.GetNjets();
  ycut  = jclust.GetYcut();

  if (gDEBUG) cerr << "Ycut = " << ycut << " Njets = " << njets << endl;

  //--
  // Make sure the No. of Jets is 2.
  //--
  if (njets != 2) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    sprintf(msg,"Njets = 2");
    strcpy(&gCutName[(Int_t)selid][0],msg);
  }

  //--
  // Find XH candidates.
  //--

  TObjArray &jets = jclust.GetJets();
  ANLJet    &j1   = *static_cast<ANLJet *>(jets[0]);
  ANLJet    &j2   = *static_cast<ANLJet *>(jets[1]);
  ANLPair z(&j1, &j2);

  Double_t zmass = z().GetMass();
  Double_t xmass = x().GetMass();
  Double_t chi2  = x.GetQuality() + TMath::Power((zmass - kMassZ)/kSigmaMz,2);
#if 1
  if (xmass < 100.) {
    cerr << " ============================ " << endl;
    cerr << " zmass = " << zmass << " xmass = " << xmass << endl;
    cerr << " &gam1 = " << &gam1 << " &gam2 = " << &gam2 << endl;
    nexttrk.Reset();
    Int_t itk = 0;
    ANLTrack *trkp;
    while ((trkp = static_cast<ANLTrack *>(nexttrk()))) {
      cerr << " trk = "  << ++itk
           << " ("       << trkp       << ") "
           << " pv  = (" << trkp->E () << ", "
                         << trkp->Px() << ", "
                         << trkp->Py() << ", "
                         << trkp->Pz() << ") " << endl;
    }
    cerr << " ============================ " << endl;
  }
#endif

  //--
  // Calculate helicity angles.
  //--
  TVector3    ez    = TVector3(0., 0., 1.);
  TVector3    exz   = x.Vect().Unit();
  TVector3    exx   = exz.Cross(ez).Unit();
  TVector3    exy   = exz.Cross(exx);

  TVector3    bstx  = TVector3(0., 0., x.Vect().Mag()/x.E());
  ANL4DVector k1h   = ANL4DVector(gam1.E(), gam1.Vect()*exx,
                                            gam1.Vect()*exy,
                                            gam1.Vect()*exz);
  k1h.Boost(-bstx);
  Double_t    csa1h = k1h.CosTheta();
  Double_t    fia1h = k1h.Phi();

  ANL4DVector k2h   = ANL4DVector(gam2.E(), gam2.Vect()*exx,
                                            gam2.Vect()*exy,
                                            gam2.Vect()*exz);
  k2h.Boost(-bstx);
  Double_t    csa2h = k2h.CosTheta();
  Double_t    fia2h = k2h.Phi();

  TVector3    ezz   = z.Vect().Unit();
  TVector3    ezx   = ezz.Cross(ez).Unit();
  TVector3    ezy   = ezz.Cross(ezx);
  TVector3    bstz  = TVector3(0., 0., z.Vect().Mag()/z.E());
  ANL4DVector p1h   = ANL4DVector(j1.E(), j1.Vect()*ezx,
                                          j1.Vect()*ezy,
                                          j1.Vect()*ezz);
  p1h.Boost(-bstz);
  Double_t    csj1h = p1h.CosTheta();
  Double_t    fij1h = p1h.Phi();

  ANL4DVector p2h   = ANL4DVector(j2.E(), j2.Vect()*ezx,
                                          j2.Vect()*ezy,
                                          j2.Vect()*ezz);
  p2h.Boost(-bstz);
  Double_t    csj2h = p2h.CosTheta();
  Double_t    fij2h = p2h.Phi();

  //--
  // Fill up Ntuple.
  //--
  Double_t data[100];
  data[ 0] = ntracks;
  data[ 1] = evis;
  data[ 2] = pt;
  data[ 3] = pl;
  data[ 4] = ngams;
  data[ 5] = njets;
  data[ 6] = ycut;
  data[ 7] = chi2;
  data[ 8] = gam1(0);
  data[ 9] = gam1(1);
  data[10] = gam1(2);
  data[11] = gam1(3);
  data[12] = gam2(0);
  data[13] = gam2(1);
  data[14] = gam2(2);
  data[15] = gam2(3);
  data[16] = j1.GetNparticles();
  data[17] = j1()(0);
  data[18] = j1()(1);
  data[19] = j1()(2);
  data[20] = j1()(3);
  data[21] = j2.GetNparticles();
  data[22] = j2()(0);
  data[23] = j2()(1);
  data[24] = j2()(2);
  data[25] = j2()(3);
  data[26] = z().GetMass();
  data[27] = x().GetMass();
  data[28] = z.Acop(x);
  data[29] = csa1h;
  data[30] = fia1h;
  data[31] = csa2h;
  data[32] = fia2h;
  data[33] = csj1h;
  data[34] = fij1h;
  data[35] = csj2h;
  data[36] = fij2h;
  data[37] = eacone[0];
  data[38] = eacone[1];

  hEvt->Fill(data);

  //--
  // End of event selection
  //--
  if (gNgoods == 0) {
    selid++;
    sprintf(msg,"END");
    strcpy(&gCutName[(Int_t)selid][0],msg);
  }
  gNgoods++;

  cerr << "------------------------------------------" << endl
       << "Event " << gJSF->GetEventNumber()           << endl
       << "------------------------------------------" << endl;
  
  //--
  // Clean up
  //--
  CleanUp(&tracks);
  CleanUp(&solutions);

  last->cd();
  return kTRUE;
}

Bool_t ZX2A2JAnalysis::Terminate()
{
  // This function is called at the end of job.
  cout << endl;
  cout << "  =============" << endl;
  cout << "   Cut Summary " << endl;
  cout << "  =============" << endl;
  cout << endl;
  cout << "  ---------------------------------------------------------" << endl;   
  cout << "  ID   No.Events     Cut Description"                        << endl;
  cout << "  ---------------------------------------------------------" << endl;
  Int_t i;
  for ( i = 0; strncmp(&gCutName[i][0],"END",4) && i < kMaxCuts ; i++ ) {
    printf("  %3d  %10d  : %s\n",i,(int)hStat->GetBinContent(i+1),&gCutName[i][0]);
  } 
  cout << "  ---------------------------------------------------------";
  return 0;
}
