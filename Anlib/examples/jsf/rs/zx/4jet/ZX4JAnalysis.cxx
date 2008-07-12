//***************************************************************************
//* ========================
//*  ZX4JAnalysis Classes
//* ========================
//* (Description)
//*   A user analysis class for JLC analyses.
//*   This reads and analyzes MC e+e- -> ZX data.
//* (Requires)
//*   library Anlib     (in LEDA and JSF)
//*   library RSZXStudy (in physsim)
//* (Provides)
//*   class ZX4JAnalysis
//*   class ZX4JAnalysisBuf
//* (Usage)
//*	...
//* (Update Record)
//*   2007/02/22 K.Fujii          Original version.
//***************************************************************************

#include "ZX4JAnalysis.h"
#include "TROOT.h"
#include "TFile.h"
#include "TNtupleD.h"
#include "TH1D.h"
#include "JSFSIMDST.h"
#include "Anlib.h"
#include "ANLTrack.h"
#include <sstream>
#include <iomanip>

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
      stringstream gCutName[kMaxCuts];   
      TH1D  *hStat = 0;

//___________________________________________________________________________
// ----------------------------
//  ZX4JAnalysis Class
// ----------------------------

ClassImp(ZX4JAnalysis)

// Constructor
ZX4JAnalysis::ZX4JAnalysis(const Char_t *name, 
                               const Char_t *title)
              : JSFModule(name, title),
                fEtrackCut( 0.10),
                fYcutCut  (0.004),
                fM2jCut   (12.000)
{
  SetBufferSize(2000);	// buffer size for event data
  cout << "ZX4JAnalysis is created... fEventBuf is "
       << (Int_t)fEventBuf << endl;
}

// Destructor
ZX4JAnalysis::~ZX4JAnalysis()
{
  cout << "ZX4JAnalysisBuf will be deleted... fEventBuf is "
       << (Int_t)fEventBuf << endl;
  delete fEventBuf; fEventBuf = 0;
}

// Functions
void ZX4JAnalysis::CleanUp(TObjArray *objs)
{
  TIter next(objs);
  TObject *obj;
  while ( (obj = next()) ) {
    objs->Remove(obj);
    delete obj;
  }
}

Bool_t ZX4JAnalysis::Process(Int_t ev)
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
    tupstr << "ntracks:evis:pt:pl:ycut:chi2"                            << ":"
           << "npg1:eg1:pg1x:pg1y:pg1z:npg2:eg2:pg2x:pg2y:pg2z"         << ":"
           << "npj1:ej1:pj1x:pj1y:pj1z:npj2:ej2:pj2x:pj2y:pj2z"         << ":"
           << "zmass:xmass:acop"                                        << ":"
           << "cosg1h:phig1h:cosg2h:phig2h:cosj1h:phij1h:cosj2h:phij2h" << ":"
           << ends; 

    hEvt = new TNtupleD("hEvt", "", tupstr.str().data());
  }

  //--
  // Analysis starts here.
  //--
  Double_t selid = -0.5;
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "No Cuts" << ends;
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
    gCutName[(Int_t)selid] << "Ntracks >= 4" << ends;
  }
  
  Double_t evis  = qsum(0);	 // E_vis
  Double_t pt    = qsum.GetPt(); // P_t
  Double_t pl    = qsum(3);	 // P_l

  if (gDEBUG) cerr << " Evis = " << evis 
                   << " Pt = " << pt << " pl = " << pl << endl;
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
  if (njets < 4) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Njets >= 4 for Ycut = " << fYcutCut << ends;
  }

  //--
  // Now force the event to be 4-jet.
  //--
  jclust.ForceNJets(4);
  njets = jclust.GetNjets();
  ycut  = jclust.GetYcut();

  if (gDEBUG) cerr << "Ycut = " << ycut << " Njets = " << njets << endl;

  //--
  // Make sure the No. of Jets is 4.
  //--
  if (njets != 4) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Njets = 4" << ends;
  }

  //--
  // Find XH candidates.
  //--

  TObjArray solutions(10);
  TObjArray &jets = jclust.GetJets();
  ANLPairCombiner zcandidates(jets,jets);
  ANLPair *xp, *zp;
  while ((zp = (ANLPair *)zcandidates())) {
    ANLPair &z = *zp;
    Double_t zmass = z().GetMass();
    if (TMath::Abs(zmass - kMassZ) > fM2jCut) continue;
    z.LockChildren();
    ANLPairCombiner xcandidates(zcandidates);
    xcandidates.Reset();
    while ((xp = (ANLPair *)xcandidates())) {
      ANLPair &x = *xp;
      if (x.IsLocked()) continue;
      Double_t xmass = x().GetMass();
      if (TMath::Abs(xmass - kMassX) > fM2jCut) continue;
      if (gDEBUG) {
        cerr << " M_z = " << zmass << " M_x = " << xmass << endl;
      }
      Double_t chi2 = TMath::Power((zmass - kMassZ)/kSigmaMz,2.)
                    + TMath::Power((xmass - kMassX)/kSigmaMx,2.);
      solutions.Add(new ANLPair(zp,xp,chi2));
    }
    z.UnlockChildren();
  }

  if (gDEBUG) cerr << " solutions = " << solutions.GetEntries() << endl;

  //--
  // Require at least one solution.
  //--
  if (solutions.GetEntries() < 1) { CleanUp(&tracks); return kFALSE; }
  hStat->Fill(++selid);
  if (gNgoods == 0) {
    gCutName[(Int_t)selid] << "Solutions > 0" << ends;
  }

  //--
  // Select the best solution.
  //--
  solutions.Sort();  
  ANLPair &sol = *static_cast<ANLPair *>(solutions.At(0));

  ANLPair &z   = *static_cast<ANLPair *>(sol[0]);
  ANLJet  &j1  = *static_cast<ANLJet  *>(z[0]);
  ANLJet  &j2  = *static_cast<ANLJet  *>(z[1]);
  ANLPair &x   = *static_cast<ANLPair *>(sol[1]);
  ANLJet  &g1  = *static_cast<ANLJet  *>(x[0]);
  ANLJet  &g2  = *static_cast<ANLJet  *>(x[1]);

  //--
  // Calculate helicity angles.
  //--
  TVector3    ez    = TVector3(0., 0., 1.);
  TVector3    exz   = x.Vect().Unit();
  TVector3    exx   = exz.Cross(ez).Unit();
  TVector3    exy   = exz.Cross(exx);

  TVector3    bstx  = TVector3(0., 0., x.Vect().Mag()/x.E());
  ANL4DVector k1h   = ANL4DVector(g1.E(), g1.Vect()*exx,
                                          g1.Vect()*exy,
                                          g1.Vect()*exz);
  k1h.Boost(-bstx);
  Double_t    csg1h = k1h.CosTheta();
  Double_t    fig1h = k1h.Phi();

  ANL4DVector k2h   = ANL4DVector(g2.E(), g2.Vect()*exx,
                                          g2.Vect()*exy,
                                          g2.Vect()*exz);
  k2h.Boost(-bstx);
  Double_t    csg2h = k2h.CosTheta();
  Double_t    fig2h = k2h.Phi();

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
  data[ 4] = ycut;
  data[ 5] = sol.GetQuality();
  data[ 6] = g1.GetNparticles();
  data[ 7] = g1()(0);
  data[ 8] = g1()(1);
  data[ 9] = g1()(2);
  data[10] = g1()(3);
  data[11] = g2.GetNparticles();
  data[12] = g2()(0);
  data[13] = g2()(1);
  data[14] = g2()(2);
  data[15] = g2()(3);
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
  data[29] = csg1h;
  data[30] = fig1h;
  data[31] = csg2h;
  data[32] = fig2h;
  data[33] = csj1h;
  data[34] = fij1h;
  data[35] = csj2h;
  data[36] = fij2h;

  hEvt->Fill(data);

  //--
  // End of event selection
  //--
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

Bool_t ZX4JAnalysis::Terminate()
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
