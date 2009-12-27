//////////////////////////////////////////////////////////////////
//
//  JetAnalysis
//
//  Sample program for Jet Analysis
//
//$Id: 
//  
//////////////////////////////////////////////////////////////////
#include <iostream>
#include <iomanip>
#include "JSFSteer.h"
#include "JSFSIMDST.h"
#include "JetAnalysis.h"
#include "ANLJetFinder.h"
#include "ANLCheatedJetFinder.h"
#include "ANLVTXTagger.h"
#include "ANLGVTXTagger.h"
#include "FlavourGetter.h"
#include "ANLTrack.h"
#include "ANLEventShape.h"
#include "JSFZVTOP3.h"
#include "JSFGeoCFit.h"

#include "ANLPairCombiner.h"
#include "TNtupleD.h"
#include <sstream>
#include <vector>
#include <utility>
#include <map>

#include "TLorentzVector.h"
#include "TVector3.h"

void BoostJet(ANLPair w, ANLJet &j1, ANLJet &j2, Double_t &abscosj1, Double_t &abscosj2);

ClassImp(JetAnalysis)

static const Double_t kMassW   = 80.22;         // W mass               
static const Double_t kMassZ   = 91.18;         // Z mass           
static const Double_t kMassH   = 120.0;         // Higgs mass          
static const Double_t kSigmaMw =   8.0;         // W mass resolution    
static const Double_t kSigmaMz =   8.0;         // Z mass resolution 
static const Double_t kSigmaMh =   8.0;         // H mass resolution     
//const int maxjet_fit = 10;

TH1D *hEcone = 0;
TH2D *hLEcone = 0;

Int_t    JetAnalysis::fForcedNJets=0;
Double_t JetAnalysis::fYcut=0.005;
Int_t    JetAnalysis::fJetFinderAlgorithm=0;

//_____________________________________________________________________________
JetAnalysis::JetAnalysis(const char *name, const char *title)
       : JSFModule(name,title)
{

  fForcedNJets=gJSF->Env()->GetValue("JetAnalysis.ForcedNJets",2);
  sscanf(gJSF->Env()->GetValue("JetAnalysis.Ycut","0.001"),"%lg",&fYcut);
  //sscanf(gJSF->Env()->GetValue("JetAnalysis.Ycut","0.004"),"%lg",&fYcut);
  fJetFinderAlgorithm=gJSF->Env()->GetValue("JetAnalysis.JetFinderAlgorithm",0);
  
}

//_____________________________________________________________________________
JetAnalysis::~JetAnalysis()
{
  //  delete fJetEvent;
}

//_____________________________________________________________________________
Bool_t JetAnalysis::Process(Int_t nev)
{

  std::cout << "################# " << nev << std::endl;
  TObjArray *particles=new TObjArray(1000);
  ShapeAnalysis(particles);
  JetClustering(particles);
  particles->Delete();
  delete particles;
  return kTRUE;

}

//_____________________________________________________________________________
void JetAnalysis::JetClustering(TObjArray *tracks)
{
  if (!hEcone) hEcone = new TH1D("hEcone", "Econe", 50, 0., 100.);
  if (!hLEcone) hLEcone = new TH2D("hLEcone", "LEcone for lepton", 120, 0., 120., 120, 0., 120.);

  Double_t wmass=-1.;
  Double_t wene=-1.;
  Double_t wmom=-1.;

  Double_t nrene  = -1.;
  Double_t nrmass = -1.;
  Double_t cosnr  = -10.;

  Double_t nb=0.;
  Double_t yvalue=0.;
  Double_t tote=0.; 
  Double_t lepe=-10.;
  Double_t coslep=-10;

  Double_t xCosCone  = 0.94;
  Double_t xEconeCut = 5.;

  Double_t lepchg = -9999;
  Double_t lepgenid = 0.;
  Double_t lepmom = -10.;

  Double_t Pxl   =0.;
  Double_t Pyl   =0.;
  Double_t Pzl   =0.;

  Double_t PxW   =0.;
  Double_t PyW   =0.;
  Double_t PzW   =0.;

  Double_t cosWl = -10.;

  Double_t misse  = -10.;
  Double_t missm  = -10.;
  Double_t missp  = -10.;
  Double_t misspt = -10.;

  //----------------------------------------------------
  // jet clustering
  //----------------------------------------------------
  // remove lepton track from track-set used for jet-clustering
  TIter next(tracks);
  Double_t  maxe   = -10.;
  ANLTrack *leptrk = 0;
  ANLTrack *trk    = 0;
  while ((trk = static_cast<ANLTrack *>(next()))) {
    Double_t econe = trk->GetConeEnergy(xCosCone, tracks);
    //    std::cout << "lep ene: " << trk->IsLepton() << ", " << trk->E() << ", " << econe << std::endl;
    //    if(trk->IsLepton()){
    if(econe < xEconeCut){
      if(trk->E() > maxe){
	maxe   = trk->E();
        leptrk = trk;
      }
    
    }
  }
  if (!leptrk) return;

  lepe   = leptrk->E();
  lepmom = leptrk->GetMag();

  Pxl = leptrk->Px();
  Pyl = leptrk->Py();
  Pzl = leptrk->Pz();

  coslep = leptrk->CosTheta();
  lepchg = leptrk->GetCharge();

  /*  ANLTrack *tp;
  if (tp->IsLepton() ==  1){
  lepgenid = leptrk->GetLTKCLTrack()->GetCDC()->GetGenID();
  }
  */
  leptrk->Lock();
  if (lepchg == 0.) cerr << "lepchg = " << lepchg << endl;

  // perform jet-clustering
  ANLDurhamJetFinder jclust;
  jclust.Initialize(*tracks);
  jclust.ForceNJets(fForcedNJets);
  Int_t njets = jclust.GetNjets();

  if ( njets == fForcedNJets ) {
    //-----------------------------------------------------
    // Does jet pairing
    //-----------------------------------------------------
    TObjArray &jets = jclust.GetJets();
    yvalue = jclust.GetYcut();

    // --- Reconstruct W-boson
    ANLPairCombiner wcandidates(jets,jets);  
    ANLPair &w = *static_cast<ANLPair *>(wcandidates());
    wmass = w.GetMass();     
    wene  = w.E();     
    wmom  = w.Vect().Mag();

    PxW   = w.Px();
    PyW   = w.Py();
    PyW   = w.Pz();

    tote = lepe + wene;

    ANL4DVector qvnr = w + *leptrk;
    nrene  = qvnr.E();
    nrmass = qvnr.GetMass();
    cosnr  = qvnr.CosTheta();

    cosWl  = (PxW*Pxl+PyW*Pyl+PzW*Pzl) /( wmom * lepmom);

   
    ANL4DVector ecm = ANL4DVector(500., 0., 0., 0.);
    ANL4DVector recoil = ecm - qvnr;
    misse = recoil.E();
    missm = recoil.GetMass();
    missp = sqrt(recoil.Px()*recoil.Px()+recoil.Py()*recoil.Py()+recoil.Pz()*recoil.Pz());
    misspt = recoil.Pt();      
  
    // --- b-tagging ------------------------------------------------------
    ANLVTXTagger btag(3,3);
    nb = btag(*(ANLJet *)w[0]) + btag(*(ANLJet *)w[1]);

  }
  

  //-------------------
  // Prepare Ntuple
  //-------------------

  static TNtupleD *hEvt = 0;
  if (!hEvt) {
    stringstream tupstr;
    tupstr << "tote:wmass:wene:wmom:coslep:nb:lepe:lepchg:nrene:nrmass:cosnr:lepgenid:cosWl:yvalue:misse:missm:missp:misspt" << ends;
    hEvt = new TNtupleD("hEvt", "", tupstr.str().data());
  }
  
  //--                                                                      
  // Fill up Ntuple.                                              
  //--                                                                    
  Double_t data[100];
  data[ 0] = tote;
  data[ 1] = wmass;
  data[ 2] = wene;
  data[ 3] = wmom;
  data[ 4] = coslep;
  data[ 5] = nb;
  data[ 6] = lepe;
  data[ 7] = lepchg;
  data[ 8] = nrene;
  data[ 9] = nrmass;
  data[10] = cosnr;
  data[11] = lepgenid;
  data[12] = cosWl;
  data[13] = yvalue;
  data[14] = misse;
  data[15] = missm;
  data[16] = missp;
  data[17] = misspt;

  hEvt->Fill(data);

  leptrk->Unlock();
}

//_____________________________________________________________________________
void JetAnalysis::ShapeAnalysis(TObjArray *tracks)
{

  JSFGenerator *gen=(JSFGenerator*)gJSF->FindModule("JSFGenerator");
  JSFGeneratorBuf *geb=(JSFGeneratorBuf*)gen->EventBuf();

  Double_t ecm=geb->GetEcm();
  ecm = 500.;
  //std::cout << "##### Ecm = " << ecm << std::endl;


  //-----------------------------------------------------
  // Save tracks for Jet clustering
  //-----------------------------------------------------

  JSFSIMDST     *sds     = (JSFSIMDST*)gJSF->FindModule("JSFSIMDST");
  JSFSIMDSTBuf  *evt     = (JSFSIMDSTBuf*)sds->EventBuf();

  Int_t          ntrks   = evt->GetNLTKCLTracks();      // No. of tracks
  TObjArray     *trks    = evt->GetLTKCLTracks();       // combined tracks

  Int_t ntracks=0;
  Int_t nc=0;
  ANL4DVector qsum;
  for ( Int_t i = 0; i < ntrks; i++ ) {
    JSFLTKCLTrack *t = (JSFLTKCLTrack*)trks->UncheckedAt(i);
    ANLTrack *qt = new ANLTrack(t);
    tracks->Add(qt);           // track 4-momentum
    qsum += *qt;               // total 4-momentum
    ntracks++;
    if( t->GetCharge() != 0 ) { nc++ ; }
  }

  //
  ANLEventShape evshape;
  evshape.Initialize(*tracks);
}

//_____________________________________________________________________________
float JetAnalysis::mass(TLorentzVector tl[], int ijet[], const int njet) {
  // calculate the invariant mass of a collection of lorentzvectors
  TLorentzVector t(0,0,0,0);
  for (int i=0; i<njet; i++) {
    t+=tl[ijet[i]];
  }
  return t.M();
}

float JetAnalysis::constraint(TLorentzVector tl[], int ijet[], const int njet, float mforce) {
  // calculate value of constraint eqn (0 when constraint is fulfilled)
  return mass(tl, ijet, njet) - mforce;
}

void JetAnalysis::dCdf(TLorentzVector tl[], const int njet, int ijet[], float a[]) {
  TLorentzVector t0(0,0,0,0);
  for (int i=0; i<njet; i++) {
    t0+=tl[ijet[i]];
    float t1=0;
    float t2=0;
    float ei = tl[ijet[i]].E();
    float pi = tl[ijet[i]].Vect().Mag();
    for (int j=0; j<njet; j++) {
      if (j==i) continue;
      float ej = tl[ijet[j]].E();
      float pj = tl[ijet[j]].Vect().Mag();
      float costhij = TMath::Cos(tl[ijet[i]].Vect().Angle(tl[ijet[j]].Vect()));
      t1+=ej;
      t2+=pj*costhij;
    }
    a[i] = t1*pi/ei - t2;
  }
  float mass = t0.M();
  for (int i=0; i<njet; i++) {
    a[i]/=mass;
  }
  return;
}


Bool_t JetAnalysis::dofit(TLorentzVector tlOrig[], float perr[],                // original lorentz vectors, energy uncertainties
	   const int njet, const int ncons,                      // # jets, # constraints
	   int ijet[][maxjet_fit], int nforce[], float mforce[], // definition of constraints
	   float& chisq_out,                                     // fit chisq (output)
	   TLorentzVector* tlFit[]) {                            // fitter 4-vectors

  //
  // this routine actually does the fitting  
  //
  // it varies the mangitude of the lor.vec momentum, 
  //    leaving the mass and direction constant
  //

  const float smallMassDiff = 0.05; // if fitted-requested mass smaller than this, good enough

  cout << "hello from dofit! njet, ncons = " << njet << " " << ncons << endl;
  
  // check not too many jets
  if (njet>maxjet_fit) {
    cout << "anatreeAnalyse::dofit ERROR increase maxjet_fit" << endl;
    return false;
  }

  bool verbose = false; // print many details if true

  const int maxiter = 15; // maximum number of iterations in "fit"

  TVector3 temp;

  // define the parameters
  float paramsOrig[maxjet_fit];
  for (int i=0; i<njet; i++) paramsOrig[i] = tlOrig[i].Vect().Mag(); // magnitude of jet momentum

  TLorentzVector tl[njet];
  for (int i=0; i<njet; i++) tl[i]=tlOrig[i]; // lorentz vectors

  float params[njet];
  for (int i=0; i<njet; i++) params[i]=paramsOrig[i];

  float a[njet];

  if (verbose)
    cout << "params, errs = " << endl;

  // check input parameters make sense
  bool isok=true;
  for (int i=0; i<njet; i++) {
    if (verbose) cout << params[i] << " " << perr[i] << endl;
    if (params[i]<=0 || perr[i]<=0) { // negative |momentum| or uncertainty
      cout << "ERROR crazy parameters! " << endl;
      isok=false;
      for (int ii=0; ii<njet; ii++) cout << ii << " " << params[ii] << " " << perr[ii] << endl;
      break;
    }
  }
  if (!isok) return false;


  float chisq_old = -999;
  float chisq = -999;

  int iter=0;
  for (iter=0; iter<maxiter; iter++) { // loop over iterations

    if (verbose) {
      cout << "--------------------------------------------------" << endl;
      cout << "iteration " << iter << " parameters ";
      for (int i=0; i<njet; i++) cout << params[i] << " ";
      cout << endl;
    }

    TMatrix C0(ncons,1);
    TMatrix dF0(njet,1);
    for (int i=0; i<ncons; i++) {
      C0(i,0) = constraint(tl, &(ijet[i][0]), nforce[i], mforce[i]);
      if (verbose)
	cout << "constraint " << i << " = " << C0(i,0) <<
	  "  mass = " << mass(tl, &(ijet[i][0]), nforce[i]) << endl;
    }
    for (int i=0; i<njet; i++) {
      dF0(i,0) = params[i]-paramsOrig[i];
    }

    // define errors on jet energies
    TMatrix var(njet,njet);
    for (int i=0; i<njet; i++) {
      for (int j=0; j<njet; j++) {
	if (i==j) var(i,j) = perr[i]*perr[i];
	else var(i,j)=0;
      }
    }
    var.Invert();

    // write them to M
    TMatrix M(njet+ncons,njet+ncons);
    for (int i=0; i<njet; i++)
      for (int j=0; j<njet; j++)
	M(i,j) = var(i,j);

    // define D matrix = dC/df
    TMatrix D(ncons, njet);
    for (int i=0; i<ncons; i++)
      for (int j=0; j<njet; j++)
	D(i,j)=0;
    for (int i=0; i<ncons; i++) {
      dCdf(tl, nforce[i], &(ijet[i][0]), a);
      for (int j=0; j<nforce[i]; j++) {
	D(i,ijet[i][j]) = a[j];
      }
    }

    // add constraint part to M
    for (int i=0; i<ncons; i++) {
      for (int j=0; j<njet; j++) {
	M(njet+i, j) = D(i,j);
	M(j, njet+i) = D(i,j);
      }
    }

    if (verbose) {
      cout << "matrix M" << endl;
      for (int i=0; i<njet+ncons; i++) {
	for (int j=0; j<njet+ncons; j++)
	  cout << std::setprecision(5) << std::setw(10) << M(i,j) << " ";
	cout << endl;
      }
    }

    M.Invert(); // invert M

    // define RHS
    TMatrix Z(njet+ncons,1);
    for (int i=0; i<njet; i++)
      Z(i,0)=0;
    TMatrix R(D, TMatrix::kMult, dF0);
    R -= C0;
    for (int i=0; i<ncons; i++) {
      Z(njet+i,0)=R(i,0);
    }
    if (verbose) {
      cout << "matrix RHS" << endl;
      for (int i=0; i<njet+ncons; i++) cout << Z(i,0) << " "; cout << endl;
    }

    // get the corrections
    TMatrix Y(M, TMatrix::kMult, Z);
    if (verbose) {
      cout << "matrix Y" << endl;
      for (int i=0; i<njet+ncons; i++) cout << Y(i,0) << " "; cout << endl;
    }

    // update parameters (momentum magnitude)
    for (int i=0; i<njet; i++) {
      params[i]=paramsOrig[i]+Y(i,0);
    }

    // check that these are ok (ie positive)
    for (int i=0; i<njet; i++) {
      if (params[i]<0) { // |momentum| has gone negative; force it positive
	if (verbose) cout << iter << " correcting " << i << " " << params[i] << endl;
	params[i]=-params[i]/2;
      }
    }

    // update 4vectors
    for (int i=0; i<njet; i++) {
      temp.SetMagThetaPhi(params[i], tlOrig[i].Vect().Theta(),
			  tlOrig[i].Vect().Phi());
      tl[i].SetVectM(temp,tlOrig[i].M());
    }
    float newcons[ncons];
    for (int i=0; i<ncons; i++) {
      newcons[i] = constraint(tl, &(ijet[i][0]), nforce[i], mforce[i]);
    }


    if (verbose) {
      cout << "new constraints ";
      for (int i=0; i<ncons; i++)
	cout << newcons[i] << " ";
      cout << endl;
    }

    // calculate chisq
    TMatrix ediff(njet,1);
    for (int i=0; i<njet; i++) ediff(i,0) = params[i] - paramsOrig[i];
    TMatrix t2(var, TMatrix::kMult, ediff);
    TMatrix t3(ediff, TMatrix::kTransposeMult, t2);
    chisq = t3(0,0);
    if (verbose) cout << "chisq = " << chisq << endl;
    bool breakit = true;
    if (TMath::Abs((chisq-chisq_old)/chisq_old)>0.005) breakit = false; // chisq hasn't changed much - reached min
    for (int i=0; i<ncons; i++) {
      if (TMath::Abs(newcons[i])>smallMassDiff) breakit = false;
    }
    if (breakit) {
      if (verbose) cout << "breaking at iteration " << iter << endl;
      break;
    }
    chisq_old=chisq;
  } // loop over iterations

  if (verbose) {
    cout << "final params = " << endl;
    for (int i=0; i<njet; i++) cout << params[i] << " " ;
    cout << endl;
  }

  for (int i=0; i<njet; i++) *(tlFit[i]) = tl[i];


  chisq_out = chisq;

  return (iter<maxiter-1);
}

void BoostJet(ANLPair w, ANLJet &j1, ANLJet &j2, Double_t &abscosj1, Double_t &abscosj2)
{
  TVector3 ez = TVector3(0., 0. ,1.);

  TVector3 ewz = w.Vect().Unit();
  TVector3 ewx = ewz.Cross(ez).Unit();
  TVector3 ewy = ewz.Cross(ewx);

  TVector3 bstw = TVector3(0., 0., w.Vect().Mag()/w.E());
  //  cout << "hoge: " << w.Vect().Mag()/w.E() << endl; 

  ANL4DVector pj1 = ANL4DVector(j1.E(), j1.Vect()*ewx, j1.Vect()*ewy, j1.Vect()*ewz);
  pj1.Boost(-bstw);
  abscosj1 = TMath::Abs(pj1.CosTheta());

  ANL4DVector pj2 = ANL4DVector(j2.E(), j2.Vect()*ewx, j2.Vect()*ewy, j2.Vect()*ewz);
  pj2.Boost(-bstw);
  abscosj2 = TMath::Abs(pj2.CosTheta());

  //  Double_t p1x = pj1.Px();  
  //  Double_t p1y = pj1.Py();  
  //  Double_t p1z = pj1.Pz();  
  //  Double_t p2x = pj2.Px();  
  //  Double_t p2y = pj2.Py();  
  //  Double_t p2z = pj2.Pz();  
  
  //  acolj = (p1x*p2x+p1y*p2y+p1z*p2z)/(sqrt(p1x*p1x+p1y*p1y+p1z*p1z)*sqrt(p2x*p2x+p2y*p2y+p2z*p2z));

  return;
}
