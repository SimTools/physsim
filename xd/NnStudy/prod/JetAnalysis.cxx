///////////////////////////////////////////////////////////////////
//
//  JetAnalysis
//
//  Sample program for Jet Analysis
//
//$Id: 
//  
//////////////////////////////////////////////////////////////////
#include <iostream>
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
#include "TMath.h"
#include "ANLPairCombiner.h"
#include "TNtupleD.h"
#include "TVector3.h"
#include <sstream>
#include <vector>
#include <utility>
#include <map>

#include <fstream>

Bool_t SolveKinematics(Double_t sqrt_s, Double_t mW, const ANL4DVector &p1, const ANL4DVector &p2, Double_t abscosWH1[2], Double_t abscosWH2[2]);
void BoostJet(ANLPair w, ANLJet &j1, ANLJet &j2, Double_t &cosj1, Double_t &cosj2);


ClassImp(JetAnalysis)

static const Double_t kMassW   = 79.9661182;         // W mass               
static const Double_t kMassZ   = 91.18;         // Z mass           
static const Double_t kMassH   = 134.0;         // Higgs mass          
//static const Double_t kSigmaMw =   8.0;         // W mass resolution    
static const Double_t kSigmaMw =   4.0;         // W mass resolution    
static const Double_t kSigmaMz =   8.0;         // Z mass resolution 
static const Double_t kSigmaMh =   8.0;         // H mass resolution     

int JetAnalysis::evn = 1; // for check Clustering

Int_t    JetAnalysis::fForcedNJets=0;
Double_t JetAnalysis::fYcut=0.005;
Int_t    JetAnalysis::fJetFinderAlgorithm=0;

//_____________________________________________________________________________
JetAnalysis::JetAnalysis(const char *name, const char *title)
       : JSFModule(name,title)
{

  fForcedNJets=gJSF->Env()->GetValue("JetAnalysis.ForcedNJets",4);
  sscanf(gJSF->Env()->GetValue("JetAnalysis.Ycut","0.001"),"%lg",&fYcut);
  //  sscanf(gJSF->Env()->GetValue("JetAnalysis.Ycut","0.004"),"%lg",&fYcut);
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
  //Double_t wn1mass=-1.;
  //Double_t wn2mass=-1.;
  //Double_t chi2wnwn=100.;
  Double_t w1mass=-1.;
  Double_t w2mass=-1.;
  Double_t chi2ww=100.;
  Int_t nb=0;

  Double_t ene1=0.;
  Double_t ene2=0.;

  Double_t recm=-1.;
  Double_t misspt=-1.;

  Double_t p1x=0.;
  Double_t p1y=0.;
  Double_t p1z=0.;
  Double_t p1=0.;
  Double_t p2x=0.;
  Double_t p2y=0.;
  Double_t p2z=0.;
  Double_t p2=0.;

  Double_t beta1=0.;
  Double_t beta2=0.;

  Double_t costheta1=0.;
  Double_t costheta2=0.;
  Double_t acop=0.;

  //Double_t cosWH1=0.;
  //Double_t cosWH2=0.;
  Double_t abscosWH1a=0.;
  Double_t abscosWH1b=0.;
  Double_t abscosWH2a=0.;
  Double_t abscosWH2b=0.;

  Double_t abscosj1w1=0.;
  Double_t abscosj2w1=0.;
  Double_t abscosj3w2=0.;
  Double_t abscosj4w2=0.;

  Double_t ene_jet[4]; // for getting Jet Energy before clustering
  Double_t ene_j1w1 = 0.; // for getting Jet Energy after clustering
  Double_t ene_j2w1 = 0.;
  Double_t ene_j3w2 = 0.;
  Double_t ene_j4w2 = 0.;


  //----------------------------------------------------
  // jet clustering
  //----------------------------------------------------
  //  if( tracks->GetEntries() <= fForcedNJets ) return;  

  /*
  ANLJetFinder *jclust=0;
  if( fJetFinderAlgorithm == 0 ) {
    jclust=new ANLJadeJetFinder(fYcut);
  }
  else if( fJetFinderAlgorithm == 1 ) {
  }
  else {
  }

  jclust->Initialize(*tracks);
  jclust->FindJets();
  Int_t njets=jclust->GetNjets();
  */

  ANLDurhamJetFinder* jclust = new ANLDurhamJetFinder();
  jclust->Initialize(*tracks);
  
  if(jclust->GetNjets() < 4)
    {
      std::cerr << "less than 4 jests" << std::endl;
      return;
    }

  Int_t njets = 4;
  jclust->ForceNJets(njets);

  ANLVTXTagger  btag(5, 3);

  if( njets == fForcedNJets ) {  // return after delete jclust object.
    //  if( njets >= fForcedNJets ) {  // return after delete jclust object.
    //-----------------------------------------------------
    // Does jet pairing
    //-----------------------------------------------------
    std::vector<ANLPair*> solvec;
    
    TObjArray solutions(10);
    ANLPair *w1p, *w2p;                                     
    TObjArray &jets = jclust->GetJets();
    ANLPairCombiner w1candidates(jets,jets);  

    //***********************************
    //* Get Jet Energy before clustering
    //***********************************
    ANLJet *jetp;
    TIter nextjet(&jets);
    int i_jet = 0;
    while ( (jetp = (ANLJet*)nextjet()) )
    {
      ANLJet &jet = *jetp;
      ene_jet[i_jet] = jet()(0);
      i_jet++;
    }


    // --- Search W1 pair                                               
    while ((w1p = (ANLPair *)w1candidates())) {                        
      ANLPair &w1 = *w1p;                                         
      w1.LockChildren();                                             
      Double_t w1mass = w1().GetMass();                              
      ANLPairCombiner w2candidates(w1candidates);                 
      w2candidates.Reset();
      // --- Search W2 pair      
      while ((w2p = (ANLPair *)w2candidates())) {              
	ANLPair &w2 = *w2p;             
	if (w2.IsLocked()) continue;                            
	Double_t w2mass = w2().GetMass();                         
	Double_t chi2w = TMath::Power((w1mass - kMassW)/kSigmaMw, 2.)      
	  + TMath::Power((w2mass - kMassW)/kSigmaMw, 2.);         
	solutions.Add(new ANLPair(w1p,w2p,chi2w));
	if(chi2w<chi2ww){
	  chi2ww=chi2w;
	}
      }                  
      w1.UnlockChildren();                                      
    }    

    if (solutions.GetEntries() > 1){
      solutions.Sort();
      ANLPair &sol = *static_cast<ANLPair *>(solutions.At(0));
      ANLPair &w1a   = *static_cast<ANLPair *>(sol[0]);
      ANLPair &w2a   = *static_cast<ANLPair *>(sol[1]);

      /*
      if(btag(*(ANLJet *)w1a[0])) nb++;
      if(btag(*(ANLJet *)w1a[1])) nb++;
      if(btag(*(ANLJet *)w2a[0])) nb++;
      if(btag(*(ANLJet *)w2a[1])) nb++;
      */

      w1mass = w1a().GetMass();
      w2mass = w2a().GetMass();

      ene1 = w1a().E();
      ene2 = w2a().E();

      ANL4DVector ecm = ANL4DVector(1000., 0., 0., 0.);
      ANL4DVector recoil = ecm - w1a - w2a;
      recm = recoil.GetMass();
      misspt = recoil.GetPt();

      p1x = w1a().Px();
      p1y = w1a().Py();
      p1z = w1a().Pz();
      p1 = w1a().GetMag();
      p2x = w2a().Px();
      p2y = w2a().Py();
      p2z = w2a().Pz();
      p2 = w2a().GetMag();

      beta1 = p1/ene1;
      beta2 = p2/ene2;

      costheta1 = w1a().CosTheta();
      costheta2 = w2a().CosTheta();

      acop = w1a().Acop(w2a);

      //======================
      // Get production angle
      //======================
      Double_t abscosWH1[2], abscosWH2[2];
      if(SolveKinematics(ecm.E(), kMassW, w1a, w2a, abscosWH1, abscosWH2)){
	abscosWH1a = abscosWH1[0];
	abscosWH2a = abscosWH2[0];
	abscosWH1b = abscosWH1[1];
	abscosWH2b = abscosWH2[1];
      }
      else{
	abscosWH1a = 12345;
	abscosWH2a = 12345;
	abscosWH1b = 12345;
	abscosWH2b = 12345;	
      }

      //===================
      // Set jet 4D vector
      //===================
      ANLJet &j1w1 = *static_cast<ANLJet *>(w1a[0]);
      ANLJet &j2w1 = *static_cast<ANLJet *>(w1a[1]);
      ANLJet &j3w2 = *static_cast<ANLJet *>(w2a[0]);
      ANLJet &j4w2 = *static_cast<ANLJet *>(w2a[1]);

      //===============================
      // Get jet angle at W rest frame
      //===============================
      BoostJet(w1a, j1w1, j2w1, abscosj1w1, abscosj2w1);
      BoostJet(w2a, j3w2, j4w2, abscosj3w2, abscosj4w2);

      //**************************************
      //* Get Jet information afte clustering
      //**************************************
      ene_j1w1 = j1w1()(0);
      ene_j2w1 = j2w1()(0);
      ene_j3w2 = j3w2()(0);
      ene_j4w2 = j4w2()(0);
    }
  
  }


//   //===================
//   //= Check clustering
//   //===================
//   char *fname = "test_out.dat";
//   std::ofstream fout;

//   //int evn = 1;

//   //================
//   //= Initilization
//   //================
//   if(evn==1){
//     fout.open(fname);
//     fout.close();
//   }

//   //==========
//   //= Writing
//   //==========
//   fout.open(fname, ios::app);
//   /*
//   fout << evn << "\t" 
//        << ene_jet[0] << "\t" << ene_jet[1] << "\t" << ene_jet[2] << "\t" << ene_jet[3] << "\t" 
//        << ene_j1w1 << "\t" << ene_j2w1 << "\t" << ene_j3w2 << "\t" << ene_j4w2 << "\t" 
//        << ene1 << "\t" << ene2 << "\t" 
//        << costheta1 << "\t" << costheta2
//        << std::endl;
//   */
//   fout << abscosWH1a <<  "\t" <<  abscosWH1b << "\t" <<  abscosWH2a << "\t" << abscosWH2b << std::endl;

//   fout.close();
//  evn = evn + 1;

  //----------------
  // Prepare Ntuple
  //----------------
  static TNtupleD *hEvt = 0;
  if (!hEvt) {
    stringstream tupstr;
    //    tupstr << "chi2h:chi2z"                                             << ":"
      //           << "h1mass:h2mass:z1mass:z2mass:njet:cale:emcale:hdcale"     << ""
    tupstr << "chi2w:w1mass:w2mass"
	   << ":"
	   << "ene1:ene2:recm:p1x:p1y:p1z:p1:p2x:p2y:p2z:p2:beta1:beta2:costheta1:costheta2:njets:nb:misspt:acop:abscosWH1a:abscosWH1b:abscosWH2a:abscosWH2b:abscosj1w1:abscosj2w1:abscosj3w2:abscosj4w2:ene_jet[0]:ene_jet[1]:ene_jet[2]:ene_jet[3]:ene_j1w1:ene_j2w1:ene_j3w2:ene_j4w2"
           << ends;
    hEvt = new TNtupleD("hEvt", "", tupstr.str().data());
  }

  //--                                                                      
  // Fill up Ntuple.                                              
  //--                                                                    
  std::multimap<Double_t, std::vector<ANLPair*> >::iterator itr;
  Double_t data[100];

  data[0] = chi2ww;
  data[1] = w1mass;
  data[2] = w2mass;

  data[3] = ene1;
  data[4] = ene2;
  data[5] = recm;

  data[6] = p1x;
  data[7] = p1y;
  data[8] = p1z;
  data[9] = p1;
  data[10] = p2x;
  data[11] = p2y;
  data[12] = p2z;
  data[13] = p2;

  data[14] = beta1;
  data[15] = beta2;

  data[16] = costheta1;
  data[17] = costheta2;

  data[18] = njets;
  data[19] = nb;

  data[20] = misspt;
  data[21] = acop;

  data[22] = abscosWH1a;
  data[23] = abscosWH1b;
  data[24] = abscosWH2a;
  data[25] = abscosWH2b;

  data[26] = abscosj1w1;
  data[27] = abscosj2w1;
  data[28] = abscosj3w2;
  data[29] = abscosj4w2;

  data[30] = ene_jet[0];
  data[31] = ene_jet[1];
  data[32] = ene_jet[2];
  data[33] = ene_jet[3];

  data[34] = ene_j1w1;
  data[35] = ene_j2w1;
  data[36] = ene_j3w2;
  data[37] = ene_j4w2;

  hEvt->Fill(data);

  delete jclust; 

}

Bool_t SolveKinematics(Double_t sqrt_s, Double_t mW, const ANL4DVector &p1, const ANL4DVector &p2, Double_t abscosWH1[2], Double_t abscosWH2[2])
{
  //=========================
  // Prepare base 3D vector
  //=========================
  ANL3DVector e1(p1.Get3D().Unit());
  ANL3DVector e2(p2.Get3D().Unit());

  ANL3DVector ex,ey,ez;
  Double_t cos12  = e1*e2;

  if (1-cos12*cos12 <= 0.) return kFALSE;

  ex = e1;
  ey = 1/TMath::Sqrt(1-cos12*cos12) * (e2 - cos12*e1);
  ez = (ex.Cross(ey)).Unit();

  //===========================
  // Project WH to base vector 
  //===========================
  Double_t mAH = 81.8515;
  Double_t mWH = 368.154;

  Double_t eneWH, pWH, betaWH;
  eneWH  = sqrt_s/2;
  pWH    = TMath::Sqrt(eneWH*eneWH - mWH*mWH);
  betaWH = pWH/eneWH;

  Double_t eWH1e1, eWH2e2;
  eWH1e1 = (sqrt_s*p1.E() - (mWH*mWH - mAH*mAH + mW*mW))/(2 * pWH * p1.GetMag());
  eWH2e2 = (sqrt_s*p2.E() - (mWH*mWH - mAH*mAH + mW*mW))/(2 * pWH * p2.GetMag());

  Double_t cx, cy, cz[2];
  cx = eWH1e1;
  cy = 1/TMath::Sqrt(1-cos12*cos12) * (-eWH2e2 - cos12*eWH1e1);

  if(1 - cx*cx - cy*cy < 0.) return kFALSE;

  cz[0] =  TMath::Sqrt(1 - cx*cx - cy*cy);
  cz[1] = -TMath::Sqrt(1 - cx*cx - cy*cy);

  ANL3DVector eWH1[2]; 
  for(Int_t i=0; i<2; i++){
    eWH1[i] = cx*ex + cy*ey + cz[i]*ez;

    abscosWH1[i] = TMath::Abs(eWH1[i](3));
    abscosWH2[i] = abscosWH1[i];
  }

  return kTRUE;
}


void BoostJet(ANLPair w, ANLJet &j1, ANLJet &j2, Double_t &abscosj1, Double_t &abscosj2)
{
  TVector3 ez = TVector3(0., 0. ,1.);

  TVector3 ewz = w.Vect().Unit();
  TVector3 ewx = ewz.Cross(ez).Unit();
  TVector3 ewy = ewz.Cross(ewx);

  TVector3 bstw = TVector3(0., 0., w.Vect().Mag()/w.E());

  ANL4DVector pj1 = ANL4DVector(j1.E(), j1.Vect()*ewx, j1.Vect()*ewy, j1.Vect()*ewz);
  pj1.Boost(-bstw);
  abscosj1 = TMath::Abs(pj1.CosTheta());

  ANL4DVector pj2 = ANL4DVector(j2.E(), j2.Vect()*ewx, j2.Vect()*ewy, j2.Vect()*ewz);
  pj2.Boost(-bstw);
  abscosj2 = TMath::Abs(pj2.CosTheta());

  return;
}


//_____________________________________________________________________________
void JetAnalysis::ShapeAnalysis(TObjArray *tracks)
{

  JSFGenerator *gen=(JSFGenerator*)gJSF->FindModule("JSFGenerator");
  JSFGeneratorBuf *geb=(JSFGeneratorBuf*)gen->EventBuf();

  Double_t ecm=geb->GetEcm();
  ecm = 1000.;
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

