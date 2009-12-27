#ifndef __JetAnalysis__
#define __JetAnalysis__

///////////////////////////////////////////////////////////////////////
//
//  JetAnalysis
//
//  Sample program for Jet Analysis
//
//$Id:
//
///////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "ANL4DVector.h"
#include "JSFModule.h"

const int maxjet_fit = 10;

// *******************************************************
class JetAnalysis : public JSFModule 
{
 protected:

  static Int_t    fForcedNJets;  // Event is clustered to this number of jets
  static Double_t fYcut;         // Initial Y Cut value
  static Int_t    fJetFinderAlgorithm;  // =0(JadeFinder), 1=(JadeEJetFinder),2=(DurhamJetFinder)

 public:
  JetAnalysis(const char *name="JetAnalysis", 
		  const char *title="JetAnalysis");
  virtual ~JetAnalysis();

  static Int_t    GetForcedNJets(){ return fForcedNJets;}
  static Double_t GetYcut(){ return fYcut; }
  static Int_t    GetJetFinderAlgorithm(){ return fJetFinderAlgorithm; }
  void     SetEtrackCut(Double_t x) { fEtrackCut = x; }

  //  virtual Bool_t Initialize();
  //  virtual Bool_t BeginRun(Int_t runno=1);
  //  virtual Bool_t EndRun();
  //  virtual Bool_t Terminate();

  virtual Bool_t Process(Int_t ev=1);
  void JetClustering(TObjArray *particles);
  void ShapeAnalysis(TObjArray *particles);
  float mass(TLorentzVector tl[], int ijet[], const int njet);
  float constraint(TLorentzVector tl[], int ijet[], const int njet, float mforce);
  void dCdf(TLorentzVector tl[], const int njet, int ijet[], float a[]);

  Bool_t dofit(TLorentzVector tlOrig[], float perr[],                // original lorentz vectors, energy uncertainties
	   const int njet, const int ncons,                      // # jets, # constraints
	   int ijet[][maxjet_fit], int nforce[], float mforce[], // definition of constraints
	   float& chisq_out,                                     // fit chisq (output)
	   TLorentzVector* tlFit[]);                           // fitter 4-vectors

 private:
  Double_t fEtrackCut;    // cut on track energy
  
  ClassDef(JetAnalysis, 1)  // JetAnalysis class


};

#endif










