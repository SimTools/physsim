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

// *******************************************************
class JetAnalysis : public JSFModule 
{
 protected:

  static Int_t    fForcedNJets;  // Event is clustered to this number of jets
  static Double_t fYcut;         // Initial Y Cut value
  static Int_t    fJetFinderAlgorithm;  // =0(JadeFinder), 1=(JadeEJetFinder),2=(DurhamJetFinder)

  static int evn; //prepare for check clustering @ "JetAnalysis.cxx"

 public:
  JetAnalysis(const char *name="JetAnalysis", 
		  const char *title="JetAnalysis");
  virtual ~JetAnalysis();

  static Int_t    GetForcedNJets(){ return fForcedNJets;}
  static Double_t GetYcut(){ return fYcut; }
  static Int_t    GetJetFinderAlgorithm(){ return fJetFinderAlgorithm; }

  //  virtual Bool_t Initialize();
  //  virtual Bool_t BeginRun(Int_t runno=1);
  //  virtual Bool_t EndRun();
  //  virtual Bool_t Terminate();

  virtual Bool_t Process(Int_t ev=1);
  void JetClustering(TObjArray *particles);
  void ShapeAnalysis(TObjArray *particles);

  ClassDef(JetAnalysis, 1)  // JetAnalysis class

};

#endif










