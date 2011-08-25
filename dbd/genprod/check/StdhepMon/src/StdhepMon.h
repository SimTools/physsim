#ifndef __StdhepMon__
#define __StdhepMon__

///////////////////////////////////////////////////////////////////////
//
//  StdhepMon
//
//  Analysis of StdHep data
//
//$Id:
//
///////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include "JSFModule.h"
#include "TNtuple.h"


// *******************************************************
class StdhepMon : public JSFModule 
{
 protected:
  Int_t   fFlag;      // An example of Int_t parameter
  Float_t fParameter; // An example of Float_t parameter
  Int_t   fLastEventNumber;

  TNtuple *fNprt;  // Particle based ntuple
  TNtuple *fNevt;  // Event based ntuple
  TNtuple *fN2j;   // Ntuple for 2 jet
  TNtuple *fN4j;   // Ntuple for 4 jet


 public:
  StdhepMon(const char *name="StdhepMon", 
		  const char *title="StdhepMon");
  virtual ~StdhepMon();

  virtual Bool_t Initialize();
  virtual Bool_t EndRun();
  virtual Bool_t Process(Int_t ev=1);
  
  ClassDef(StdhepMon, 1)  // StdhepMon class

};

#endif










