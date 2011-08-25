//*LastUpdate : v.01.01 undefined by undefined
//*-- Author  : undefined undefined

///////////////////////////////////////////////////////////////////
//
//  StdhepMon
//
//  Analysis of StdHep data
//
//$Id: 
//  
//////////////////////////////////////////////////////////////////

#include <iostream>

#include "JSFSteer.h"
#include "StdhepMon.h"
#include "JSFReadStdHep.h"
#include "ANL4DVector.h"
#include "ANLEventShape.h"
#include "ANLJetFinder.h"

ClassImp(StdhepMon)

//_____________________________________________________________________________
  StdhepMon::StdhepMon(const char *name, const char *title)
       : JSFModule(name,title)
{
  //  Example to get parameter  
  fFlag = gJSF->Env()->GetValue("StdhepMon.Flag",3);
  sscanf(gJSF->Env()->GetValue("StdhepMon.Parameter","3.14"),"%g",&fParameter);
  fLastEventNumber=0;

}

//_____________________________________________________________________________
StdhepMon::~StdhepMon()
{
}

//_____________________________________________________________________________
Bool_t StdhepMon::Initialize()
{
  // Initialize
  fNprt=new TNtuple("ntp","Particle based ntuple",
		    "id:e:pt:cs:vr:vz");
  fNevt=new TNtuple("nte","Event based ntuple",
		    "etot:evis:pt:pz:nchg:ngam:nhad:echg:egam:ehad:thrust:csth:njets");
//		    "etot:evis:pt:pz:nchg,ngam,nhad,echg,egam,ehad,thrust:csth:njets:yc2j:csj2:ej2:yc4j:csj4:ej4");
  fN2j=new TNtuple("nj2","2 jet variables",
	           "yc2j:j0e:j0cs:j0m:j1e:j1cs:j1m");
  fN4j=new TNtuple("nj4","4 jet variables",
	           "yc4j:j0e:j0cs:j1e:j1cs:j2e:j2cs:j3e:j3cs:m01:m23:m02:m13:m03:m12"); 

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StdhepMon::Process(Int_t nev)
{
// 
  Int_t ievt=gJSF->GetEventNumber();
  if ( ievt < 50 || ( ievt<1000&&ievt%100==1) || ( ievt<10000&&ievt%1000==1) || ievt%5000==1) {
    std::cout << "Event No. " << ievt << std::endl;
  }
  fLastEventNumber=ievt;
  if( ievt < 10) {
    JSFReadStdHep *stdhep=(JSFReadStdHep*)gJSF->FindModule("JSFReadStdHep");
    JSFReadStdHepBuf *buf=(JSFReadStdHepBuf*)stdhep->EventBuf();
    buf->GetHEPEV4().Print();
  }


  JSFGenerator *gen=(JSFGenerator*)gJSF->FindModule("JSFGenerator");
  JSFGeneratorBuf *gbuf=(JSFGeneratorBuf*)gen->EventBuf();
  JSFGeneratorParticle *p=0;
  TClonesArray *gp=gbuf->GetParticles();
  if ( ievt < 5 ) {
    std::cerr << " GeneratorParticles=" << gbuf->GetNParticles() << std::endl;
    for( Int_t i=0;i<gbuf->GetNParticles();i++) {
      p=(JSFGeneratorParticle*)gp->UncheckedAt(i);
      if( i==0 ) { p->ls("form3,title"); }
      else {
        p->ls("form3");
      }
    }
  }

// ================================
  ANL4DVector psum(0.0,0.0,0.0,0.0);
  ANL4DVector psum0(0.0,0.0,0.0,0.0);
  Double_t echg=0.0;
  Double_t egam=0.0;
  Double_t ehad=0.0;
  Int_t nchg=0;
  Int_t ngam=0;
  Int_t nhad=0;

  TObjArray particles(1000);
  Bool_t dprint=false;

  for(Int_t i=0;i<gbuf->GetNParticles();i++) {
    p=(JSFGeneratorParticle*)gp->UncheckedAt(i);
    if( p->GetNDaughter() == 0 && p->GetStatus() == 1 ) {
        ANL4DVector pa(p->GetE(),p->GetPx(),p->GetPy(), p->GetPz());
        psum0+=pa;
        Int_t aid=TMath::Abs(p->GetID());
        if(   aid != 22 && aid != 11 && aid != 12 && aid != 13 && aid !=14 && aid != 16 
  	   && aid != 211 && aid != 130 && aid != 310 && aid != 321 && aid != 2212 && aid != 2112 ) {
             std::cerr << "=== Found illegal PID .. PID=" << p->GetID() << std::endl;
          dprint=true;
        }

        if( TMath::Abs(p->GetCharge())> 0.5 ) {
          nchg++;
          echg+=p->GetE();
        }
        else if( aid == 22 ) {
          ngam++;
          egam+=p->GetE();
        }
        else if ( aid !=12 && aid != 14 && aid != 16 ) {
          nhad++;
          ehad+=p->GetE();
        }
        Float_t vr=TMath::Sqrt(p->GetX()*p->GetX() + p->GetY()*p->GetY());
        fNprt->Fill(p->GetID(),pa.E(),pa.Pt(),pa.CosTheta(),vr,p->GetZ());
        if( aid ==12 || aid==14 || aid==16 ) continue;


        psum+=pa;
        particles.Add(new ANL4DVector(pa));
        
     }  // for stable particles
    
  }
  ANLEventShape evshape;
  evshape.Initialize(particles);
  Float_t thrust=evshape.GetThrust();
  Float_t csthrust=(evshape.thrustAxis())->CosTheta();
  Int_t njets=0;
  Float_t yc2j=0.0;
  Float_t yc4j=0.0;

  Double_t ycut=0.001;
  
  if( particles.GetEntries() > 4 ) {

    ANLJadeJetFinder jclust(ycut);
    jclust.Initialize(particles);
    jclust.FindJets();
    njets=jclust.GetNjets();
    if( njets >= 2 ) { 
      jclust.ForceNJets(2);
      yc2j=jclust.GetYmax();
      TObjArray &jets=jclust.GetJets();
      ANLJet *jet0=(ANLJet*)jets.UncheckedAt(0);
      ANLJet *jet1=(ANLJet*)jets.UncheckedAt(1);
      fN2j->Fill(yc2j,jet0->E(),jet0->CosTheta(),jet0->M(), jet1->E(),jet1->CosTheta(),jet1->M());

      ANLJadeJetFinder jclust4(ycut);
      jclust4.Initialize(particles);
      jclust4.FindJets();
      if( jclust4.GetNjets() >= 4 ) {
        jclust4.ForceNJets(4);
        yc4j=jclust4.GetYmax();
        TObjArray &jets4=jclust4.GetJets();
        ANLJet *j0=(ANLJet*)jets4.UncheckedAt(0);
        ANLJet *j1=(ANLJet*)jets4.UncheckedAt(1);
        ANLJet *j2=(ANLJet*)jets4.UncheckedAt(2);
        ANLJet *j3=(ANLJet*)jets4.UncheckedAt(3);
	Double_t m01=(*j0+*j1).M();
	Double_t m02=(*j0+*j2).M();
	Double_t m03=(*j0+*j3).M();
	Double_t m12=(*j1+*j2).M();
	Double_t m13=(*j1+*j3).M();
	Double_t m23=(*j2+*j3).M();

        fN4j->Fill(yc4j,j0->E(),j0->CosTheta(),j1->E(),j1->CosTheta(),j2->E(),j2->CosTheta(),
	j3->E(),j3->CosTheta(),m01,m23,m02,m13,m03,m12);
      }	

    }
  }


  fNevt->Fill(Float_t(psum0.E()),Float_t(psum.E()),Float_t(psum.Pt()),Float_t(psum.Pz()),float(nchg),float(ngam),
        float(nhad),Float_t(echg), Float_t(egam),Float_t(ehad),thrust,csthrust,float(njets));
//  fNevtp->Fill(thrust, csthrust, yc2j,cs2j,e2j,yc4j,cs4j,e4j);


  particles.Delete();

// === Forced printout of particle list of this event.
  static int noutevt=0;
  if( dprint ) { 
      std::cerr << "=== Forced print of particle contents ... event number is " << ievt << std::endl;
      std::cerr << "=== Nout count is " << noutevt << " ( no print if noutevt > 100 ) " << std::endl;
      if( noutevt < 100 ) {
        noutevt++;
        std::cerr << " GeneratorParticles=" << gbuf->GetNParticles() << std::endl;
        for( Int_t i=0;i<gbuf->GetNParticles();i++) {
          p=(JSFGeneratorParticle*)gp->UncheckedAt(i);
          p->ls();
        }
	std::cerr << "=== End of particle list output === " << std::endl;
     } 
  }

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t StdhepMon::EndRun()
{
  std::cout << "Last event number was " << fLastEventNumber << std::endl;
  return kTRUE;
}



