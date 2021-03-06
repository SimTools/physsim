//****************************************************
//*
//*  Dummy UserAnalysis Script to be used with gui.C 
//*  Functions in this scripts are called by GUIMainMacro.C  
//*  
//*  In this script, global functions, UserInitialize() nad UserAnalysis()
//*  Must be defined.  UserInitialize() is called at job initialization 
//*  to initialize histograms, etc.  It is also called when reset historgam
//*  menu is selected from JSF Control Panel.  UserAnalysis() is called 
//*  during the event loop after executing Process functions of JSF modules
//*  and display event daia but before display histogram.  See GetEvent()
//*  function in GUIMainMacro.C, where UserAnalysis() is called.
//*
//*  UserSetOptions() and UserTerminate() may be defined in this file.
//*  UserSetOptions() are called before declaration of JSF modules.
//*  It can be used to set parameters optional for user analsis.
//*  UserTerminate() is called during the JSF termination process.
//*
//*  When runmode=4, UserModuleDefine() must be defined in this file.
//*  It is used to define JSF modules specific to user analysis.  
//*  
//*$Id$
//*
//****************************************************

// JSFWriteStdHep *stdhep;

TNtuple *nT;


//______
void UserAnotherModules()
{
  cerr << "UserModuleDefine was called." <<endl;
  stdhep=new JSFWriteStdHep();

}

//_________________________________________________________
void UserInitialize()
{
  nT=new TNtuple("nT","ntuple","esum:eh");
}

//_________________________________________________________
void UserAnalysis()
{
  Int_t ievt=gJSF->GetEventNumber();
  if( ievt < 100 || ( ievt < 1000 && ievt%100 == 1 ) 
	|| ( ievt < 10000 && ievt%1000 == 1 ) || ievt%10000 ==1 ) {
     std::cout << "Event No. " << gJSF->GetEventNumber() << std::endl;
   }
   Bool_t DODUMP=false;
//   if ( ievt < 2 ) { DODUMP=true; }
//   if ( ievt < 6 || ievt == 44  ) { DODUMP=true; }
//   if ( ievt < 10  ) { DODUMP=true; }
   Bool_t DOMON=true; 

//  if( gJSF->GetEventNumber() < 5 ) {
    Double_t sesum=0.0;
    Double_t spxsum=0.0;
    Double_t spysum=0.0;
    Double_t spzsum=0.0;
    Int_t spncall=0;
    JSFSpring *spr=(JSFSpring*)gJSF->FindModule("JSFSpring");
    JSFSpringBuf *sbuf=(JSFSpringBuf*)spr->EventBuf();
    TClonesArray *spa=sbuf->GetPartons();
    JSFSpringParton *sp=0;
    Float_t ehiggs=0;
    if ( DODUMP ) {
      std::cout << "Printout SpringPartons in UserAnalysis ievt=" << ievt << std::endl;
    }
    for(Int_t i=0;i<sbuf->GetNpartons() ; i++) {
	sp=(JSFSpringParton*)spa->UncheckedAt(i);
	if( DODUMP ) {
          if( i==0 ) { sp->ls("form2,title"); }
          else { sp->ls("form2"); }
        }
        if( sp->GetID() == 25 ) {
          ehiggs=sp->GetE();
        }

        if( sp->GetNDaughter() == 0 ) {
          sesum+=sp->GetE();
	  spxsum+=sp->GetPx();
	  spysum+=sp->GetPy();
	  spzsum+=sp->GetPz();
          spncall++;
        }
    }
    if ( DODUMP ) {
    std::cout << " Spring ncall=" << spncall 
	<< "Esum=" << sesum << " Pxsum=" << spxsum
	<< " Pysum=" << spysum << " Pzsum=" << spzsum << " ehiggs=" << ehiggs
	<< std::endl;
    }
    if ( sesum > 1005.0 ) {
      std::cout << "Spring particles esum > 1005.0  ievt=" << ievt ;
      std::cout 
	<< "Esum=" << sesum << " Pxsum=" << spxsum
	<< " Pysum=" << spysum << " Pzsum=" << spzsum << " ehiggs=" << ehiggs
	<< std::endl;
    } 

    nT->Fill(sesum,ehiggs);

  if ( DOMON ) {
    JSFGenerator *gen=(JSFGenerator*)gJSF->FindModule("JSFGenerator");
    JSFGeneratorBuf *gbuf=(JSFGeneratorBuf*)gen->EventBuf();
    Double_t esum=0.0;
    Double_t pxsum=0.0;
    Double_t pysum=0.0;
    Double_t pzsum=0.0;
    Int_t nsum=0;
    TClonesArray *gp=gbuf->GetParticles();
    JSFGeneratorParticle *p=0;
//    if( DODUMP ) {
//    std::cout << "Printout GeneratorParticles in UserAnalysis ievt=" << ievt 
//	<< " # of generator particles=" << gbuf->GetNParticles() << std::endl;
//    }
    for( Int_t i=0;i<gbuf->GetNParticles();i++) {
      p=(JSFGeneratorParticle*)gp->UncheckedAt(i);
      if( DODUMP ) {
        if( i==0 ) { p->ls("form3,title"); }
        else { p->ls("form3"); }
      }
//      if( p->GetNDaughter() == 0 ) {
        if( p->GetStatus()==1 ) {
	nsum++;
        esum+=p->GetE();
        pxsum+=p->GetPx();
        pysum+=p->GetPy();
        pzsum+=p->GetPz();      
//        p->ls();
      }

   }
   if( esum > 1010.0 || DODUMP ) {
       if ( esum > 1010.0 ) {
         std::cout << " Esum of generator particles > 1010.0 " ;
       }
       std::cout << " at ievt=" << ievt ;
       std::cout << " nsum=" << nsum << " Sum E=" << esum << " SumPx=" << pxsum 
	<< " SumPy=" << pysum << " SumPz=" << pzsum << std::endl;
   }
   if( abs(sesum-esum) >1.0 ) {
       std::cout << " Evt=" << ievt << " Esum of gen.particles=" << esum << " Spring Esum=" << sesum << " difference is " << sesum-esum << std::endl;
   }

  }

}

//_________________________________________________________
void DrawHist()
{

}





