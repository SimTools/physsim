//****************************************************
//*
//*  Sample UserAnalysis Script 
//*  
//****************************************************

TCanvas *cHist;
TDirectory *cDir;
TH1F *hNCDC;
TH1F *hNVTX;
TH1F *hNGen;
TNtupleD *hEvt;

//_________________________________________________________
void UserInitialize()
{
  //  This function is called at the begining of the job or when
  //  "reset hist" action is selected in the gui menu.
  //  This is used to define/reset histograms.

  if( !hNCDC ) delete hNCDC; 
  if( !hNVTX ) delete hNVTX; 
  if( !hNGen ) delete hNGen; 
  if( !hEvt  ) delete hEvt; 
  hNCDC=new TH1F("hNCDC","Number of CDC Tracks",100,0,100);
  hNVTX=new TH1F("hNVTX","Number of VTX Hits",100,0,200);
  hNGen=new TH1F("hNGen","Number of Generator Tracks",100,0,200);
  hEvt=new TNtupleD("hEvt", "", "helm:help:erhom:erhop:epim:epip");
  cDir=gDirectory;
}

//_________________________________________________________
void UserAnalysis()
{
  // This function is called when the processing of one event is completed.
  // Any data processing of the event can be performed in this function.
  // 

#if 1
  JSFSpring    *sp =(JSFSpring*)jsf->FindModule("FFSpring");
  JSFSpringBuf *spb=(JSFSpringBuf*)sp->EventBuf();
  Int_t nsps = spb->GetNpartons();
  //cerr << " --- Npartons = " << nsps << " -------- " << endl;
  ANL4DVector pcm;
  TClonesArray *partons = spb->GetPartons();
#if 1
  Int_t hel[2];
  Int_t np = 0;
#endif
  TIter next(partons);
  JSFSpringParton *parton;
  while ((parton = (JSFSpringParton *)next())) {
#if 0
    cerr << parton->GetSerial() << " "
         << " PID=" << parton->GetID()
         << " Q="   << parton->GetCharge()
         << " M="   << parton->GetMass()
         << " H="   << parton->GetHelicity()
         << " C="   << parton->GetColorID()
         << " S="   << parton->GetShowerInfo()
	 << " p=("  << parton->GetE() << ","
	            << parton->GetPx() << ","
	            << parton->GetPy() << ","
	            << parton->GetPz() << ")" << endl;
#endif
    if (parton->GetNDaughter()) continue;
    pcm += ANL4DVector(parton->GetPV());
#if 1
    if (np<2) hel[np++] = parton->GetHelicity();
#endif
  }
#if 0
  cerr << " SpringParton : ";
  pcm.DebugPrint();
#endif
#endif


  JSFSIMDSTBuf *sdb=(JSFSIMDSTBuf*)simdst->EventBuf();
  //  Accumulate information in the histogram
  hNCDC->Fill((Float_t)sdb->GetNCDCTracks());
  hNVTX->Fill((Float_t)sdb->GetNVTXHits());
  hNGen->Fill((Float_t)sdb->GetNGeneratorParticles());

  /*
  **  If these comments are removed, generator particle information
  **  are printed.
  */
  ANL4DVector qcm;
#if 1
  printf(" # Generator Particle is %d\n",sdb->GetNGeneratorParticles());
  TClonesArray *gen=sdb->GetGeneratorParticles();
  for(Int_t i=0;i<sdb->GetNGeneratorParticles();i++){
    JSFGeneratorParticle *g=gen->UncheckedAt(i);
    Int_t ndau=g->GetNDaughter();
#if 0
    g->ls();
#endif
#if 1
    Double_t erhom=0;
    Double_t erhop=0;
    Double_t epim=0;
    Double_t epip=0;
    Int_t pid=g->GetID();
    Int_t idau = g->GetFirstDaughter();
    if (pid==15 && ndau == 2) {
       JSFGeneratorParticle *gd=(JSFGeneratorParticle *)gen->UncheckedAt(idau);
       Int_t pidd=gd->GetID();
       if (pidd == -211) epim = gd->GetE();
       if (pidd == -213) {
         erhom = gd->GetE(); // rho
         Int_t idd = gd->GetFirstDaughter();
         JSFGeneratorParticle *gdd=(JSFGeneratorParticle *)gen->UncheckedAt(idd-1);
	 epim = gdd->GetE();
       }
    } else if (pid==-15 && ndau == 2) {
       JSFGeneratorParticle *gd=gen->UncheckedAt(idau);
       Int_t pidd=gd->GetID();
       if (pidd == 211) epip = gd->GetE();
       if (pidd == 213) {
         erhop = gd->GetE(); // rho
         Int_t idd = gd->GetFirstDaughter();
         JSFGeneratorParticle *gdd=(JSFGeneratorParticle *)gen->UncheckedAt(idd-1);
	 epip = gdd->GetE();
       }
    }
    hEvt->Fill(hel[0],hel[1],erhom,erhop,epim,epip);
#endif
    if( ndau != 0 ) continue;
    // printf(" ndau=%d\n",ndau);
    qcm += ANL4DVector(g->GetPV());
    //g->ls();
  }
#if 0
  cerr << " GeneratorParticle : ";
  qcm.DebugPrint();
#endif
#endif
}

//_________________________________________________________
void DrawHist()
{
  //  This function is called to draw histograms during the interactive 
  //  session.  Thus you can see the accumulation of the histogram
  //  interactively.  

  TDirectory *old=gDirectory;
  if( !cHist ) {
    cHist=new TCanvas("cHist","Canvas 1",100, 100, 800, 800);
    cHist->Divide(2,2);
  } 
  else {
    cHist->cd();
  }
  cHist->cd(1);
  hNCDC->Draw();
  cHist->cd(2);
  hNVTX->Draw();
  cHist->cd(3);
  hNGen->Draw();
  cHist->Update();

  old->cd();
}

//_________________________________________________________
void UserSetOptions()
{
  // This function is called only once, soon after jsf is started.
  // This function can be used to define parameters which is not 
  // defined in jsf.conf file.

}

//_________________________________________________________
void UserTerminate()
{
  // This function is called at the end of job.

}





