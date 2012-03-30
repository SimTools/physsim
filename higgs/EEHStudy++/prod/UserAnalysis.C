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

//_________________________________________________________
void UserInitialize()
{
  //  This function is called at the begining of the job or when
  //  "reset hist" action is selected in the gui menu.
  //  This is used to define/reset histograms.

  if( !hNCDC ) delete hNCDC; 
  if( !hNVTX ) delete hNVTX; 
  if( !hNGen ) delete hNGen; 
  hNCDC=new TH1F("hNCDC","Number of CDC Tracks",100,0,100);
  hNVTX=new TH1F("hNVTX","Number of VTX Hits",100,0,200);
  hNGen=new TH1F("hNGen","Number of Generator Tracks",100,0,200);
  cDir=gDirectory;
}

//_________________________________________________________
void UserAnalysis()
{
  // This function is called when the processing of one event is completed.
  // Any data processing of the event can be performed in this function.
  // 

#if 1
  JSFSpring    *sp =(JSFSpring*)jsf->FindModule("EEHSpring");
  JSFSpringBuf *spb=(JSFSpringBuf*)sp->EventBuf();
  Int_t nsps = spb->GetNpartons();
  cerr << " --- Npartons = " << nsps << " -------- " << endl;
  ANL4DVector pcm;
  TClonesArray *partons = spb->GetPartons();
  TIter next(partons);
  JSFSpringParton *parton;
  while ((parton = (JSFSpringParton *)next())) {
    cerr << parton->GetSerial() << " "
         << " PID=" << parton->GetID()
         << " Q="   << parton->GetCharge()
         << " M="   << parton->GetMass()
         << " M="   << parton->GetColorID()
         << " M="   << parton->GetShowerInfo()
	 << " p=("  << parton->GetE() << ","
	            << parton->GetPx() << ","
	            << parton->GetPy() << ","
	            << parton->GetPz() << ")" << endl;
    if (parton->GetNDaughter()) continue;
    pcm += ANL4DVector(parton->GetPV());
  }
  cerr << " SpringParton : ";
  pcm.DebugPrint();
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
    if( ndau != 0 ) continue;
    // printf(" ndau=%d\n",ndau);
    qcm += ANL4DVector(g->GetPV());
    // g->ls();
  }
  cerr << " GeneratorParticle : ";
  qcm.DebugPrint();
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
