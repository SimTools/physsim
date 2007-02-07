//*----------------------------------------------------------------------------
//* ================
//*  Macro Plot.C
//* ================
//*
//* (Description)
//*	A sample macro to display ntuple contents.
//* (Usage)
//*	$ root -l
//*	......
//*	root [0] .L Plot.C
//*	root [1] Plot("jsf.root"); 
//*
//*	This will display distributions of some useful observables 
//*
//* (Update Record)
//*	2007/02/04	K.Fujii		Original version.
//*
//*----------------------------------------------------------------------------

TH2D *hMxMz      = new TH2D("hMxMz"    , "", 120,   90.,  150.,
                                             120,   60.,  120.);
TH2D *hCosXCosZ  = new TH2D("hCosXCosZ", "", 100,   -1.,   +1.,
                                             100,   -1.,   +1.); 
TH1D *hEvis      = new TH1D("hEvis"    , "", 120,    0.,  600.);
TH1D *hCosJh     = new TH1D("hCosJh"   , "", 100,   -1.,   +1.);
TH1D *hPhiJh     = new TH1D("hPhiJh"   , "", 180, -180., +180.);
TH1D *hCosAh     = new TH1D("hCosAh"   , "", 100,   -1.,   +1.);
TH1D *hPhiAh     = new TH1D("hPhiAh"   , "", 180, -180., +180.);

void Plot(Char_t *filen = "jsf.root")
{
  //gROOT->Reset();
  gSystem->Load("libS4Utils.so");
  gSystem->Load("libAnlib.so");

  Int_t nZoneX = 4;
  Int_t nZoneY = 3;

  TCanvas *c1 = new TCanvas("c1","",0,0, 300*nZoneX, 300*nZoneY);
  c1->SetHighLightColor(5);
  c1->SetFillColor(19);
  c1->Divide(nZoneX,nZoneY);

  gStyle->SetOptFit();

  cerr << filen << endl;
  TFile   *filep = new TFile(filen);
  TNtupleD *tup   = (TNtupleD *)gROOT->FindObject("hEvt");

  //--
  // Preselection
  //--
  TString sel("evis>0.");
  tup->Draw(">>elist",sel.Data(),"goff");

  TEventList *elist = static_cast<TEventList*>(gROOT->FindObject("elist"));
  Int_t nlist = elist->GetN();
  cerr << "Nevent = " << nlist << endl;

  //--
  // Loop over preselected events
  //--
  Double_t kToDeg = 180./TMath::Pi();
  for (Int_t i=0; i<nlist; i++) {
#if 0
        Int_t event = elist->GetEntry(i);
        cerr << "i = " << i << " event = " << event << endl;
#else
        Int_t event = i;
#endif
        Double_t ntracks;
        Double_t evis;
        Double_t pt;
        Double_t pl;
        Double_t ngams;
        Double_t njets;
        Double_t ycut;
        Double_t chi2;
        Double_t ea1;
        Double_t pa1x;
        Double_t pa1y;
        Double_t pa1z;
        Double_t ea2;
        Double_t pa2x;
        Double_t pa2y;
        Double_t pa2z;
        Double_t npj1;
        Double_t ej1;
        Double_t pj1x;
        Double_t pj1y;
        Double_t pj1z;
        Double_t npj2;
        Double_t ej2;
        Double_t pj2x;
        Double_t pj2y;
        Double_t pj2z;
        Double_t cosa1h;
        Double_t phia1h;
        Double_t cosa2h;
        Double_t phia2h;
        Double_t cosj1h;
        Double_t phij1h;
        Double_t cosj2h;
        Double_t phij2h;
        Double_t econe1;
        Double_t econe2;
        Double_t xmass;
        Double_t zmass;

        tup->SetBranchAddress("ntracks",&ntracks);
        tup->SetBranchAddress("evis",&evis);
        tup->SetBranchAddress("pt",&pt);
        tup->SetBranchAddress("pl",&pl);
        tup->SetBranchAddress("ngams",&ngams);
        tup->SetBranchAddress("njets",&njets);
        tup->SetBranchAddress("ycut",&ycut);
        tup->SetBranchAddress("chi2",&chi2);
        tup->SetBranchAddress("ea1",&ea1);
        tup->SetBranchAddress("pa1x",&pa1x);
        tup->SetBranchAddress("pa1y",&pa1y);
        tup->SetBranchAddress("pa1z",&pa1z);
        tup->SetBranchAddress("ea2",&ea2);
        tup->SetBranchAddress("pa2x",&pa2x);
        tup->SetBranchAddress("pa2y",&pa2y);
        tup->SetBranchAddress("pa2z",&pa2z);
        tup->SetBranchAddress("npj1",&npj1);
        tup->SetBranchAddress("ej1",&ej1);
        tup->SetBranchAddress("pj1x",&pj1x);
        tup->SetBranchAddress("pj1y",&pj1y);
        tup->SetBranchAddress("pj1z",&pj1z);
        tup->SetBranchAddress("npj2",&npj2);
        tup->SetBranchAddress("ej2",&ej2);
        tup->SetBranchAddress("pj2x",&pj2x);
        tup->SetBranchAddress("pj2y",&pj2y);
        tup->SetBranchAddress("pj2z",&pj2z);
        tup->SetBranchAddress("cosa1h",&cosa1h);
        tup->SetBranchAddress("phia1h",&phia1h);
        tup->SetBranchAddress("cosa2h",&cosa2h);
        tup->SetBranchAddress("phia2h",&phia2h);
        tup->SetBranchAddress("cosj1h",&cosj1h);
        tup->SetBranchAddress("phij1h",&phij1h);
        tup->SetBranchAddress("cosj2h",&cosj2h);
        tup->SetBranchAddress("phij2h",&phij2h);
        tup->SetBranchAddress("econe1",&econe1);
        tup->SetBranchAddress("econe2",&econe2);
        tup->SetBranchAddress("xmass",&xmass);
        tup->SetBranchAddress("zmass",&zmass);

        tup->GetEntry(event);
#if 0
  cerr << " pa1 = (" << ea1  << ", "
                     << pa1x << ", "
                     << pa1y << ", "
                     << pa1z << ") " << endl;
#endif
        ANL4DVector pa1(ea1, pa1x, pa1y, pa1z);
        ANL4DVector pa2(ea2, pa2x, pa2y, pa2z);
        ANL4DVector pj1(ej1, pj1x, pj1y, pj1z);
        ANL4DVector pj2(ej2, pj2x, pj2y, pj2z);
        ANL4DVector px = pa1 + pa2;
        ANL4DVector pz = pj1 + pj2;
#if 0 
  cerr << " xmass = " << px.GetMass() << " zmass = " << pz.GetMass() << endl;
#endif

        phij1h *= kToDeg;
        phij2h *= kToDeg;
        phia1h *= kToDeg;
        phia2h *= kToDeg;
        hMxMz    ->Fill(xmass, zmass, 1.);
        hCosXCosZ->Fill(px.CosTheta(), pz.CosTheta(), 1.);
        hEvis    ->Fill(evis  , 1.);
        hCosJh   ->Fill(cosj1h, 1.);
        hCosJh   ->Fill(cosj2h, 1.);
        hPhiJh   ->Fill(phij1h, 1.);
        hPhiJh   ->Fill(phij2h, 1.);
        hCosAh   ->Fill(cosa1h, 1.);
        hCosAh   ->Fill(cosa2h, 1.);
        hPhiAh   ->Fill(phia1h, 1.);
        hPhiAh   ->Fill(phia2h, 1.);
  }
  Int_t id = 0;
  id++; c1->cd(id); hStat->Draw();
  id++; c1->cd(id); hMxMz->Draw();
  id++; c1->cd(id); hMxMz->ProjectionX()->Draw();
  id++; c1->cd(id); hMxMz->ProjectionY()->Draw();
  id++; c1->cd(id); hEvis->Draw();
  id++; c1->cd(id); hCosXCosZ->Draw();
  id++; c1->cd(id); hCosXCosZ->ProjectionX();
                    hCosXCosZ_px->SetMinimum(0.);
                    hCosXCosZ_px->Draw();
  id++; c1->cd(id); hCosXCosZ->ProjectionY();
                    hCosXCosZ_py->SetMinimum(0.);
                    hCosXCosZ_py->Draw();
  id++; c1->cd(id); hCosJh->SetMinimum(0.); hCosJh->Draw();
  id++; c1->cd(id); hPhiJh->SetMinimum(0.); hPhiJh->Draw();
  id++; c1->cd(id); hCosAh->SetMinimum(0.); hCosAh->Draw();
  id++; c1->cd(id); hPhiAh->SetMinimum(0.); hPhiAh->Draw();
}
