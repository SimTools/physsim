#include <iostream>
#include <vector>

namespace std

string root_dir=string("root/");
string logfile_dir=string("log/");
string gprocinfo;
string ginfile;
string gprocid;
string gDBSdir("/home/ilc/miyamoto/DBS/runinfo/");

vector<string> gFigs;

// Plot Event based Ntuple
Int_t Plotfig_event()
{

  TFile *f=new TFile(ginfile.c_str());
  TNtuple *nte=(TNtuple*)f->Get("nte");
  TCanvas *cev1=new TCanvas("cev1","Event plot",1000,800);
  cev1->Divide(3,2);

  NTuplePlot(cev1, 1, 1, nte, "etot", "Total Energy (inc. #nu): ","Total Energy(GeV)");
  NTuplePlot(cev1, 2, 1, nte, "etot", "Total Energy (inc. #nu): ","Total Energy(GeV)","etot>340.0","etotc");
  NTuplePlot(cev1, 3, 1, nte, "evis", "Total Energy (w/o #nu): ","Energy(GeV)");
  NTuplePlot(cev1, 4, 1, nte, "pt", "Pt Sum: ","Pt(sum) (GeV)");
  NTuplePlot(cev1, 5, 1, nte, "pz", "Pz Sum: ","Pz(sum) (GeV)");

  NTuplePlot(cev1, 6, 0, nte, "thrust", "Thrust ","Thrust");
  PrintCanvas(cev1, "-event1");

  TCanvas *cev2=new TCanvas("cev2","Event plot2",1000,800);
  cev2->Divide(3,2);
  NTuplePlot(cev2,1,0,nte,"nchg","No. of Charged particles : ","No. of charged particles");
  NTuplePlot(cev2,2,0,nte,"ngam","No. of Gamma : ","No. of gamma");
  NTuplePlot(cev2,3,0,nte,"nhad","No. of Hadronic particles : ","No. of hadronic particles");
  NTuplePlot(cev2,4,0,nte,"echg","Total charged energy : ","Total charged energy (GeV)");
  NTuplePlot(cev2,5,0,nte,"egam","Total gamma energy : ","Total gamma energy (GeV)");
  NTuplePlot(cev2,6,0,nte,"ehad","Total hadronic energy : ","Total gamma energy (GeV)");
  PrintCanvas(cev2, "-event2");

// 2-Jets ntuple
  TCanvas *cj2=new TCanvas("cj2","2 Jets ", 1000, 800);
  cj2->Divide(3,2);
  NTuplePlot(cj2, 1, 0, nte, "csth", "cos#Theta(thrust) ","Thrust");
  NTuplePlot(cj2, 2,0,nte,"njets","No. jets (ycut=0.001) : ","No. of jets");

  TNtuple *nj2=(TNtuple*)f->Get("nj2");
  Int_t nj2ent=nj2->GetEntries();
  if( nj2ent > 0 ) { 
    NTuplePlot(cj2,3,0,nj2,"yc2j","Ycut for 2-jets :","ycut for forced 2-jets");
    NTuplePlot2Hist(cj2,4,0,nj2,"j0e","j1e","energy for 2-jets :","Jet energy (GeV)");	
    NTuplePlot2Hist(cj2,5,0,nj2,"j0cs","j1cs","cos#theta for 2-jets :","cos#theta(jet)");	
    NTuplePlot2D(cj2,6,nj2,"j0m","j1m","Jet mass for forced 2-jets : ","Mass(jet1)", "Mass(jet2)");
  }
  else {
    cj2->cd(3);
    string msg=gprocinfo;
    TLatex *tx2=new TLatex(0.2,0.5,msg.data());
    tx2->Draw();
    string msg2b=string(": There are no 2-jets ");
    TLatex *tx2b=new TLatex(0.2,0.4,msg2b.data());
    tx2b->Draw();
  }


  PrintCanvas(cj2, "-2jets");

// 4-Jets ntuple
  TCanvas *cj4=new TCanvas("cj4","4 Jets ", 1000, 800);
  TNtuple *nj4=(TNtuple*)f->Get("nj4");
  Int_t nj4ent=nj4->GetEntries();
  if( nj4ent > 0 ) {
    cj4->Divide(3,2);
    NTuplePlot(cj4,1,1,nj4,"yc4j","Ycut for 4-jets :","ycut for forced 4-jets");
    NTuplePlotAdd(cj4, 2, 1, nj4, "Jet energies(All jets)","E(Jet)",
     "j0e","j1e","j2e","j3e");
    NTuplePlotAdd(cj4, 3, 1, nj4, "Jet cos#theta (All jets)","cos#theta",
     "j0cs","j1cs","j2cs","j3cs");
    NTuplePlotAdd(cj4, 4, 1, nj4, "2jets mass(All comb)","M(J_{i}J_{k}) (GeV)",
     "m01","m02","m03","m12","m13","m23");
  }
  else {
    string msg=gprocinfo+string("  : There are no 4-jets ");
    TLatex *tx4=new TLatex(0.2,0.5,msg.data());
    tx4->Draw();
  }
  PrintCanvas(cj4, "-4jets");
}

// ======================================================================
void PrintCanvas(TCanvas *c, char *pref)
{
  string pngfile=logfile_dir+gprocid+string(pref)+string(".png");
  c->Print(pngfile.data());
  gFigs.push_back(gprocid+string(pref));

  string epsfile=logfile_dir+gprocid+string(pref)+string(".eps");
  c->Print(epsfile.data());
}

// ======================================================================
void NTuplePlot(TCanvas *cv, Int_t loc, Int_t logy, TNtuple *nt, char *vname, char *title, char *xaxis, char *pcut="" , char *hname="")
{
  cv->cd(loc);
  gPad->SetLogy(logy); 
  string hvname=string(vname);
  if( string(hname) != string("")) { hvname = string(hname); }
  string pltcmd=vname+string(" >>h_")+hvname;
  nt->Draw(pltcmd.data(),pcut);
  string hobj=string("h_")+hvname;
//  std::cerr << " hobj=" << hobj << std::endl;
  TH1 *h=(TH1*)gROOT->FindObject(hobj.data());
  string thistitle=string(title)+gprocinfo;
  if( h ) { 
	h->SetTitle(thistitle.data());
	h->GetXaxis()->SetTitle(xaxis);
	h->Draw();
  }
  else { 
    cerr << " Can not find " << thistitle << " histogram" << endl;
  }
}

// ======================================================================
void NTuplePlot2Hist(TCanvas *cv, Int_t loc, Int_t logy, TNtuple *nt, char *vname1, char* vname2, char *title, char *xaxis )
{
// Overlay two histogram
  cv->cd(loc);
  gPad->SetLogy(logy);
  string pltcmd=vname1+string(" >>h_")+vname1;
  nt->Draw(pltcmd.data());
  string hobj1=string("h_")+vname1;
  string pltcmd2=vname2+string(" >>h_")+vname2;
  nt->Draw(pltcmd2.data(),"","same");
  string hobj2=string("h_")+vname2;

  TH1 *h1=(TH1*)gROOT->FindObject(hobj1.data());
  TH1 *h2=(TH1*)gROOT->FindObject(hobj2.data());
  h1->Add(h2);
  string thistitle=string(title)+gprocinfo;
  h1->SetTitle(thistitle.data());
  h1->GetXaxis()->SetTitle(xaxis);
  h1->Draw();
//  h2->Draw("same");

}

// ================================================================================
void NTuplePlotAdd(TCanvas *cv, Int_t loc, Int_t logy, TNtuple *nt, char *title, char *xaxis ,
 char *vname1, char *vname2, char *vname3="", char *vname4="", char *vname5="", char *vname6="",
 char *vname7="", char *vname8="")
{
// Overlay two histogram
  cv->cd(loc);
  gPad->SetLogy(logy);

  vector<string> vnames;
  vnames.push_back(string(vname1));
  vnames.push_back(string(vname2));
  if( string(vname3) != string("") ) { vnames.push_back(string(vname3)); }  
  if( string(vname4) != string("") ) { vnames.push_back(string(vname4)); }  
  if( string(vname5) != string("") ) { vnames.push_back(string(vname5)); }  
  if( string(vname6) != string("") ) { vnames.push_back(string(vname6)); }  
  if( string(vname7) != string("") ) { vnames.push_back(string(vname7)); }  
  if( string(vname8) != string("") ) { vnames.push_back(string(vname8)); }  

  Int_t nvar=vnames.size();
  string pltcmd=string(vname1)+string(" >>ha_")+string(vname1);
  nt->Draw(pltcmd.data());
  string hobj1=string("ha_")+string(vname1);
  TH1 *h1=(TH1*)gROOT->FindObject(hobj1.data());
  for(Int_t i=1;i<nvar;i++) {
    string pltcmdi=vnames[i]+string(">>ha_")+vnames[i];
    nt->Draw(pltcmdi.data(),"","same");
  }
  for(Int_t i=1;i<nvar;i++) {
    string hobji=string("ha_")+vnames[i];
    TH1 *h1i=(TH1*)gROOT->FindObject(hobji.data());
    h1->Add(h1i);
  } 

  string thistitle=string(title)+gprocinfo;
  h1->SetTitle(thistitle.data());
  h1->GetXaxis()->SetTitle(xaxis);
  h1->Draw();

}


// ======================================================================
void NTuplePlot2D(TCanvas *cv, Int_t loc, TNtuple *nt, char *vname1, char *vname2,  char *title, char *xaxis, char *yaxis)
{
  cv->cd(loc);
  string pltcmd=vname1+string(":")+vname2+string(" >>h_")+vname1+vname2;
  nt->Draw(pltcmd.data());
  string hobj=string("h_")+vname1+vname2;
//  std::cerr << " hobj=" << hobj << std::endl;
  TH2 *h=(TH2*)gROOT->FindObject(hobj.data());
  string thistitle=string(title)+gprocinfo;
  h->SetTitle(thistitle.data());
  h->GetXaxis()->SetTitle(xaxis);
  h->GetYaxis()->SetTitle(yaxis);
  h->Draw();

}





// ======================================================================
// Plot Particle based Ntuple
Int_t Plotfig_part()
{

  TFile *f=new TFile(ginfile.c_str());
  TNtuple *ntp=(TNtuple*)f->Get("ntp");

  TCanvas *c1=new TCanvas("c1","Particle plot",1000,800);
  c1->Divide(3,2);


  NTuplePlot(c1, 1, 1, ntp, "id", "Particle ID : ","Particle ID");
  NTuplePlot(c1, 2, 1, ntp, "e", "Particle Energy : ","Energy(GeV)");
  NTuplePlot(c1, 3, 1, ntp, "pt", "Particle Pt : ","Pt(GeV)");
  NTuplePlot(c1, 4, 1, ntp, "cs", "Particle cos#theta : ","cos#theta");
  NTuplePlot(c1, 5, 1, ntp, "vr", "Start point R : ","Start point Radius");
  NTuplePlot(c1, 6, 1, ntp, "vz", "Start point Z : ","Start point ");

  PrintCanvas(c1,"-part");

}

// ==================================================================
void ShowTitle()
{
  string runinfo=gDBSdir+gprocid+string(".txt");
  ifstream ifs(runinfo.data());
  string line;
  TCanvas *ct=new TCanvas("ct","Title",1000,800);

  string titlemsg=string("Generator file information: ") + gprocinfo;
  TLatex *tx0=new TLatex(0.05,0.95,titlemsg.data());
  tx0->Draw();

  Double_t xpos=0.1;
  Double_t ypos=0.90;
  Double_t dy=0.03;
  while ( ifs && getline(ifs,line) ) {
    TLatex *txt=new TLatex(xpos,ypos,line.data());
    txt->SetTextSize(0.03); 
    txt->Draw();
    ypos-=dy;
  }
  PrintCanvas(ct,"-title");

}

// ==================================================================
// Int_t plotfig(const char* procid="w100689", const char* procname="ssbb_o")
Int_t plotfig(const char* procid="w100375", const char* procname="n1e1e1n1_o",
	      const char* fnpref="")
{

  JSFSteer *gJSF=JSFSteer::Instance();

  gDBSDir=std::string(gJSF->Env()->GetValue("RunInfoDir","/home/ilc/miyamoto/DBS/runinfo"));

//   ginfile=root_dir+string(procid)+string(".root");
  ginfile=root_dir+string(fnpref)+string(".root");
  gprocinfo=string(procid)+string(" ")+string(procname);
  gprocid=string(procid);

  ShowTitle();
  Plotfig_part();
  Plotfig_event();

#if 0
  Int_t nfigs=gFigs.size();
  for(Int_t i=0;i<nfigs;i++) {
    string ccmd=string("/usr/bin/convert ")+logfile_dir+gFigs[i]+string(".png ")
		+logfile_dir+gFigs[i]+string(".pdf");
    cerr << " ccmd=" << ccmd << endl;
    gSystem->Exec(ccmd.data());
  }
  string ccmd=string("/home/ilc/miyamoto/bin/pdftk ");
  for(Int_t i=0;i<nfigs;i++) {
    ccmd+= logfile_dir+gFigs[i]+string(".pdf ");
  }
  ccmd+=string(" output ")+logfile_dir+gprocid+string(".pdf");
  gSystem->Exec(ccmd.data());
  cerr << ccmd << " executed." << endl;
#endif
}


