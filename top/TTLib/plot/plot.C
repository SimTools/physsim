{
const int kNpt = 1000;
double xdat[kNpt], ydat[10][kNpt];
ifstream in("sigma.dat");
int npt = 0;

const int kMax = 1024;
char dummy[kMax];
while (!in.eof()) {
   in.getline(dummy, kMax);
   stringstream ins(dummy);
   if (dummy[0] == '(') continue;
   cerr << dummy << endl;
   ins >> xdat[npt]
       >> ydat[0][npt]
       >> ydat[1][npt]
       >> ydat[2][npt]
       >> ydat[3][npt]
       >> ydat[4][npt]
       >> ydat[5][npt]
       >> ydat[6][npt]
       >> ydat[7][npt];
   npt++;
}
TH2D h("h","",10,330.,400.,10,0.,1500.);
h.Draw();
h.SetStats(0);
TGraph *gp[10];
for (int i=0; i<8; i++) {
   gp[i] = new TGraph(npt-1,xdat,&ydat[i][0]);
   gp[i]->SetLineWidth(2);
   if (i != 3) gp[i]->Draw("c same");
}
}
