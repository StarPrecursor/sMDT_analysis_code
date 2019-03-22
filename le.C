#include "Riostream.h"
#include"TFile.h"
#include"TH2.h"
#include"math.h"
#include"cmath"
void le(){
Int_t i;
Double_t x,y;
Int_t x1,y1,z1;
Double_t k,o,d;
ifstream in("./fileinfo4.txt");
char filename[50];
in>>filename;
ifstream in1(filename);
TFile *f=new TFile("muondist.root");
TCanvas *c1=new TCanvas("c1","hittime");
c1->cd();
TH2D* h=new TH2D("h","legendre space",1000,0,3.2,1000,-44,44);
while(in1>>x>>y){
h->Fill(x,y);
}
Int_t MaxBin=h->GetMaximumBin();
h->GetBinXYZ(MaxBin,x1,y1,z1);
printf("maxbin's position is (%d,%d,%d)\n",x1,y1,z1);
printf("maxbin is %d\n",MaxBin);
o=x1/1000*3.2;
k=-1/tan(o);
d=y1/1000*88-44;
d=d/sin(o);
cout<<"Track's function is"<<"y="<<k<<"*x+"<<d;



}
