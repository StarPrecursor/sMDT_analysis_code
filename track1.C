#include "TGraph.h"
#include "TCanvas.h"
#include "Riostream.h"
#include "TF1.h"
void track1(){
TFile *f=new TFile("trackdemo.root");
TCanvas *c1=new TCanvas("c1","track",900,900);
c1->cd();
Int_t n=2;
Int_t x[2]={0,40};
Int_t y[2]={0,40};
    Int_t x1[2]={0,1};
    Int_t y1[2]={0,1};
Int_t t[3];
Int_t i2;
Int_t i3;
Double_t dt;
Double_t xx[4][4][6][2];
Double_t r1;
Double_t x2,y2;
TGraph *gr = new TGraph(n,x,y);
    gr->SetLineColor(1);
    gr->SetName("gr");
    gr->SetTitle("event id:62651(run213785)");
gr->Draw("AP");


    Int_t i,j,k;
    TEllipse *c[4][4][6];
    TEllipse *c3[16];
    Double_t r=1.46*sqrt(3);
    for(i=0;i<4;++i){
        if((i==1)||(i==0)){
            for(j=0;j<4;++j){
                if(j==0){
                    for(k=0;k<6;++k){
                        c[i][j][k]=new TEllipse(2.92+2.92*k+17.52*(i),1.46,1.46,1.46);
                        xx[i][j][k][0]=2.92+2.92*k+17.52*(i);
                        xx[i][j][k][1]=1.46;
                        c[i][j][k]->Draw();
                    }
                }

                if(j==1){
                    for(k=0;k<6;++k){
                        c[i][j][k]=new TEllipse(1.46+2.92*k+17.52*(i),1.46+r,1.46,1.46);
                        xx[i][j][k][0]=1.46+2.92*k+17.52*(i);
                        xx[i][j][k][1]=1.46+r;
                        c[i][j][k]->Draw();
                    }
                }
                if(j==2){
                    for(k=0;k<6;++k){
                        c[i][j][k]=new TEllipse(2.92+2.92*k+17.52*(i),1.46+2*r,1.46,1.46);
                        xx[i][j][k][0]=2.92+2.92*k+17.52*(i);
                        xx[i][j][k][1]=1.46+2*r;
                        c[i][j][k]->Draw();
                    }
                }
                if(j==3){
                    for(k=0;k<6;++k){
                        c[i][j][k]=new TEllipse(1.46+2.92*k+17.52*(i),1.46+3*r,1.46,1.46);
                        xx[i][j][k][0]=1.46+2.92*k+17.52*i;
                        xx[i][j][k][1]=1.46+3*r;
                        c[i][j][k]->Draw();
                    }
                }
                }
            }
        if((i==2)||(i==3)){
            for(j=0;j<4;++j){
                if(j==0){
                    for(k=0;k<6;++k){
                        c[i][j][k]=new TEllipse(2.92+2.92*k+17.52*(i-2),1.46+27,1.46,1.46);
                        xx[i][j][k][0]=2.92+2.92*k+17.52*(i-2);
                        xx[i][j][k][1]=1.46+27;
                        c[i][j][k]->Draw();
                    }
                }
                if(j==1){
                    for(k=0;k<6;++k){
                        c[i][j][k]=new TEllipse(1.46+2.92*k+17.52*(i-2),1.46+r+27,1.46,1.46);
                        xx[i][j][k][0]=1.46+2.92*k+17.52*(i-2);
                        xx[i][j][k][1]=1.46+r+27;
                        c[i][j][k]->Draw();
                    }
                }
                if(j==2){
                    for(k=0;k<6;++k){
                        c[i][j][k]=new TEllipse(2.92+2.92*k+17.52*(i-2),1.46+2*r+27,1.46,1.46);
                        xx[i][j][k][0]=2.92+2.92*k+17.52*(i-2);
                        xx[i][j][k][1]=1.46+2*r+27;
                        c[i][j][k]->Draw();
                    }
                }
                if(j==3){
                    for(k=0;k<6;++k){
                        c[i][j][k]=new TEllipse(1.46+2.92*k+17.52*(i-2),1.46+3*r+27,1.46,1.46);
                        xx[i][j][k][0]=1.46+2.92*k+17.52*(i-2);
                        xx[i][j][k][1]=1.46+3*r+27;
                        c[i][j][k]->Draw();
                    }
                }

            }
        }
    }

    ifstream in("./fileinfo3.txt");
    char filename[50];
    in>>filename;
    ifstream in1(filename);
    i3=0;
    while(in1>>t[0]>>t[1]>>t[2]>>dt>>r1){
        if(t[0]==0){
            i2=1;
        }
        if(t[0]==1){
            i2=3;
        }
        if(t[0]==2){
            i2=0;
        }
        if(t[0]==3){
            i2=2;
        }
        x2=xx[i2][t[1]][t[2]-1][0];
        y2=xx[i2][t[1]][t[2]-1][1];
        c3[i3]=new TEllipse(x2,y2,r1,r1);
        c3[i3]->Draw();
        c3[i3]->SetLineColor(3);
        ++i3;
        cout<<i3<<" "<<i2<<" "<<t[1]<<" "<<t[2]-1<<endl;

    }
    TF1 *fa1=new TF1("fa1","-3.922862635*x+116.8505531",0,40);
    fa1->Draw("SAME");






}
