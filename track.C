#include "TGraph.h"
#include "TCanvas.h"
#include "Riostream.h"

void track() {
    Int_t n=2;
    Int_t x[2]={0,40};
    Int_t y[2]={0,40};
    Int_t x1[2]={0,1};
    Int_t y1[2]={0,1};
    Int_t t[3];
    Int_t i2;
    Double_t dt;
    Double_t r1;

    TFile *f = new TFile("trackdemo.root");
    TCanvas *c1 = new TCanvas("c1", "track", 900, 900);
    c1->cd();

    TGraph *gr = new TGraph(n, x, y);
    gr->SetLineColor(1);
    gr->SetName("gr");
    gr->Draw("AP");

    // prepare tubes distribution graph without fill color
    Int_t i, j, k;
    TEllipse *c[4][4][6];
    Double_t r = 1.46 * sqrt(3);
    for (i = 0; i < 4; ++i) {
        if ((i == 1) || (i == 0)) {
            for(j = 0; j < 4; ++j) {
                if (j == 0) {
                    for (k = 0; k < 6; ++k) {
                        c[i][j][k] = new TEllipse(2.92 + 2.92 * k + 17.52 * (i),
                                                  1.46, 1.46, 1.46);
                        c[i][j][k]->Draw();
                    }
                }
                if (j == 1) {
                    for (k = 0; k < 6; ++k) {
                        c[i][j][k] = new TEllipse(1.46 + 2.92 * k + 17.52 * (i),
                                                  1.46 + r, 1.46, 1.46);
                        c[i][j][k]->Draw();
                    }
                }
                if (j == 2) {
                    for (k = 0; k <6 ; ++k) {
                        c[i][j][k] = new TEllipse(2.92 + 2.92 * k + 17.52 * (i),
                                                  1.46 + 2 * r, 1.46, 1.46);
                        c[i][j][k]->Draw();
                    }
                }
                if(j == 3) {
                    for (k = 0; k < 6; ++k) {
                        c[i][j][k] = new TEllipse(1.46 + 2.92 * k + 17.52 * (i),
                                                  1.46 + 3 * r, 1.46, 1.46);
                        c[i][j][k]->Draw();
                    }
                }
            }
        }
        if ((i == 2) || (i == 3)) {
            for (j = 0; j < 4; ++j) {
                if (j == 0) {
                    for (k = 0; k < 6; ++k) {
                        c[i][j][k] = new TEllipse(2.92 + 2.92 * k + 17.52 * (i-2), 1.46 + 27, 1.46, 1.46);
                        c[i][j][k]->Draw();
                    }
                }
                if (j == 1) {
                    for (k = 0; k < 6; ++k) {
                        c[i][j][k] = new TEllipse(1.46 + 2.92 * k + 17.52 * (i-2), 1.46 + r + 27, 1.46, 1.46);
                        c[i][j][k]->Draw();
                    }
                }
                if (j == 2) {
                    for (k = 0; k < 6; ++k) {
                        c[i][j][k] = new TEllipse(2.92 + 2.92 * k + 17.52 * (i-2), 1.46 + 2 * r + 27, 1.46, 1.46);
                        c[i][j][k]->Draw();
                    }
                }
                if (j == 3) {
                    for (k = 0; k < 6; ++k) {
                        c[i][j][k] = new TEllipse(1.46 + 2.92 * k + 17.52 * (i-2), 1.46 + 3 * r + 27, 1.46, 1.46);
                        c[i][j][k]->Draw();
                    }
                }
            }
        }
    }

    // fill the tube with color in hit position to show the track
    ifstream in("./fileinfo2.txt");
    char filename[50];
    in >> filename;
    ifstream in1(filename);
    while (in1 >> t[0] >> t[1] >> t[2] >> dt) {
        if (t[0] == 0) {
            i2 = 1;
        }
        if(t[0] == 1) {
            i2 = 3;
        }
        if(t[0] == 2) {
            i2 = 0;
        }
        if(t[0] == 3) {
            i2 = 2;
        }
        c[i2][t[1]][t[2]-1]->SetFillColor(1);
        cout << " " << i2 << " " << t[1] << " " << t[2]-1 << endl;
    }
}
