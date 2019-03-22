#include <fstream>
#include <stdio.h>
#include <TMath.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TPad.h>
#include <TH1.h>
#include <TH2F.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TGaxis.h>

// Global variables can be set by the program including resiplot.h 
const char chend[4]="png";     // File extension for plots
char filename[255];            // prepended to plot file names
char plottitle[255];           // prepended to plot title

// ResiFitData is used to save fit parameters returned by various fit functions.
// To save ResiFitData to a file include statsfile.h and write out a stats file.
const int FITDATAMAX=20;             //maximum number of possible fit items  
Double_t ResiFitData[FITDATAMAX];    //filled by DrawGaus, DrawResi, DrawResiML, Drawh1h1, Drawh1h1h1
Double_t ResiFitData2[FITDATAMAX];   //filled by DrawResiML

//    DrawResiParams = parameters for gaussian fit, used by DrawResi, DrawResiML
//     DrawResiParams[0] => scale factor for residuals 
//                 set to 1000 if residual plots are in mm.
//     DrawResiParams[1] = 0 fit full range of plot 
//               > 0 number of RMSs for limit of DrawResi fit
//               < 0 +-ranges of fit [mm]  
Double_t DrawResiParams[2] = { 1., 5. };  

//single gaussian fit of histogram
int DrawGaus( TH1 *hp );                
TH1 *DrawGaus( const char *hname );

//check for 2 gaussians separated by 25ns in T0refit plot
TH1 *DrawT0refit( const char *hname );  
int  DrawT0refit(        TH1 *hp );

// Double gaussian fit, 2 gaussians constrained to same mean
TH1 *DrawResi(const char *hname );
Double_t  DrawResi(       TH1 *hp );

// Double gaussian fit, 2 independent gaussians 
TH1 *DrawResiInd(const char *hname);
int  DrawResiInd(       TH1 *hp   );
//  Gaussian functions for ROOT fits
Double_t DoubleG(Double_t *x,Double_t *par);
Double_t G1     (Double_t *x,Double_t *par);
Double_t G2     (Double_t *x,Double_t *par);

// Version for fitting 2 residual distributions and putting on same plot (e.g. to compare MLs)
int DrawResiML(const char *hname1, const char *hname2, const char *pname="", const char *fname="" );
int DrawResiML(          TH1 *hp1,           TH1 *hp2, const char *pname="", const char *fname="" );

// This just does the double gaussion fit, does not draw plot (used by DrawResiML)
int FitDoubleGaussian(TH1 *hp, Double_t values[FITDATAMAX]);

// print TH1F histogram
void Drawh1( const char *hname, const char *opt="", const char *stat="eMR");

// print 2 TH1F histogram2 overlaid (e.g. raw & seg or ML1 & ML2)
void Drawh1h1(const char *hname1, const char *hname2, const int iflag=1, const char *pname="", const char *fname="");

// print 3 TH1F histogram2 overlaid   
void Drawh1h1h1(const char *hname1, const char *hname2, const char *hname3, const char *pname="" );
void Drawh1h1h1(       TH1 *h1,            TH1 *h2,     const char *hname3, const char *pname="" );

// print TH2F histogram to file
TH2F* Drawh2(const char *hname, const char *chopt="", const char *xlab="", const char *ylab="", const char *chstat="" );
void  Drawh2(      TH2F *hp   , const char *chopt="", const char *xlab="", const char *ylab="", const char *chstat="" );

// utility program to make a plot file from whatever has just been plotted on the current canvas
void makeplot(const char *pname="plot" );

/************************************************************
*  DrawResiInd => Do double gaussian fit to a residual distribution, gaussians constrained to same mean
*   Parameters:
*     hp     => Pointer to histogram
*     DrawResiParams[0] => scale factor for residuals 
*                 set to 1000 if residual plots are in mm.
*     DrawResiParams[1] = 0 fit full range of plot 
*               > 0 number of RMSs for limit of DrawResi fit
*               < 0 +-ranges of fit [mm]  
*     ResiFitData => Contains fit parameters, see below
************************************************************/
TH1 *DrawResi(const char *hname ) {
//  By using TH1 this works for TH1F, TH1D, etc
  TH1 *hp = (TH1 *) gDirectory->Get(hname);
  if( hp == NULL ) {
    printf("DrawResi ERROR: No histogram %s\n",hname);
    return NULL;
  }
  DrawResi( hp );
  return hp;
}

/************************************************************
*  DrawResi => Do double gaussian fit to a residual distribution
*  Version fiting 2 gaussians constrained to have the same mean
*  ResiFitData[0] = Number of entries
*  ResiFitData[1] = Histogram mean
*  ResiFitData[2] = Error of mean = RMS/sqrt(entries)
*  ResiFitData[3] = RMS
*  ResiFitData[4] = Gaussian mean   (G1+G2 common mean)
*  ResiFitData[5] = Gaussian mean error  
*  ResiFitData[6] = Gaussian1 sigma  (narrow) 
*  ResiFitData[7] = Gaussian1 sigma error  (narrow)
*  ResiFitData[8] = Gaussian2 sigma  (wide)
*  ResiFitData[9] = Gaussian2 sigma error  (wide)
*  ResiFitData[10] = Ratio of area of gaussians (narrow/wide)
*  ResiFitData[11] = SigmaHM (==FWHM/2.3548 of the combined G1+G2 function)
************************************************************/
Double_t DrawResi(TH1 *hp) {

// Zero ResiFitData array
  for( int i=0; i<FITDATAMAX; i++ ) ResiFitData[i] = 0.;

//  Do not plot if no entries or RMS==0 (entries are all underflow/overflow)
  if( hp == NULL ) {
    printf("DrawResi WARNING: histogram is NULL\n");
    return -1;
  }

//  Do not plot if no entries or RMS==0 (entries are all underflow/overflow)
  if( hp->GetEntries() == 0 || hp->GetRMS() == 0 ) {
    printf("DrawResi WARNING: histogram %s no entries, skipping plot\n", hp->GetName());
    return -1;
  }

  printf("DrawResi: Fitting %s\n",hp->GetName());

// Rename histogram
  if( plottitle[0] != '\0' ) {
    hp->SetTitle(Form("%s %s",plottitle,hp->GetTitle()));
  }

  gStyle->SetOptStat("emr");
  gStyle->SetOptFit(0);
  //  gStyle->SetTitleX(0.55);
  //  gStyle->SetTitleY(1.05);
  //  gStyle->SetTitleH(0.2);
  //  gStyle->SetTitleW(0.75);
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.35);
  gStyle->SetStatY(0.9);

//  Set initial guesses of parameters
  Int_t    imax = hp->GetMaximumBin();
  Double_t mean = hp->GetMean();
  Float_t  ampl = (hp->GetBinContent(imax-1)+hp->GetBinContent(imax)+hp->GetBinContent(imax+1))/3.;
  Float_t  rms  = hp->GetRMS();

// Set fit limits according to param[1]
  Double_t xmin=hp->GetXaxis()->GetBinLowEdge(1);
  Double_t xmax=hp->GetXaxis()->GetBinUpEdge( hp->GetXaxis()->GetNbins() );
  if(      DrawResiParams[1] > 0. ) {
    xmin = mean - DrawResiParams[1]*rms;
    xmax = mean + DrawResiParams[1]*rms;
  } else if( DrawResiParams[1] < 0. ) {
    xmin =  DrawResiParams[1];
    xmax = -DrawResiParams[1];
  }

  TF1 *fun3 = new TF1("fun3",DoubleG,xmin,xmax,5);
  fun3->SetParNames("Const1","Sigma1","Mean","Const2","Sigma2");
  fun3->SetParameters(ampl,rms/4,mean,ampl/10.,rms*2); // constant,mean,sigma

  hp->Fit("fun3","RQ");

// Check results of fit.  If failed try again with restricted range
/*   if( fun3->GetParError(1) > 2.*fun3->GetParameter(1)  || */
/*       fun3->GetParError(4) > 2.*fun3->GetParameter(4) ) { */
/*     fun3->SetParameters(ampl,mean,rms/4,ampl/10.,mean,rms*2);    //constant,mean,sigma */
/*     hp->Fit("fun3","","",mean-2*rms,mean+2*rms); */
/*   } */
//  Determine SigmaHM.  Do binary search to find half max points on both sides
  Double_t Xmax = fun3->GetMaximumX(mean-rms,mean+rms);
  Double_t Hmax = fun3->Eval(Xmax)/2;   //GetMaximum(mean-rms,mean+rms);
  Double_t X1 = Xmax-2*rms;
  Double_t X2 = Xmax;
  Double_t Xlo = (X1+X2)/2.;
  for( int ii=0; ii<25; ii++ ) {
    if( fun3->Eval(Xlo) < Hmax ) X1 = Xlo;
    else                         X2 = Xlo;
    Xlo = (X1+X2)/2.;
  }
  X1 = Xmax;
  X2 = Xmax+2*rms;
  Double_t Xhi = (X1+X2)/2.;
  for( int ii=0; ii<25; ii++ ) {
    if( fun3->Eval(Xhi) < Hmax ) X2 = Xhi;
    else                         X1 = Xhi;
    Xhi = (X1+X2)/2.;
  }
  Double_t sigmaHM = (Xhi-Xlo)/2.3548*DrawResiParams[0];
  //  Xmax *= DrawResiParams[0];
  //  printf("Xmax=%lf Hmax=%lf Xlo=%lf XloF=%lf\n",Xmax,Hmax,Xlo,fun3->Eval(Xlo));
  //  printf("Xhi=%lf XhiF=%lf sigma=%lf\n",Xhi,fun3->Eval(Xhi),sigmaHM);
//  Draw fit results on plot
//  TPaveText *sigpt = new TPaveText(0.1,0.5,0.38,0.9,"NDC");
  TPaveText *sigpt = new TPaveText(0.1,0.4,0.38,0.9,"NDC");
  sigpt->SetTextAlign(12);
  Float_t sig1   = TMath::Abs(fun3->GetParameter(1))*DrawResiParams[0];
  Float_t sig2   = TMath::Abs(fun3->GetParameter(4))*DrawResiParams[0];
  Float_t mean1  = fun3->GetParameter(2)*DrawResiParams[0];
  Float_t const1 = fun3->GetParameter(0);
  Float_t const2 = fun3->GetParameter(3);
  Float_t ratio  = 0.;
  if( sig1 < sig2 ) {
    sigpt->AddText(Form("#sigma_{1}=%.1lf #mum",sig1));
    sigpt->AddText(Form("#sigma_{2}=%.1lf #mum",sig2));
    sigpt->AddText(Form("#mu=%.1lf #mum",mean1));
    if( sig2 > 0. && const2 > 0. ) {
      ratio = sig1*const1/sig2/const2;
      sigpt->AddText(Form("A_{1}/A_{2}=%.3lf",ratio));
    } 
    ResiFitData[6] = sig1;
    ResiFitData[7] = fun3->GetParError(1)*DrawResiParams[0];
    ResiFitData[8] = sig2;
    ResiFitData[9] = fun3->GetParError(4)*DrawResiParams[0];
 } else {
    sigpt->AddText(Form("#sigma_{1}=%.1lf #mum",sig2));
    sigpt->AddText(Form("#sigma_{2}=%.1lf #mum",sig1));
    sigpt->AddText(Form("#mu=%.1lf #mum",mean1));
    if( sig1 > 0. && const1 > 0. ) {
      ratio = sig2*const2/sig1/const1;
      sigpt->AddText(Form("A_{1}/A_{2}=%.3lf",ratio));
    }
    ResiFitData[6] = sig2;
    ResiFitData[7] = fun3->GetParError(4)*DrawResiParams[0];
    ResiFitData[8] = sig1;
    ResiFitData[9] = fun3->GetParError(1)*DrawResiParams[0];
  }
  ResiFitData[0] = hp->GetEntries();
  ResiFitData[1] = hp->GetMean()*DrawResiParams[0];
  ResiFitData[3] = hp->GetRMS()*DrawResiParams[0];
  if( ResiFitData[0] > 0. ) ResiFitData[2] = ResiFitData[3]/TMath::Sqrt(ResiFitData[0]);
  ResiFitData[4]  = mean1;
  ResiFitData[5]  = fun3->GetParError(2)*DrawResiParams[0];
  ResiFitData[10] = ratio;
  ResiFitData[11] = sigmaHM;

  sigpt->AddText(Form("#sigma_{#epsilon}=%.1lf #mum",sigmaHM));
  sigpt->Draw();

  makeplot(hp->GetName());

  return sigmaHM;
} //end DrawResi ==========================================================

/************************************************************
*  DrawResiInd => Do double gaussian fit to a residual distribution, independent gaussians
*   Parameters:
*     hp     => Pointer to histogram
*     DrawResiParams[0] => scale factor for residuals 
*                 set to 1000 if residual plots are in mm.
*     DrawResiParams[1] = 0 fit full range of plot 
*               > 0 number of RMSs for limit of DrawResi fit
*               < 0 +-ranges of fit [mm]  
*  ResiFitData => 2 versions, see below
************************************************************/
TH1 *DrawResiInd(const char *hname ) {
//  By using TH1 this works for TH1F, TH1D, etc
  TH1 *hp = (TH1 *) gDirectory->Get(hname);
  if( hp == NULL ) {
    printf("DrawResiInd ERROR: No histogram %s\n",hname);
    return NULL;
  }
  DrawResiInd( hp );
  return hp;
}

/************************************************************
*  DrawResiInd => Do double gaussian fit to a residual distribution
*  Version with two independent gaussians (NO FINAL dummy integer)
*  ResiFitData[0] = Number of entries
*  ResiFitData[1] = Histogram mean
*  ResiFitData[2] = RMS
*  ResiFitData[3] = Gaussian1 mean   (narrow)
*  ResiFitData[4] = Gaussian1 sigma  (narrow) 
*  ResiFitData[5] = Gaussian2 mean   (wide)
*  ResiFitData[6] = Gaussian2 sigma  (wide)
*  ResiFitData[7] = Ratio of area of gaussians (narrow/wide)
*  ResiFitData[8] = Error of mean = RMS/sqrt(entries)
*  ResiFitData[9] = Gaussian1 mean error  (narrow)
*  ResiFitData[10] = Gaussian2 mean error  (wide)
*  ResiFitData[11] = Gaussian1 sigma error  (narrow)
*  ResiFitData[12] = Gaussian2 sigma error  (wide)
*  ResiFitData[13] = X value of maximum of combined function G1+G2
*  ResiFitData[14] = SigmaHM (==FWHM/2.3548 of the combined G1+G2 function)
*  ResiFitData[15] = X value of maximum of combined function G1+G2 error
*  ResiFitData[16] = SigmaHM error
************************************************************/
int DrawResiInd(TH1 *hp) {

// Zero ResiFitData array
  for( int i=0; i<FITDATAMAX; i++ ) ResiFitData[i] = 0.;

//  Do not plot if no entries or RMS==0 (entries are all underflow/overflow)
  if( hp == NULL ) {
    printf("DrawResiInd WARNING: histogram is NULL\n");
    return -1;
  }

//  Do not plot if no entries or RMS==0 (entries are all underflow/overflow)
  if( hp->GetEntries() == 0 || hp->GetRMS() == 0 ) {
    printf("DrawResiInd WARNING: histogram %s no entries, skipping plot\n",hp->GetName());
    return -1;
  }

  printf("DrawResiInd: Fitting %s\n",hp->GetName());

// Rename histogram
  if( plottitle[0] != '\0' ) {
    hp->SetTitle(Form("%s %s",plottitle,hp->GetTitle()));
  }

  //  printf("Create TF1\n");
  TF1 *fun3=new TF1("fun3","gaus(0)+gaus(3)"); 
  fun3->SetParNames("Const1","Mean1","Sigma1","Const2","Mean2","Sigma2");
//  Set initial guesses of parameters
  Int_t       i = hp->GetMaximumBin();
  Double_t mean = hp->GetMean();
  Float_t  ampl = (hp->GetBinContent(i-1) + hp->GetBinContent(i) + hp->GetBinContent(i+1))/3.;
  Float_t  rms  = hp->GetRMS();
  fun3->SetParameters(ampl,mean,rms/4,ampl/10.,mean,rms*2);    //constant,mean,sigma

  gStyle->SetOptStat("emr");
  gStyle->SetOptFit(0);
  gStyle->SetStatY(0.9);
  //  gStyle->SetTitleH(0.1);
  //  gStyle->SetTitleW(0.5);
  gStyle->SetStatW(0.38);
  gStyle->SetStatH(0.35);

// Set fit limits according to param[1]
  if(      DrawResiParams[1] > 0. ) hp->Fit("fun3","","",mean-DrawResiParams[1]*rms,mean+DrawResiParams[1]*rms);
  else if( DrawResiParams[1] < 0. ) hp->Fit("fun3","","",DrawResiParams[1],-DrawResiParams[1]);
  else                      hp->Fit("fun3");

// Check results of fit.  If failed try again with restricted range
/*   if( fun3->GetParError(1) > 2.*fun3->GetParameter(1)  || */
/*       fun3->GetParError(4) > 2.*fun3->GetParameter(4) ) { */
/*     fun3->SetParameters(ampl,mean,rms/4,ampl/10.,mean,rms*2);    //constant,mean,sigma */
/*     hp->Fit("fun3","","",mean-2*rms,mean+2*rms); */
/*   } */
//  Determine SigmaHM.  Do binary search to find half max points on both sides
  Double_t Xmax = fun3->GetMaximumX(mean-rms,mean+rms);
  Double_t Hmax = fun3->Eval(Xmax)/2;   //GetMaximum(mean-rms,mean+rms);
  Double_t X1 = Xmax-2*rms;
  Double_t X2 = Xmax;
  Double_t Xlo = (X1+X2)/2.;
  for( int ii=0; ii<25; ii++ ) {
    if( fun3->Eval(Xlo) < Hmax ) X1 = Xlo;
    else                         X2 = Xlo;
    Xlo = (X1+X2)/2.;
  }
  X1 = Xmax;
  X2 = Xmax+2*rms;
  Double_t Xhi = (X1+X2)/2.;
  for( int ii=0; ii<25; ii++ ) {
    if( fun3->Eval(Xhi) < Hmax ) X2 = Xhi;
    else                         X1 = Xhi;
    Xhi = (X1+X2)/2.;
  }
  Xmax *= DrawResiParams[0];
  Double_t sigmaHM = (Xhi-Xlo)/2.3548*DrawResiParams[0];
  //  printf("Xmax=%lf Hmax=%lf Xlo=%lf XloF=%lf\n",Xmax,Hmax,Xlo,fun3->Eval(Xlo));
  //  printf("Xhi=%lf XhiF=%lf sigma=%lf\n",Xhi,fun3->Eval(Xhi),sigmaHM);
//  Draw fit results on plot
//  TPaveText *sigpt = new TPaveText(0.1,0.5,0.38,0.9,"NDC");
  TPaveText *sigpt = new TPaveText(0.1,0.2,0.38,0.9,"NDC");
  sigpt->SetTextAlign(12);
  Float_t sig1   = TMath::Abs(fun3->GetParameter(2))*DrawResiParams[0];
  Float_t sig2   = TMath::Abs(fun3->GetParameter(5))*DrawResiParams[0];
  Float_t mean1  = fun3->GetParameter(1)*DrawResiParams[0];
  Float_t mean2  = fun3->GetParameter(4)*DrawResiParams[0];
  Float_t const1 = fun3->GetParameter(0);
  Float_t const2 = fun3->GetParameter(3);
  Float_t ratio  = 0.;
  if( sig1 < sig2 ) {
    sigpt->AddText(Form("#sigma_{1}=%.1lf #mum",sig1));
    sigpt->AddText(Form("#sigma_{2}=%.1lf #mum",sig2));
    sigpt->AddText(Form("#mu_{1}=%.1lf #mum",mean1));
    sigpt->AddText(Form("#mu_{2}=%.1lf #mum",mean2));
    if( sig2 > 0. && const2 > 0. ) {
      ratio = sig1*const1/sig2/const2;
      sigpt->AddText(Form("A_{1}/A_{2}=%.3lf",ratio));
    } 
    ResiFitData[3] = mean1;
    ResiFitData[4] = sig1;
    ResiFitData[5] = mean2;
    ResiFitData[6] = sig2;
    ResiFitData[7] = ratio;
    ResiFitData[9] = fun3->GetParError(1)*DrawResiParams[0];
    ResiFitData[10] = fun3->GetParError(4)*DrawResiParams[0];
    ResiFitData[11] = fun3->GetParError(2)*DrawResiParams[0];
    ResiFitData[12] = fun3->GetParError(5)*DrawResiParams[0];
 } else {
    sigpt->AddText(Form("#sigma_{1}=%.1lf #mum",sig2));
    sigpt->AddText(Form("#sigma_{2}=%.1lf #mum",sig1));
    sigpt->AddText(Form("#mu_{1}=%.1lf #mum",mean2));
    sigpt->AddText(Form("#mu_{2}=%.1lf #mum",mean1));
    if( sig1 > 0. && const1 > 0. ) {
      ratio = sig2*const2/sig1/const1;
      sigpt->AddText(Form("A_{1}/A_{2}=%.3lf",ratio));
    }
    ResiFitData[3] = mean2;
    ResiFitData[4] = sig2;
    ResiFitData[5] = mean1;
    ResiFitData[6] = sig1;
    ResiFitData[7] = ratio;
    ResiFitData[9] = fun3->GetParError(4)*DrawResiParams[0];
    ResiFitData[10] = fun3->GetParError(1)*DrawResiParams[0];
    ResiFitData[11] = fun3->GetParError(5)*DrawResiParams[0];
    ResiFitData[12] = fun3->GetParError(2)*DrawResiParams[0];
  }
  ResiFitData[0] = hp->GetEntries();
  ResiFitData[1] = hp->GetMean()*DrawResiParams[0];
  ResiFitData[2] = hp->GetRMS()*DrawResiParams[0];
  if( ResiFitData[0] > 0. ) ResiFitData[8] = ResiFitData[2]/TMath::Sqrt(ResiFitData[0]);
  ResiFitData[13] = Xmax;
  ResiFitData[14] = sigmaHM;
// Xmax, sigmaHM errors are errors of G1+G2 added in quadrature 
// These are too big.  Just use narrow gaussian errors
//   ResiFitData[15] = TMath::Sqrt(ResiFitData[9]*ResiFitData[9]+ResiFitData[10]*ResiFitData[10]);
//   ResiFitData[16] = TMath::Sqrt(ResiFitData[11]*ResiFitData[11]+ResiFitData[12]*ResiFitData[12]);
  ResiFitData[15] = ResiFitData[9];
  ResiFitData[16] = ResiFitData[11];
 
  sigpt->AddText(Form("#mu_{#epsilon}=%.1lf #mum",Xmax));
  sigpt->AddText(Form("#sigma_{#epsilon}=%.1lf #mum",sigmaHM));
  sigpt->Draw();

  makeplot(hp->GetName());

  return 0;
}   //end DrawResiInd()

Double_t DoubleG(Double_t *x,Double_t *par) {
  return G1(x,par)+G2(x,&par[2]);
} // DoubleG ================================================================

Double_t G1(Double_t *x,Double_t *par) {
  // 0=peak  1=sigma  2=mean
  Double_t t = (x[0]-par[2])/par[1];
  return par[0]*TMath::Exp(-0.5*t*t);
} // G1 ================================================================

Double_t G2(Double_t *x,Double_t *par) {
  // 0=mean  1=peak  2=sigma  
  Double_t t = (x[0]-par[0])/par[2];
  return par[1]*TMath::Exp(-0.5*t*t);
} // G2 ================================================================

//  Version to do fit 2 plots at once to display ML data
int DrawResiML(const char *hname1, const char *hname2, const char *pname, const char *fname ) {
//  By using TH1 this works for TH1F, TH1D, etc
  TH1 *hp1 = (TH1 *) gDirectory->Get(hname1);
  if( hp1 == 0 ) {
    printf("DrawResiML ERROR: No histogram %s\n",hname1);
    return -1;
  }
  TH1 *hp2 = (TH1 *) gDirectory->Get(hname2);
  if( hp2 == 0 ) {
    printf("DrawResiML ERROR: No histogram %s\n",hname2);
    return -1;
  }
  return DrawResiML( hp1, hp2, pname, fname );
}

//  Version fiting 2 gaussians constrained to have the same mean
int DrawResiML(TH1 *hp1, TH1 *hp2, const char *pname, const char *fname ) {
// Zero ResiFitData array
  for( int i=0; i<FITDATAMAX; i++ ) {
    ResiFitData[i] = 0.;
    ResiFitData2[i] = 0.;
  }
//  Do not plot if no entries or RMS==0 (entries are all underflow/overflow)
  if( hp1->GetEntries() == 0 || hp1->GetRMS() == 0 ) {
    printf("DrawResi WARNING: histogram %s no entries, skipping plot\n",hp1->GetName());
    return -1;
  }

  printf("DrawResi double fit: Fitting %s %s\n",hp1->GetName(),hp2->GetName());

// Add filename to title of histogram
// Rename histogram
  if( plottitle[0] != '\0' ) {
    if( pname[0] != '\0' ) {
      hp1->SetTitle(Form("%s %s",plottitle,pname));
    } else {
      hp1->SetTitle(Form("%s %s",plottitle,hp1->GetTitle()));
    }
  } else if( pname[0] != '\0' ) {  
    hp1->SetTitle(pname);
  }

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  //  gStyle->SetTitleX(0.5);
  //  gStyle->SetTitleH(0.25);
  //  gStyle->SetTitleW(0.95);

  // Reset hp1 maximum if hp2 max is larger so hp2 does not go offscale
  if( hp2!=0 && hp1->GetMaximum() < hp2->GetMaximum() ) {
    hp1->SetMaximum( 1.05*hp2->GetMaximum() );
  }

  // Do the double gaussian fits
  FitDoubleGaussian(hp1,ResiFitData);
  FitDoubleGaussian(hp2,ResiFitData2);

  TF1 *fun = hp1->GetFunction("fun3");
  fun->SetLineColor(kBlue);
  hp1->SetLineColor(kBlue);
  hp1->Draw();

  TPaveText *sigpt1 = new TPaveText(0.1,0.53,0.40,0.9,"NDC");
  sigpt1->SetTextAlign(12);
  sigpt1->SetTextColor(kBlue);
  sigpt1->AddText(Form("ML1 %.0lf",hp1->GetEntries()));
  sigpt1->AddText(Form("#sigma_{1}=%.1lf #mum",ResiFitData[4]));
  if( ResiFitData[6] > 1000. ) {
    sigpt1->AddText(Form("#sigma_{2}=%.1lf mm",ResiFitData[6]/1000.));
  } else {
    sigpt1->AddText(Form("#sigma_{2}=%.1lf #mum",ResiFitData[6]));
  }
  sigpt1->AddText(Form("#mu=%.1lf #mum",ResiFitData[3]));
  sigpt1->Draw();

  if( hp2 && hp2->GetEntries() > 0. ) {
    hp2->SetLineColor(kRed);
    hp2->Draw("same");
    TPaveText *sigpt2 = new TPaveText(0.64,0.53,0.94,0.9,"NDC");
    sigpt2->SetTextAlign(12);
    sigpt2->SetTextColor(kRed);
    sigpt2->AddText(Form("ML2 %.0lf",hp2->GetEntries()));
    sigpt2->AddText(Form("#sigma_{1}=%.1lf #mum",ResiFitData2[4]));
    if( ResiFitData2[6] > 1000. ) {
      sigpt1->AddText(Form("#sigma_{2}=%.1lf mm",ResiFitData2[6]/1000.));
    } else {
      sigpt2->AddText(Form("#sigma_{2}=%.1lf #mum",ResiFitData2[6]));
    }
    sigpt2->AddText(Form("#mu=%.1lf #mum",ResiFitData2[3]));
    sigpt2->Draw();
  }

// Print plot to file
  if( filename[0] != '\0' ) {
    if( fname[0] != '\0' ) {
      gPad->Print(Form("%s_%s.png",filename,fname));
    } else {
      gPad->Print(Form("%s_%s.png",filename,hp1->GetTitle()));
    }
  } else if( fname[0] != '\0' ) {  
    gPad->Print(Form("%s.png",hp1->GetTitle()));
  }

  //  if( fname[0] == '\0' ) gPad->Print(Form("%s.png",hp1->GetName()));
  //else                   gPad->Print(Form("%s.png",fname));

  return 0;
} //end DrawResiML


// Do a DoubleGaussian fit on a histogram.  Mean of 2 gaussians constrained to be the same
int FitDoubleGaussian(TH1 *hp, Double_t values[FITDATAMAX]) {

  if( hp == 0 ) return -1;

//  Set initial guesses of parameters
  Int_t    imax = hp->GetMaximumBin();
  Double_t mean = hp->GetMean();
  Float_t  ampl = (hp->GetBinContent(imax-1)+hp->GetBinContent(imax)+hp->GetBinContent(imax+1))/3.;
  Float_t  rms  = hp->GetRMS();

// Set fit limits according to param[1]
  Double_t xmin=hp->GetXaxis()->GetBinLowEdge(1);
  Double_t xmax=hp->GetXaxis()->GetBinUpEdge( hp->GetXaxis()->GetNbins() );
  if(      DrawResiParams[1] > 0. ){
    xmin=mean-DrawResiParams[1]*rms;
    xmax=mean+DrawResiParams[1]*rms;
  } else if( DrawResiParams[1] < 0. ){
    xmin= DrawResiParams[1];
    xmax=-DrawResiParams[1];
  }

  TF1 *fun3 = new TF1("fun3",DoubleG,xmin,xmax,5);
  fun3->SetParNames("Const1","Sigma1","Mean","Const2","Sigma2");
  fun3->SetParameters(ampl,rms/4,mean,ampl/10.,rms*2); // constant,mean,sigma

  hp->Fit("fun3","RQ");

//  Determine SigmaHM.  Do binary search to find half max points on both sides
  Double_t Xmax = fun3->GetMaximumX(mean-rms,mean+rms);
  Double_t Hmax = fun3->Eval(Xmax)/2;   //GetMaximum(mean-rms,mean+rms);
  Double_t X1 = Xmax-2*rms;
  Double_t X2 = Xmax;
  Double_t Xlo = (X1+X2)/2.;
  for( int ii=0; ii<25; ii++ ) {
    if( fun3->Eval(Xlo) < Hmax ) X1 = Xlo;
    else                         X2 = Xlo;
    Xlo = (X1+X2)/2.;
  }
  X1 = Xmax;
  X2 = Xmax+2*rms;
  Double_t Xhi = (X1+X2)/2.;
  for( int ii=0; ii<25; ii++ ) {
    if( fun3->Eval(Xhi) < Hmax ) X2 = Xhi;
    else                         X1 = Xhi;
    Xhi = (X1+X2)/2.;
  }
  Xmax *= DrawResiParams[0];
  Double_t sigmaHM = (Xhi-Xlo)/2.3548*DrawResiParams[0];
  Float_t sig1   = TMath::Abs(fun3->GetParameter(1))*DrawResiParams[0];
  Float_t sig2   = TMath::Abs(fun3->GetParameter(4))*DrawResiParams[0];
  Float_t mean1  = fun3->GetParameter(2)*DrawResiParams[0];
  Float_t mean2  = fun3->GetParameter(2)*DrawResiParams[0];
  Float_t const1 = fun3->GetParameter(0);
  Float_t const2 = fun3->GetParameter(3);
  if( sig1 < sig2 ) {
    values[3] = mean1;
    values[4] = sig1;
    values[5] = mean2;
    values[6] = sig2;
    values[7] = sig1*const1/sig2/const2;
    values[9]  = fun3->GetParError(2)*DrawResiParams[0];
    values[10] = fun3->GetParError(2)*DrawResiParams[0];
    values[11] = fun3->GetParError(1)*DrawResiParams[0];
    values[12] = fun3->GetParError(4)*DrawResiParams[0];
 } else {
    values[3] = mean2;
    values[4] = sig2;
    values[5] = mean1;
    values[6] = sig1;
    values[7] = sig2*const2/sig1/const1;
    values[9]  = fun3->GetParError(2)*DrawResiParams[0];
    values[10] = fun3->GetParError(2)*DrawResiParams[0];
    values[11] = fun3->GetParError(4)*DrawResiParams[0];
    values[12] = fun3->GetParError(1)*DrawResiParams[0];
  }
  values[0] = hp->GetEntries();
  values[1] = hp->GetMean()*DrawResiParams[0];
  values[2] = hp->GetRMS()*DrawResiParams[0];
  if( values[0] > 0. ) values[8] = values[2]/TMath::Sqrt(values[0]);
  values[13] = Xmax;
  values[14] = sigmaHM;
  values[15] = values[9];
  values[16] = values[11];

  return 0;
} //end FitDoubleGaussian


/************************************************************
*  Draw individual 1D histograms with a single Gaussian fit
*  ResiFitData[0] = Number of entries
*  ResiFitData[1] = Histogram mean
*  ResiFitData[2] = RMS
*  ResiFitData[3] = Gasssian mean;
*  ResiFitData[4] = Gaussian sigma;
*  ResiFitData[5] = Histogram mean error = rms/sqrt(n)
*  ResiFitData[6] = Gaussian mean error
*  ResiFitData[7] = Gaussian sigma error
************************************************************/
TH1 * DrawGaus(const char *hname ) {
//  By using TH1 this works for TH1F, TH1D, etc
  TH1 *hp = (TH1 *) gDirectory->Get(hname);
  if( hp == NULL ) {
    printf("DrawGaus ERROR: No histogram %s\n",hname);
    return NULL;
  }
  DrawGaus( hp );
  return hp;
}
int DrawGaus( TH1 *hp ) {
// Zero ResiFitData array
  for( int i=0; i<FITDATAMAX; i++ ) ResiFitData[i] = 0.;

//  By using TH1 this works for TH1F, TH1D, etc
  if( hp == NULL ) {
    printf("DrawGaus ERROR: Null histogram pointer passed\n");
    return -1;
  }
//  Do not plot if no entries or if RMS=0
  if( hp->GetEntries() == 0 || hp->GetRMS() == 0 ) {
    printf("DrawGaus WARNING: histogram %s no entries, skipping plot\n",hp->GetName());
    return -1;
  }

  printf("DrawGaus: Fitting %s\n",hp->GetName());

// Add filename to title of histogram
  if( plottitle[0] != '\0' ) {
    hp->SetTitle(Form("%s %s",plottitle,hp->GetTitle()));
  }

  gStyle->SetOptStat("emr"); 
//   gStyle->SetStatX(1.0);
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.35);

  TF1 *fun1=new TF1("fun1","gaus");
  fun1->SetParNames("Const","Mean","Sigma");
//  Set initial guesses of parameters
  Int_t       i = hp->GetMaximumBin();
  Double_t ampl = (hp->GetBinContent(i-1) + hp->GetBinContent(i) + hp->GetBinContent(i+1))/3.;
  Double_t mean = hp->GetMean();
  Double_t rms  = hp->GetRMS();
  Double_t peak = hp->GetBinCenter(hp->GetMaximumBin());
//  Find narrow peak on a background.  
  Double_t  lo, hi;
//   Double_t halfmax = hp->GetMaximum()/3.;     //set to half max
//   lo = hp->GetBinCenter(1);
//   hi = hp->GetBinCenter(hp->GetNbinsX()-1);
//  Find lo
//   for( int i=hp->GetMaximumBin(); i>0; i-- ) {
//     if( hp->GetBinContent(i) < halfmax ) {
//       lo = hp->GetBinCenter(i);
//       break;
//     } 
//   }
//  Find hi
//   for( int i=hp->GetMaximumBin(); i<hp->GetNbinsX(); i++ ) {
//     if( hp->GetBinContent(i) < halfmax ) {
//       hi = hp->GetBinCenter(i);
//       break;
//     } 
//   }
/*   lo = mean-rms; */
/*   hi = mean+rms;   */
  lo = peak-6;
  hi = peak+6;  
  fun1->SetParameters(ampl,peak,rms);    //constant,mean,sigma
  hp->Fit("fun1","","",lo,hi);
//  Draw fit results on plot
//  TPaveText *sigpt = new TPaveText(0.16,0.79,0.45,0.99,"NDC");  //without histo title 
//  TPaveText *sigpt = new TPaveText(0.16,0.74,0.45,0.94,"NDC");    //with histo title
  TPaveText *sigpt = new TPaveText(0.1,0.7,0.4,0.9,"NDC"); 
  sigpt->SetFillColor(0);
  sigpt->SetTextAlign(12);
  Float_t mean1 = fun1->GetParameter(1);
  Float_t sig1  = TMath::Abs(fun1->GetParameter(2));
  //  sigpt->AddText(Form("Peak=%.1lf ns",peak));
  sigpt->AddText(Form("#mu=%.1lf ns",mean1));
  sigpt->AddText(Form("#sigma=%.1lf ns",sig1));
  sigpt->Draw();

  makeplot(hp->GetName());

//  We already know entries>0
  ResiFitData[0] = hp->GetEntries();
  ResiFitData[1] = mean;
  ResiFitData[2] = rms;
  if( ResiFitData[0] > 0. ) ResiFitData[5] = ResiFitData[2]/TMath::Sqrt(ResiFitData[0]);
// Make sure gaussian mean falls within limits of the histogram range
// i.e. check that fit is reasonable.
  if( mean < hp->GetXaxis()->GetXmax() && mean > hp->GetXaxis()->GetXmin() ) {
    ResiFitData[3] = mean1;
    ResiFitData[4] = sig1;
    ResiFitData[6] = fun1->GetParError(1);
    ResiFitData[7] = fun1->GetParError(2);
  }
  
  return 0;
}  //end DrawGaus()

/************************************************************
*  Draw individual 1D histograms with a single Gaussian fit
*  ResiFitData[0] = Number of entries
*  ResiFitData[1] = Histogram mean
*  ResiFitData[2] = RMS
*  ResiFitData[3] = Gasssian mean;
*  ResiFitData[4] = Gaussian sigma;
*  ResiFitData[5] = Histogram mean error = rms/sqrt(n)
*  ResiFitData[6] = Gaussian mean error
*  ResiFitData[7] = Gaussian sigma error
************************************************************/
TH1 * DrawT0refit(const char *hname ) {
  ResiFitData[0] = 0;
//  By using TH1 this works for TH1F, TH1D, etc
  TH1 *hp = (TH1 *) gDirectory->Get(hname);
  if( hp == NULL ) {
    printf("DrawT0refit ERROR: No histogram %s\n",hname);
    return NULL;
  }
  DrawT0refit( hp );
  return hp;
}

int DrawT0refit( TH1 *hp ) {

// Zero ResiFitData array
  for( int i=0; i<12; i++ ) ResiFitData[i] = 0;

//  By using TH1 this works for TH1F, TH1D, etc
  if( hp == NULL ) {
    printf("DrawT0refit ERROR: Null histogram pointer passed\n");
    return -1;
  }
//  Do not plot if no entries or if RMS=0
  if( hp->GetEntries() == 0 || hp->GetRMS() == 0 ) {
    printf("DrawT0refit WARNING: histogram %s no entries, skipping plot\n",hp->GetName());
    return -1;
  }

  printf("DrawT0refit: Fitting %s\n",hp->GetName());

// Add filename to title of histogram
  if( plottitle[0] != '\0' ) {
    hp->SetTitle(Form("%s %s",plottitle,hp->GetTitle()));
  }

  gStyle->SetOptStat("emr"); 
//   gStyle->SetStatX(1.0);
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.35);

  TF1 *gaus1=new TF1("gaus1","gaus");
  gaus1->SetParNames("Const","Mean","Sigma");
//  Set initial guesses of parameters
  Int_t       i = hp->GetMaximumBin();
  Double_t ampl = (hp->GetBinContent(i-1) + hp->GetBinContent(i) + hp->GetBinContent(i+1))/3.;
  Double_t mean = hp->GetMean();
  Double_t rms  = hp->GetRMS();
  Double_t peak = hp->GetBinCenter(hp->GetMaximumBin());
//  Find narrow peak on a background.  
  Double_t  lo, hi;
  lo = peak-6;
  hi = peak+6;  
  gaus1->SetParameters(ampl,peak,rms);    //constant,mean,sigma
  hp->Fit("gaus1","","",lo,hi);
//  Draw fit results on plot
//  TPaveText *sigpt = new TPaveText(0.16,0.79,0.45,0.99,"NDC");    //without histo title 
//  TPaveText *sigpt = new TPaveText(0.16,0.74,0.45,0.94,"NDC");    //with histo title
  TPaveText *sigpt = new TPaveText(0.1,0.7,0.4,0.9,"NDC"); 
  sigpt->SetFillColor(0);
  sigpt->SetTextAlign(12);
  sigpt->SetTextColor(2);
  Float_t mean1 = gaus1->GetParameter(1);
  Float_t sig1  = TMath::Abs(gaus1->GetParameter(2));
  //  sigpt->AddText(Form("Peak=%.1lf ns",peak));
  sigpt->AddText(Form("#mu=%.1lf ns",mean1));
  sigpt->AddText(Form("#sigma=%.1lf ns",sig1));
  sigpt->Draw();

//  We already know entries>0
  ResiFitData[0] = hp->GetEntries();
  ResiFitData[1] = mean;
  ResiFitData[2] = rms;
// Make sure gaussian mean falls within limits of the histogram range
// i.e. check that fit is reasonable.
//  if( mean < hp->GetXaxis()->GetXmax() && mean > hp->GetXaxis()->GetXmin() ) {
  ResiFitData[3] = mean1;
  ResiFitData[4] = sig1;
//  }
  
//  Look for a second peak 25ns below first peak
  TF1 *gaus2=new TF1("gaus2","gaus");
  gaus2->SetParNames("Const","Mean","Sigma");
  ampl /= 2.;
  peak = mean1 - 25.;
  rms  = sig1;
  lo = peak-12;
  hi = peak+12;  
  gaus2->SetParameters(ampl,peak,rms);    //constant,mean,sigma
  gaus2->SetParLimits(1,-100.,peak+12.);

  gaus2->SetLineColor(4);
  hp->Fit("gaus2","+","",lo,hi);

  printf("mean1=%.2f peak=%0.2f limit=%.2f\n",mean1,peak,peak+12.);
  printf("LO FIT %.2f %.2f %.2f %.2f %.2f\n",gaus1->GetParameter(0), gaus2->GetParameter(0), gaus2->GetParameter(1), gaus2->GetParameter(2), gaus2->GetParameter(1)-gaus1->GetParameter(1));

  TPaveText *sigpt2 = new TPaveText(0.1,0.5,0.4,0.7,"NDC"); 
  sigpt2->SetFillColor(0);
  sigpt2->SetTextAlign(12);
  sigpt2->SetTextColor(4);

  Double_t offset = TMath::Abs(gaus2->GetParameter(1)-gaus1->GetParameter(1));
  // if 2nd peak is not significant look for peak above main peak
  if( gaus2->GetParameter(0) < 0.05*gaus1->GetParameter(0) || offset < 15. ) {
    delete hp->GetListOfFunctions()->FindObject("gaus2"); 
    peak = mean1 + 25.;
    lo = peak-6;
    hi = peak+6;  
    gaus2->SetParameters(ampl,peak,rms);    //constant,mean,sigma
    gaus2->SetParLimits(1,peak-12.,100.);
    offset = TMath::Abs(gaus2->GetParameter(1)-gaus1->GetParameter(1));
    hp->Fit("gaus2","+","",lo,hi);
    printf("HI FIT %.2f %.2f %.2f %.2f\n",gaus1->GetParameter(0), gaus2->GetParameter(0), gaus2->GetParameter(1), gaus2->GetParameter(2));
  }

  // If 2nd peak is significant write results on plot, if not delete function from histogram.
  if(  gaus2->GetParameter(0) > 0.05*gaus1->GetParameter(0) && offset > 15. ) {
    Float_t mean2 = gaus2->GetParameter(1);
    Float_t sig2  = TMath::Abs(gaus2->GetParameter(2));
    sigpt2->AddText(Form("#mu_{2}=%.1lf ns",mean2));
    sigpt2->AddText(Form("#sigma_{2}=%.1lf ns",sig2));
    sigpt2->Draw();
    ResiFitData[5]  = mean2;
    ResiFitData[6]  = sig2;
    ResiFitData[7]  = gaus2->GetParameter(1)-gaus1->GetParameter(1);
  } else {
    delete hp->GetListOfFunctions()->FindObject("gaus2"); 
  } 

  makeplot(hp->GetName());
  return 0;
}  //end DrawT0refit()

/************************************************************
*  Draw individual 1D histograms
*  Histogram hname is read from the open ROOT file
************************************************************/
void Drawh1( const char *hname, const char *opt, const char *stat) {
//  By using TH1 this works for TH1F, TH1D, etc
  TH1 *hp = (TH1 *) gDirectory->Get(hname);
  if( hp == 0 ) {
    printf("Drawh1 ERROR: No histogram %s \n",hname);
    return;
  }

//  Do not plot if no entries or RMS==0 (entries are all underflow/overflow)
  if( hp->GetEntries() == 0 ) {
    printf("Drawh1 WARNING: histogram %s no entries, skipping plot\n",hname);
    return;
  }

// Rename histogram
  if( plottitle[0] != '\0' ) {
    hp->SetTitle(Form("%s %s",plottitle,hp->GetTitle()));
  }

  gStyle->SetOptStat(stat); 
//   gStyle->SetStatX(1.0);
//   gStyle->SetStatW(0.25);
//   gStyle->SetStatY(1.0);
//   gStyle->SetStatH(0.2);

//  Reset various parameters of the histogram
  hp->Draw(opt);
  makeplot(hname);
  return;
}  //end Drawh1

/************************************************************
*  Draw individual 2D histograms
************************************************************/
TH2F* Drawh2(const char *hname, const char *chopt, const char *xlab, const char *ylab, const char *chstat ) {
//  By using TH2 this works for TH2F, TH2D, etc
  TH2F *hp = (TH2F *) gDirectory->Get(hname);
  if( hp == 0 ) {
    printf("Drawh2 ERROR: No histogram %s\n",hname);
    return 0;
  }
  Drawh2( hp, chopt, xlab, ylab, chstat );
  return hp;
}

/************************************************************
*  Draw individual 2D histograms
************************************************************/
void Drawh2( TH2F *hp, const char *chopt, const char *xlab, const char *ylab, const char *chstat ) {
// Skip non-existant & empty histograms
  if( hp == 0 || hp->GetEntries() == 0. ) return;

// Rename histogram
  if( plottitle[0] != '\0' ) {
    hp->SetTitle(Form("%s %s",plottitle,hp->GetTitle()));
  }

  gStyle->SetOptStat(chstat); 
//   gStyle->SetStatX(1.0);
//   gStyle->SetStatW(0.25);
//   gStyle->SetStatY(1.0);
//   gStyle->SetStatH(0.2);

//  Reset various parameters of the histogram
  if( xlab[0] != '\0' ) hp->GetXaxis()->SetTitle(xlab); 
  if( ylab[0] != '\0' ) hp->GetYaxis()->SetTitle(ylab); 
  hp->Draw(chopt);
  makeplot(hp->GetName());
  return;
}

// void TwoPlot( TH1 *h[], int inorm, const char *title ) {
//   printf("TwoPlot inorm = %i \n",inorm);
//   for( int i=0; i<2; i++ ) {
//     if( h[i] ) {
//       printf("TwoPlot segfault now histo %s\n",h[0]->GetName());
//     } else {
//       printf("h[%i] is null\n",i);
//     }
//   }
//   gStyle->SetOptStat(0);
//   gPad->SetTicks(0,0);

//   if( h[0] == NULL || h[1] == NULL ) {
//     printf("ERROR: Twoplot null pointer\n");
//     return;
//   }
//   h[0]->SetFillColor(3);
//   h[0]->Draw();

//   gPad->Update();   //I don't know what this does but if you don't do it the axis will not work
//   Float_t rightmax = 1.1*h[0]->GetMaximum();
//   if( inorm == 1 && h[1]->GetEntries()>0) {
//     printf("inorm = 1\n");
//     rightmax *= h[0]->GetEntries()/h[1]->GetEntries();
//     h[1]->Scale(h[0]->GetEntries()/h[1]->GetEntries());
//   } else if (inorm == 2 && h[1]->GetMaximum()>0) {
//     printf("inorm = 2\n");
//     rightmax = 1.1*h[1]->GetMaximum();
//     h[1]->Scale(gPad->GetUymax()/rightmax);
//   }
// //draw an axis on the right side
//   h[1]->Draw("SAME");
// //  Using NDC coordinates, as per ROOT manual example, but which does not work
//   TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
// 			    gPad->GetUymax(),0,rightmax,510,"NDC+L");
//   axis->Draw();

// //  Make a title from the titles of the 2 histograms else use the one supplied
//   h[0]->SetBit(TH1::kNoTitle);   //suppress the original histogram title
// //  TPaveText *titlept = new TPaveText(0.02,0.92,0.75,0.99,"NDC");   //with histo title
//   TPaveText *titlept = new TPaveText(0.02,0.82,0.75,0.89,"NDC");       //without histo title
//   if( title[0] != '\0' ) {
//     titlept->AddText(title);                        
//   } else {
//     titlept->AddText(Form("%s (L) & %s (R)",h[0]->GetTitle(),h[1]->GetTitle()));
//   }
//   titlept->Draw();     

// //  Print to a file
//   gPad->Print(Form("%s_%s_TP.png",h[0]->GetName(),h[1]->GetName()));
//   h[0]->ResetBit(TH1::kNoTitle);   //enable the histogram title
// }     //End TwoPlot
/*########################################################
  Draw 2 histograms to TCanvas c1 and save to a file.
   hname1, hname2 = Names of histograms to be plotted (assumed to be in gDirectory)
   iflag = Indicates type of pair of histograms being plotted (which controls where ResiFitData is stored)
         =0  for plotting RAW & Segment plots overlaid
             ResiFitData[0,1,2] = events, mean, rms of hname1, 
             ResiFitData[3,4,5] = events, mean, rms of hname2
         =1  for plotting RAW ML plots
             ResiFitData[0,1,2]  = events, mean, rms of hname1, 
             ResiFitData2[0,1,2] = events, mean, rms of hname2
         =2  for plotting Segment ML plots
             ResiFitData[3,4,5]  = events, mean, rms of hname1, 
             ResiFitData2[3,4,5] = events, mean, rms of hname2
   pname = Plot title (overwrites hname1 title)
   fname = File name for image file, including extension (png, gif, etc)
       Default filename is <Histogram Name>.png
   ResiFitData = Array to write histogram stats data 
  Stats box for second histogram is drawn below first histo stats box
  ########################################################*/
void Drawh1h1(const char *hname1, const char *hname2, const int iflag, const char *pname, const char *fname ) {
  for( int i=0; i<6; i++ ) {
    if( iflag == 0 ) ResiFitData[i]=0.;
    if( (iflag==1 && i<3) || (iflag==2 && i>2) ) {
      ResiFitData[i]=0.;  
      ResiFitData2[i]=0.;  
    }
  }

//  By using TH1 this works for TH1F, TH1D, etc
  TH1 *h1 = (TH1 *) gDirectory->Get(hname1);
  if( h1 == 0 ) {
    printf("Drawh1h1 ERROR: No histogram %s\n",hname1);
    return;
  }
// Get the other histogram
  TH1 *h2 = (TH1 *) gDirectory->Get(hname2);

//  Do not plot if no entries 
  if( h1->GetEntries()==0. && h2!=0 && h2->GetEntries()==0. ) {
    printf("Drawh1h1 WARNING: histograms %s and/or %s no entries, skipping plot\n",hname1,hname2);
    return;
  }

// Name for PNG file (do before changing histogram title)
  char pngname[80] = "file.png";
  if( fname[0] != '\0' ) {
    if( filename[0] != '\0' ) {
      sprintf(pngname,"%s_%s.png",filename,fname);
    } else {
      sprintf(pngname,"%s.png",fname);
    }
  } else if( filename[0] != '\0' ) {  
    sprintf(pngname,"%s_%s.png",filename,h1->GetName());
  } else {
    sprintf(pngname,"%s.png",h1->GetName());
  }

// Rename histogram
  if( plottitle[0] != '\0' ) {
    if( pname[0] != '\0' ) {
      h1->SetTitle(Form("%s %s",plottitle,pname));
    } else {
      h1->SetTitle(Form("%s %s",plottitle,h1->GetTitle()));
    }
  } else if( pname[0] != '\0' ) {  
    h1->SetTitle(pname);
  }

  // Reset h1 maximum if h2 max is larger so h2 does not go offscale
  if( h2!=0 && h1->GetMaximum() < h2->GetMaximum() ) {
    h1->SetMaximum( 1.05*h2->GetMaximum() );
  }

  // Set h1 min to 0 to prevent h2 from being cut off if it has many fewer events 
  h1->SetMinimum(0.);
  h1->Draw();

  if( iflag==0 || iflag==1 ) {
    ResiFitData[0] = h1->GetEntries();
    ResiFitData[1] = h1->GetMean();
    ResiFitData[2] = h1->GetRMS();
  } else if( iflag==2 ) {
    ResiFitData[3] = h1->GetEntries();
    ResiFitData[4] = h1->GetMean();
    ResiFitData[5] = h1->GetRMS();
  }

  if( h2!=0 && h2->GetEntries()!=0. ) {
    if( iflag==0 ) {
      ResiFitData[3] = h2->GetEntries();
      ResiFitData[4] = h2->GetMean();
      ResiFitData[5] = h2->GetRMS();
    } else if( iflag==1 ) {
      ResiFitData2[0] = h2->GetEntries();
      ResiFitData2[1] = h2->GetMean();
      ResiFitData2[2] = h2->GetRMS();
    } else if( iflag==2 ) {
      ResiFitData2[3] = h2->GetEntries();
      ResiFitData2[4] = h2->GetMean();
      ResiFitData2[5] = h2->GetRMS();
    }
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);
    h2->SetLineColor(kRed);
    h2->Draw("sames");
    gPad->Update();
//  Change color of h2 stats box and move it below the other stats box
    TPaveStats *s = (TPaveStats*) h2->FindObject("stats");
    if( s == 0 ) {
      printf("Drawh1h1 ERROR finding stats box %s\n",hname2);
    } else {
      float ymin = s->GetY1NDC();
      float ymax = s->GetY2NDC();
      s->SetY2NDC(ymin);
      s->SetY1NDC(2*ymin-ymax);
      s->SetTextColor(kRed);
    }
  }

// Print plot to file
  gPad->Print(pngname);
  return;
}  //Drawh1h1

void Drawh1h1h1(const char *hname1, const char *hname2, const char *hname3, const char *pname ) {
//  By using TH1 this works for TH1F, TH1D, etc
  TH1 *h1 = (TH1 *) gDirectory->Get(hname1);
  if( h1 == 0 ) {
    printf("Drawh1h1h1 ERROR: No histogram %s\n",hname1);
    return;
  }

  //  Do not plot if no entries 
  if( h1->GetEntries() == 0. ) {
    printf("Drawh1h1h1 WARNING: histogram %s no entries, skipping plot\n",hname1);
    return;
  }

  TH1 *h2 = (TH1 *) gDirectory->Get(hname2);
  if( h2==0 ) {
    printf("Drawh1h1h1 WARNING: NO histogram %s\n",hname2);
  }

  Drawh1h1h1( h1, h2, hname3, pname );
  return;
}   //end Drawh1h1h1

/*########################################################
  Draw 3 histograms and save to a file.
   fname = File name for image file, including extension (png, gif, etc)
       Default filename is <Histogram Name>.png
  Second histogram is drawn on first histogram in different color set by fillcolor
  Stats box for second histogram is drawn below first histo stats box
  ########################################################*/
void Drawh1h1h1(TH1 *h1, TH1 *h2, const char *hname3, const char *pname ) {
  for( int i=0; i<6; i++ ) ResiFitData[i]=0.;

//  By using TH1 this works for TH1F, TH1D, etc
  if( h1 == 0 ) return;

//  Do not plot if no entries 
  if( h1->GetEntries() == 0. ) {
    printf("Drawh1h1h1 WARNING: histogram %s no entries, skipping plot\n",h1->GetName());
    return;
  }

// Rename histogram
  if( plottitle[0] != '\0' ) {
    if( pname[0] != '\0' ) {
      h1->SetTitle(Form("%s %s",plottitle,pname));
    } else {
      h1->SetTitle(Form("%s %s",plottitle,h1->GetTitle()));
    }
  }

  h1->SetMinimum(0.);
  h1->SetLineWidth(2);
  h1->Draw();
  gPad->Update();

  ResiFitData[0] = h1->GetEntries();
  ResiFitData[1] = h1->GetMean();
  ResiFitData[2] = h1->GetRMS();

  if( h2!=0 && h2->GetEntries()!=0. ) {
    h2->SetLineColor(kGreen);
    h2->SetLineWidth(2);
    h2->Draw("sames");
    gPad->Update();
//  Change color of h2 stats box and move it below the other stats box
    TPaveStats *s = (TPaveStats*) h2->FindObject("stats");
    if( s == 0 ) {
      printf("Drawh1h1h1 ERROR finding stats box %s\n",h2->GetName());
    } else {
      float ymin = s->GetY1NDC();
      float ymax = s->GetY2NDC();
      s->SetY2NDC(ymin);
      s->SetY1NDC(2*ymin-ymax);
      s->SetTextColor(kGreen);
    }
  } 

  TH1 *h3 = (TH1 *) gDirectory->Get(hname3);
  if( h3!=0 && h3->GetEntries()!=0. ) {
    ResiFitData[3] = h3->GetEntries();
    ResiFitData[4] = h3->GetMean();
    ResiFitData[5] = h3->GetRMS();
    h3->SetLineColor(kRed);

//scale h3 to the pad coordinates
    Float_t rightmax = 1.1*h3->GetMaximum();
    Float_t scale = gPad->GetUymax()/rightmax;
    h3->SetLineWidth(2);
    h3->Scale(scale);
    //   h3->Draw("same");
//draw an axis on the right side
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
			      gPad->GetUymax(),0,rightmax,510,"-L");
    axis->SetLineColor(kRed);
    axis->SetLabelColor(kRed);
    axis->Draw();

    h3->Draw("sames");
    gPad->Update();
//  Change color of h3 stats box and move it below the other stats box
    TPaveStats *s = (TPaveStats*) h3->FindObject("stats");
    if( s == 0 ) {
      printf("Drawh1h1h1 ERROR finding stats box %s\n",hname3);
    } else {
      float ymin = s->GetY1NDC();
      float ymax = s->GetY2NDC();
      s->SetY2NDC(2*ymin-ymax);
      s->SetY1NDC(3*ymin-2*ymax);
      s->SetTextColor(kRed);
    }
  } else {
    printf("Drawh1h1h1 WARNING: NO histogram %s\n",hname3);
  }

  makeplot(h1->GetName());
  return;
}

// utility program to make a plot file from whatever has just been plotted on the current canvas
void makeplot(const char *pname ) {
  if( filename[0] == '\0' ) gPad->Print(Form("%s.%s",pname,chend));
  else                      gPad->Print(Form("%s_%s.%s",filename,pname,chend));
}
