/*******************************************************************************
  file name: FitResidual.cxx
  author: Zhe Yang
  created: 03/15/2019
  last modified: 03/15/2019

  description:
  Plot double-gaussian fitting result of residual distribution

  remark:
  -search window is set to (0, 100ns)

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

#include "AtlasStyle.C"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLine.h"
#include "TTree.h"

#include "plotresi.h"

#define _USE_MATH_DEFINES

string s_input_filename = "run00187691_20190301.dat.dir/noise_level_output/noise_level_output031901_run00187691_20190301.dat.noise.root";

using namespace std;

int FitResidual() {
  //SetAtlasStyle();

  // read input rootfile and set branch
  TFile *p_input_rootfile = new TFile(s_input_filename.c_str());

  // prepare canvas
  char hist_name[256];
  char output_name[256];
  strcpy(hist_name, "residual_distribution");
  strcpy(output_name, hist_name);
  strcat(output_name, ".png");
  TCanvas *p_canvas_residual = new TCanvas(hist_name, hist_name);
  TCanvas *p_canvas_residual_vs_radius[8];
  TCanvas *p_canvas_sigma = new TCanvas("sigma_vs_radius", "sigma_vs_radius");
  for (Int_t radius_id = 0; radius_id != 8; radius_id++) {
    sprintf(hist_name, "residual_distribution_%dmm", radius_id);
    p_canvas_residual_vs_radius[radius_id] = new TCanvas(hist_name, hist_name);
  }

  // get histogram and do the fitting
  strcpy(hist_name, "residual_distribution");
  TH1F *p_residual_hist;
  TH1F *p_residual_hist_vs_radius[8];
  TH1F *p_sigma_vs_radius = new TH1F("sigma_vs_radius", "sigma_vs_radius", 8, 0, 8);
  
  Double_t sigma[8];
  p_input_rootfile->GetObject(hist_name, p_residual_hist);
  p_canvas_residual->cd();
  p_residual_hist->GetXaxis()->SetTitle("radius/mm");
  p_residual_hist->GetYaxis()->SetTitle("entries");
  p_residual_hist->Draw();
  DrawResi(p_residual_hist);
  p_canvas_residual->SaveAs(output_name);

  for (Int_t radius_id = 0; radius_id != 8; radius_id++) {
    sprintf(hist_name, "residual_distribution_%dmm", radius_id);
    strcpy(output_name, hist_name);
    strcat(output_name, ".png");
    p_input_rootfile->GetObject(hist_name, p_residual_hist_vs_radius[radius_id]);
    p_canvas_residual_vs_radius[radius_id]->cd();
    p_residual_hist_vs_radius[radius_id]->GetXaxis()->SetTitle("radius/mm");
    p_residual_hist_vs_radius[radius_id]->GetYaxis()->SetTitle("entries");
    p_residual_hist_vs_radius[radius_id]->Draw();
    sigma[radius_id] = DrawResi(p_residual_hist_vs_radius[radius_id]);
    p_canvas_residual_vs_radius[radius_id]->SaveAs(output_name);
    
    p_sigma_vs_radius->Fill(radius_id, sigma[radius_id]);
    cout << "sigma between " << radius_id << " mm and " << radius_id + 1 << " mm is: " << sigma[radius_id] << endl;
  }

  p_canvas_sigma->cd();
  p_sigma_vs_radius->GetYaxis()->SetRange(0, 300);
  p_sigma_vs_radius->SetStats(0);
  p_sigma_vs_radius->Draw("hist");
  sprintf(output_name, "sigma_vs_radius.png");
  p_canvas_sigma->SaveAs(output_name);

  return 0;
}