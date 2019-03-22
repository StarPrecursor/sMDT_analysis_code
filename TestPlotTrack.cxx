/*******************************************************************************
  file name: TestPlotTrack.cxx
  author: Zhe Yang
  created: 02/27/2019
  last modified: 02/27/2019

  description:
  -use toy MC data to test PlotTrack.cxx's validity

*******************************************************************************/

#include "PlotTrack.h"

#include <stdio.h>
#include <math.h>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TRandom.h"
#include "TTree.h"

#include "DecodeRawData.h"
#include "GetHitInfo.cxx"
#include "GetHitLayerColumn.cxx"
#include "rtfunction.h"

#define _USE_MATH_DEFINES

const unsigned long selected_event_id = 20;
//string s_input_filename = "run00187573_20190131.dat.dir/run00187573_20190131.dat.out.root";
//string s_input_filename = "run00187603_20190211.dat.dir/run00187603_20190211.dat.out.root";
//string s_input_filename = "run00187614_20190214.dat.dir/run00187614_20190214.dat.out.root";
//string s_input_filename = "run00187617_20190215.dat.dir/run00187617_20190215.dat.out.root";
//string s_input_filename = "run00187631_20190218.dat.dir/run00187631_20190218.dat.out.root";
string s_input_filename = "run00187664_20190222.dat.dir/run00187664_20190222.dat.out.root";
const Double_t FITTING_PARA0 = 25000;
const Double_t FITTING_PARA1 = -80;
const Double_t FITTING_PARA2 = 10;
const Double_t FITTING_PARA3 = -1;

using namespace std;

Double_t TestFunction(Double_t *x, Double_t *par) {
  return par[0] + par[1] * x[0];
}// used for debugging only

Double_t GetPointLineDistance(Double_t x, Double_t y, Double_t k, Double_t b);

int TestPlotTrack() {
  // create toy MC signal data
  Double_t mc_track_k = 3.59101, mc_track_b = -1394.14;
  Int_t event_signal_length = 8;
  Int_t event_signal[event_signal_length][5];

  event_signal[0][0] = 4;
  event_signal[0][1] = 9;
  event_signal[0][2] = 21;

  event_signal[1][0] = 4;
  event_signal[1][1] = 10;
  event_signal[1][2] = 0;

  event_signal[2][0] = 4;
  event_signal[2][1] = 11;
  event_signal[2][2] = 15;

  event_signal[3][0] = 4;
  event_signal[3][1] = 9;
  event_signal[3][2] = 23;

  event_signal[4][0] = 4;
  event_signal[4][1] = 11;
  event_signal[4][2] = 16;

  event_signal[5][0] = 4;
  event_signal[5][1] = 11;
  event_signal[5][2] = 13;

  event_signal[6][0] = 4;
  event_signal[6][1] = 9;
  event_signal[6][2] = 22;

  event_signal[7][0] = 4;
  event_signal[7][1] = 11;
  event_signal[7][2] = 18;

  // get toy MC signal drift distance
  Double_t drift_distance[100];
  TRandom *error_generator = new TRandom();
  for (Int_t signal_id = 0; signal_id < event_signal_length; signal_id++) {
    Double_t x, y, k, b;
    GetHitInfo(event_signal[signal_id][1], event_signal[signal_id][2], &x, 
               &y);
    drift_distance[signal_id] = GetPointLineDistance(x, y, mc_track_k, 
                                                     mc_track_b);
    drift_distance[signal_id] += error_generator->Gaus(0, 0.060);
  }

  // display event's data
  cout << "selected event's signal: " << endl
       << "header type | tdc id | channel id | drift distance" << endl;
  for (int signal_id = 0; signal_id < event_signal_length; signal_id++) {
    cout << event_signal[signal_id][0] << "\t" 
         << event_signal[signal_id][1] << "\t" 
         << event_signal[signal_id][2] << "\t" 
         << drift_distance[signal_id] << endl;
  }

  // prepare base of output track display
  const Double_t layer_distance = 13.0769836;
  const Double_t column_distance = 15.1;
  const Double_t radius = 7.5;

  TCanvas *track_base = new TCanvas("track base", "track base", 0, 0, 1200, 480);
  track_base->cd();
  Double_t center_x, center_y;
  Double_t track_corner_x[2] = {0, 800};
  Double_t track_corner_y[2] = {0, 320};
  TGraph * track_baseline = new TGraph(2, track_corner_x, track_corner_y);
  char track_name[256];
  sprintf(track_name, "selected_event_id_%lu", selected_event_id);
  track_baseline->SetNameTitle(track_name, track_name);
  track_baseline->Draw("AP");

  TEllipse *tube_model[54][8];
  TEllipse *hit_model[512];
  for (Int_t layer_id = 0; layer_id != 4; layer_id++) {
    for (Int_t column_id = 0; column_id != 54; column_id++) {
    center_x = 7.5 + column_id * column_distance + ((layer_id + 1) % 2) *
               column_distance / 2.0;
    center_y = 7.5 + layer_id * layer_distance;
    tube_model[layer_id][column_id] = new TEllipse(center_x, center_y,
                                                   radius, radius);
    if ((column_id / 6) % 2 == 0) {
      tube_model[layer_id][column_id]->SetFillColor(kGray);
    }
    tube_model[layer_id][column_id]->Draw();
    }
  }
  for (Int_t layer_id = 4; layer_id != 8; layer_id++) {
    for (Int_t column_id = 0; column_id != 54; column_id++) {
    center_x = 7.5 + column_id * column_distance + ((layer_id + 1) % 2) *
               column_distance / 2.0;
    center_y = 7.5 + (layer_id - 4) * layer_distance + 224.231;
    tube_model[layer_id][column_id] = new TEllipse(center_x, center_y,
                                                   radius, radius);
    if ((column_id / 6) % 2 == 0) {
      tube_model[layer_id][column_id]->SetFillColor(kGray);
    }
    tube_model[layer_id][column_id]->Draw();
    }
  }
  for (Int_t signal_id = 0; signal_id < event_signal_length; signal_id++) {
    if (event_signal[signal_id][0] == 4) {
      Double_t x_0, y_0;
      GetHitInfo(event_signal[signal_id][1], event_signal[signal_id][2], &x_0, 
                 &y_0);
      hit_model[signal_id] = new TEllipse(x_0, y_0, drift_distance[signal_id], 
                                          drift_distance[signal_id]);
      hit_model[signal_id]->SetLineColor(kBlue);
      hit_model[signal_id]->SetFillColor(kBlue);
      hit_model[signal_id]->Draw();
    }
  }

  // plot legendre curve, find the intersection
  Double_t theta, r, x_0, y_0; // parameters used in legendre curve
  TH2F *plot_map = new TH2F("plot_map", "plot_map", 800, 0, 4, 180, -900, 900);
  track_base->cd();
  for (Int_t signal_id = 0; signal_id < event_signal_length; signal_id++) {
    if (event_signal[signal_id][0] == 4) {
      GetHitInfo(event_signal[signal_id][1], event_signal[signal_id][2], &x_0, 
                 &y_0);
      for (Int_t theta_id = 0; theta_id != 10000; theta_id++) {
        theta = M_PI * theta_id / 10000;
        r = LegendreUpperCurve(theta, x_0, y_0, drift_distance[signal_id]);
        plot_map->Fill(theta, r);
        r = LegendreLowerCurve(theta, x_0, y_0, drift_distance[signal_id]);
        plot_map->Fill(theta, r);
      }
    }
  }
  Int_t max_bin_theta, max_bin_r, max_bin_z;
  plot_map->GetMaximumBin(max_bin_theta, max_bin_r, max_bin_z);
  Double_t line_para_k, line_para_b; // line parameters for track line
  line_para_k = -1 / tan(max_bin_theta * 4.0 / 800.0);
  line_para_b = (-900 + max_bin_r * 1800 / 180.0) / 
                sin(max_bin_theta * 4.0 / 800.0);
  //cout << max_bin_theta << " " << max_bin_r << " " << max_bin_z << endl;
  cout << "line parameter b = " << line_para_b << endl 
       << "line parameter k = " << line_para_k << endl;
  TLine *track_line = new TLine(0,  line_para_b,
                                900, line_para_b + line_para_k * 900);
  track_line->SetLineColor(kGreen);
  track_line->Draw("same");

  // find the accurate intersection
  theta = 0;
  r = 0;
  x_0 = 0;
  y_0 = 0;
  Double_t min_theta_limit, max_theta_limit;
  Double_t min_r_limit, max_r_limit;
  min_theta_limit = (max_bin_theta - 10) * 4.0 / 800.0;
  if (min_theta_limit < 0) min_theta_limit = 0;
  max_theta_limit = (max_bin_theta + 10) * 4.0 / 800.0;
  if (max_theta_limit > 4) max_theta_limit = 4;
  min_r_limit = -900 + (max_bin_r - 2) * 1800 / 180.0;
  if (min_r_limit < -800) min_theta_limit = -800;
  max_r_limit = -900 + (max_bin_r + 2) * 1800 / 180.0;
  if (max_r_limit > 800) max_r_limit = 800;
  TH2F *plot_map_accurate = new TH2F("plot_map_accurate", "plot_map_accurate", 
                                     40, min_theta_limit, max_theta_limit, 
                                     50, min_r_limit, max_r_limit);
  for (Int_t signal_id = 0; signal_id < event_signal_length; signal_id++) {
    if (event_signal[signal_id][0] == 4) {
      for (Int_t theta_id = 0; theta_id != 10000; theta_id++) {
        theta = min_theta_limit + 
                theta_id * (max_theta_limit - min_theta_limit) / 10000;
        GetHitInfo(event_signal[signal_id][1], event_signal[signal_id][2], &x_0,
                   &y_0);
        r = LegendreUpperCurve(theta, x_0, y_0, drift_distance[signal_id]);
        plot_map_accurate->Fill(theta, r);
        r = LegendreLowerCurve(theta, x_0, y_0, drift_distance[signal_id]);
        plot_map_accurate->Fill(theta, r);
      }
    }
  }
  max_bin_theta = 0;
  max_bin_r = 0;
  max_bin_z = 0;
  plot_map_accurate->GetMaximumBin(max_bin_theta, max_bin_r, max_bin_z);
  line_para_k = -1 / tan(min_theta_limit + max_bin_theta * (max_theta_limit - 
                min_theta_limit) / 40);
  line_para_b = (min_r_limit + max_bin_r * (max_r_limit - min_r_limit) / 50) / 
                sin(min_theta_limit + max_bin_theta * (max_theta_limit - min_theta_limit) / 40);

  TLine *track_line_accurate = new TLine(0, line_para_b, 900, 
                                         line_para_b + line_para_k * 900);
  track_line_accurate->SetLineColor(kViolet);
  track_line_accurate->Draw("same");

  // find the point near the line found by legendremethod to perform linear fitting
  Double_t hit_position_x[100];
  Double_t hit_position_y[100];
  bool hit_position_acception[100];
  Double_t x_sum = 0;
  Double_t y_sum = 0;
  Double_t xy_sum = 0;
  Double_t xx_sum = 0;
  Double_t quantity = 0;
  for (Int_t signal_id = 0; signal_id < event_signal_length; signal_id++) {
    Double_t possible_hit_position_x1;
    Double_t possible_hit_position_x2;
    Double_t possible_hit_position_y1;
    Double_t possible_hit_position_y2;
    Double_t x0, y0;
    GetHitInfo(event_signal[signal_id][1], event_signal[signal_id][2], &x0, 
               &y0);
    possible_hit_position_x1 = x0 - drift_distance[signal_id] * line_para_k / sqrt(line_para_k * line_para_k + 1);
    possible_hit_position_x2 = x0 + drift_distance[signal_id] * line_para_k / sqrt(line_para_k * line_para_k + 1);
    possible_hit_position_y1 = y0 + drift_distance[signal_id] / sqrt(line_para_k * line_para_k + 1);
    possible_hit_position_y2 = y0 - drift_distance[signal_id] / sqrt(line_para_k * line_para_k + 1);

    if (GetPointLineDistance(possible_hit_position_x1, possible_hit_position_y1, line_para_k, line_para_b) > GetPointLineDistance(possible_hit_position_x2, possible_hit_position_y2, line_para_k, line_para_b)) {
      hit_position_x[signal_id] = possible_hit_position_x2;
      hit_position_y[signal_id] = possible_hit_position_y2;
    } else if (GetPointLineDistance(possible_hit_position_x1, possible_hit_position_y1, line_para_k, line_para_b) < GetPointLineDistance(possible_hit_position_x2, possible_hit_position_y2, line_para_k, line_para_b)) {
      hit_position_x[signal_id] = possible_hit_position_x1;
      hit_position_y[signal_id] = possible_hit_position_y1;
    }

    if (fabs(GetPointLineDistance(hit_position_x[signal_id], hit_position_y[signal_id], line_para_k, line_para_b) - drift_distance[signal_id]) < 2) {
      hit_position_acception[signal_id] = true;
      quantity += 1.0;
      x_sum += hit_position_x[signal_id];
      y_sum += hit_position_y[signal_id];
      xy_sum += hit_position_x[signal_id] * hit_position_y[signal_id];
      xx_sum += hit_position_x[signal_id] * hit_position_x[signal_id];
    } else {
      hit_position_acception[signal_id] = false;
    }
  }
  
  Double_t x_mean = 0;
  Double_t y_mean = 0;
  Double_t xy_mean = 0;
  Double_t xx_mean = 0;
  if (quantity != 0) {
    x_mean =  x_sum / quantity;
    y_mean =  y_sum / quantity;
    xy_mean =  xy_sum / quantity;
    xx_mean =  xx_sum / quantity;
  } else if (quantity == 0) {
    x_mean =  x_sum;
    y_mean =  y_sum;
    xy_mean =  xy_sum;
    xx_mean =  xx_sum;
  }
  
  Double_t final_k = (xy_mean - x_mean * y_mean) / (xx_mean - x_mean * x_mean);
  Double_t final_b = y_mean - final_k * x_mean;

  cout << "accurate line parameter b = " << final_b << endl 
       << "accurate line parameter k = " << final_k << endl;
  TLine *track_line_final = new TLine(0, final_b, 900, 
                                         final_b + final_k * 900);
  track_line_final->SetLineColor(kPink);
  track_line_final->Draw("same");

  // calculate segment residual
  cout << endl << "// Segment Residual ///////////////////////" << endl;
  for (Int_t signal_id = 0; signal_id < event_signal_length; signal_id++) {
    if (event_signal[signal_id][0] == 4) {
      Double_t x_0, y_0, distance, residual;
      Int_t hit_layer, hit_column;
      GetHitInfo(event_signal[signal_id][1], event_signal[signal_id][2], &x_0, 
                 &y_0);
      GetHitLayerColumn(event_signal[signal_id][1], event_signal[signal_id][2],
                        &hit_layer, &hit_column);
      distance = GetPointLineDistance(x_0, y_0, final_k, final_b);
      residual = drift_distance[signal_id] - distance;
      cout << "hit layer: " << hit_layer << " | " << "hit column: " << hit_column << " | residual: " << residual << endl;
    }
  }
  cout << "///////////////////////////////////////////" << endl;

  // draw the plots for debugging
  TCanvas *track_map_canvas = new TCanvas("track", "track", 0, 550, 600, 600);
  track_map_canvas->cd();
  plot_map->SetStats(0);
  plot_map->Draw("COLZ");

  TCanvas *track_map_accurate_canvas = new TCanvas("track_accurate", "track_accurate", 600, 550, 600, 600);
  track_map_accurate_canvas->cd();
  plot_map_accurate->SetStats(0);
  plot_map_accurate->Draw("COLZ");

  return 0;
}
// end PlotTrack ///////////////////////////////////////////////////////////////

Double_t GetPointLineDistance(Double_t x, Double_t y, Double_t k, Double_t b) {
  return fabs((k * x - y +b) / sqrt(k * k + 1.0));
}