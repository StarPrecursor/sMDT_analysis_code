/*******************************************************************************
  file name: CircleDrawingTest.cxx
  author: Zhe Yang
  created: 01/25/2019
  last modified: 02/05/2019

  description: 
  Test code for tube ploting function in DecodeRawData.cxx
*******************************************************************************/

#include <stdio.h>
#include <iostream>
#include <bitset>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TH1.h"

void CircleDrawingTest() {
  
  const Double_t layer_distance = 13.077;
  const Double_t column_distance = 15.1;
  const Double_t radius = 7.5;

  TFile *test_rootfile = new TFile("test.root", "RECREATE");
  TCanvas *track_base = new TCanvas("track base", "track base", 1200, 480);
  track_base->cd();
  //track_base->Range(0, 0, 820, 400);
  Double_t center_x, center_y;
  Double_t track_corner_x[2] = {0, 800}; 
  Double_t track_corner_y[2] = {0, 320}; 
  TGraph * track_baseline = new TGraph(2, track_corner_x, track_corner_y);
  track_baseline->SetTitle("run");
  track_baseline->Draw("AP");
  
      TEllipse *tube_model[53][8];
      for (Int_t layer_id = 0; layer_id != 4; layer_id++) {
          for (Int_t column_id = 0; column_id != 53; column_id++) {
            center_x = 7.5 + column_id * column_distance 
                       + (layer_id % 2) * column_distance / 2.0;
            center_y = 7.5 + layer_id * layer_distance;
            tube_model[layer_id][column_id] = new TEllipse(center_x, center_y,
                                                           radius, radius);
            tube_model[layer_id][column_id]->Draw();
          }
      }
      for (Int_t layer_id = 4; layer_id != 8; layer_id++) {
          for (Int_t column_id = 0; column_id != 53; column_id++) {
            center_x = 7.5 + column_id * column_distance 
                       + (layer_id % 2) * column_distance / 2.0;
            center_y = 7.5 + (layer_id - 4) * layer_distance + 224.231;
            tube_model[layer_id][column_id] = new TEllipse(center_x, center_y,
                                                           radius, radius);
            tube_model[layer_id][column_id]->Draw();
          }
      }

      TEllipse *hit_model[128];
    for (int i = 0 ; i < 3; i++) {
      
      test_rootfile->cd();
      track_base->cd();

      for (int j = 0; j < 20; j++) {
        hit_model[j] = new TEllipse(center_x - 15.1 * j, center_y, radius, radius);
        hit_model[j]->SetFillColor(kRed);
        hit_model[j]->Draw();
      }

      
      
      TCanvas *new_canvas = new TCanvas("new canvas","new canvas",1200,480);
      new_canvas->Divide(1, 1);
      new_canvas->cd(1);
      track_base->DrawClonePad();



      test_rootfile->WriteTObject(new_canvas);
      delete new_canvas;
      
      for (int j = 0; j < 20; j++) {
        delete hit_model[j];
      }

      //tube_model[30][4]->Draw();
      //hit_model[31][5]->Draw();
    }
}