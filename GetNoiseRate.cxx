/*******************************************************************************
  file name: GetNoiseRate.cxx
  author: Zhe Yang
  created: 02/19/2019
  last modified: 03/05/2019

  description:
  -Calculate the noise rate for each tube using the formula:
  noise_rate  = hits_in_searchwindow / (searchwindow * trigger_number) 

  remark:
  -search window is set to (0, 70ns)

*******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TLine.h"
#include "TTree.h"

#include "DecodeRawData.h"
#include "FitTDCSlope.cxx"
#include "GetHitInfo.h"
#include "GetTrack.h"

#define _USE_MATH_DEFINES
//#define GETNOISERATE_DEBUG

//string s_input_filename = "run00187664_20190222.dat";
string s_input_filename = "run00187691_20190301.dat";

using namespace std;

int GetNoiseRate() {
  // stop canvas from popping up
  Bool_t batch_mode_flag;
  if (gROOT->IsBatch()) {
    batch_mode_flag = true;
  } else {
    batch_mode_flag = false;
  }
  gROOT->SetBatch(kTRUE);

  // read input rootfile and set branch
  char input_path_name[256];
  char input_directory_name[256];
  sprintf(input_path_name, s_input_filename.c_str());
  strcat(input_path_name, ".dir/");
  strcpy(input_directory_name, input_path_name);
  strcat(input_path_name, s_input_filename.c_str());
  strcat(input_path_name, ".out.root");
  TFile *p_input_rootfile = new TFile(input_path_name);
  TTree *p_input_tree = (TTree*)p_input_rootfile->Get("HPTDCData");
  Float_t type, tdc, data0, data1, data2, signal_flag;
  p_input_tree->SetBranchAddress("type", &type);
  p_input_tree->SetBranchAddress("tdc", &tdc);
  p_input_tree->SetBranchAddress("data0", &data0);
  p_input_tree->SetBranchAddress("data1", &data1);
  p_input_tree->SetBranchAddress("data2", &data2);
  p_input_tree->SetBranchAddress("signal_flag", &signal_flag);

  // get drift time distribution histogram and calculate T0
  Long_t total_noise_count = 0;
  Long_t test_noise_count = 0;
  Long_t total_trigger_count = 0;
  Long_t noise_count[MAX_TDC_QUANTITY][MAX_TDC_CHANNEL_QUANTITY] = {{0}};
  Double_t t0_value[MAX_TDC_QUANTITY][MAX_TDC_CHANNEL_QUANTITY];
  Double_t fit_output_value[2][5];
  TH1F *p_drift_time_hist[MAX_TDC_QUANTITY];
  char directory_name[256];
  char hist_name[256];
  char temp_name[256];
  for (Int_t tdc_id = 8; tdc_id < 12; tdc_id++) {
    for (Int_t ch_id = 0; ch_id != MAX_TDC_CHANNEL_QUANTITY; ch_id++) {
      sprintf(directory_name, "TDC_%02d_of_%02d_Time_Spectrum", tdc_id, MAX_TDC_QUANTITY);
      //sprintf(temp_name, "tdc_%d_channel_%d_tdc_time_spectrum", tdc_id, ch_id);
      sprintf(temp_name, "tdc_%d_channel_%d_tdc_time_spectrum_corrected", tdc_id, ch_id);
      //sprintf(temp_name, "tdc_%d_tdc_time_spectrum", tdc_id);
      //sprintf(temp_name, "tdc_%d_tdc_time_spectrum_corrected", tdc_id);
      strcpy(hist_name, directory_name);
      strcat(hist_name, "/");
      strcat(hist_name, temp_name);
      p_input_rootfile->GetObject(hist_name, p_drift_time_hist[tdc_id]);
      FitTDCSlope(p_drift_time_hist[tdc_id], fit_output_value, nullptr, nullptr, kFALSE);
      t0_value[tdc_id][ch_id] = fit_output_value[0][1];
      cout << fit_output_value[0][0] << " | " << fit_output_value[0][1] << " | " << fit_output_value[0][2] << " | " << fit_output_value[0][3] << " | " << fit_output_value[0][4] << endl;
      
      for (Int_t bin_number = 150; bin_number != 180; bin_number++) {
        test_noise_count += p_drift_time_hist[tdc_id]->GetBin(bin_number);
      }
      cout << endl << " count: " << test_noise_count << endl << endl;
    }
    
  }
  for (Int_t tdc_id = 0; tdc_id < MAX_TDC_QUANTITY; tdc_id++) {
    cout << "// #" << tdc_id << " //////////////" << endl;
    for (Int_t ch_id =0; ch_id < MAX_TDC_CHANNEL_QUANTITY; ch_id++) {
      cout << "t0 = " << t0_value[tdc_id][ch_id] << endl;
    }
  }

  //// calculate noise rate accurately
  // prepare output file
  char output_directory_name[256];
  strcpy(output_directory_name, input_directory_name);
  strcat(output_directory_name, "noise_level_output");
  if (mkdir(output_directory_name, 0777) == -1) {
    cerr << strerror(errno) << endl;
  }
  char output_path_name[256];
  strcpy(output_path_name, output_directory_name);
    strcat(output_path_name, s_input_filename.c_str());
    strcat(output_path_name, ".noise.root");
  TFile *p_output_rootfile = new TFile(output_path_name, "RECREATE");
  TNtuple *p_noise_data = new TNtuple("NoiseData", "NoiseData", "event_id:tdc_id:channel_id:coarse:fine");
  TH1F *p_radius_hist = new TH1F("radius_distribution", "radius_distribution", 50, 0, 8);
  TH1F *p_residual_hist = new TH1F("residual_distribution", "residual_distribution", 100, -1000, 1000);
  TH1F *p_residual_hist_vs_radius[8];
  for (Int_t i = 0; i != 8; i++) {
    char hist_name[128];
    sprintf(hist_name, "residual_distribution_%dmm", i);
    p_residual_hist_vs_radius[i] = new TH1F(hist_name, hist_name, 100, -1000, 1000);
  }

  // find the event and get event's data
  Int_t total_entries = p_input_tree->GetEntries();
  Int_t event_id = 0;
  Int_t current_event_id = -1;
  Int_t event_length = 0;
  Int_t event_trigger_length = 0;
  Int_t event_signal_length = 0;
  Int_t event_trigger[128][6];
  Int_t event_signal[128][6];
  Double_t output_line_parameter_k;
  Double_t output_line_parameter_b;
  Double_t current_time = 0;
  bool event_start_flag = false;
  bool good_track_flag = false;
  Int_t output_good_hit_flag[128];
  Double_t output_residual[128];
  Double_t output_drift_time[128];
  Double_t output_drift_distance[128];
  #ifdef GETNOISERATE_DEBUG
  gROOT->SetBatch(kFALSE);
  #endif
  for (Int_t entry_id = 0; entry_id < /*500000*/ total_entries; entry_id++) {
    p_input_tree->GetEntry(entry_id);

    if (type == 0 || type == 2 || type == 3) {
      event_id = data1;
    } else if (type == 4 || type == 5) {
      if (event_length < 100) {
        if (tdc == 1) {
          event_trigger[event_trigger_length][0] = type;
          event_trigger[event_trigger_length][1] = tdc;
          event_trigger[event_trigger_length][2] = data0;
          event_trigger[event_trigger_length][3] = data1;
          event_trigger[event_trigger_length][4] = data2;
          event_trigger[event_trigger_length][5] = signal_flag;
          event_trigger_length++;
        } else {
          if (signal_flag) {
            event_signal[event_signal_length][0] = type;
            event_signal[event_signal_length][1] = tdc;
            event_signal[event_signal_length][2] = data0;
            event_signal[event_signal_length][3] = data1;
            event_signal[event_signal_length][4] = data2;
            event_signal[event_signal_length][5] = signal_flag;
            event_signal_length++;
          }
        }
        event_length++;
      }
    }

    if (current_event_id != event_id) {
      cout << "Current event id: " << current_event_id << endl;
      good_track_flag = GetTrack(event_trigger_length, event_trigger,
                                 event_signal_length, event_signal, t0_value,
                                 &output_line_parameter_k,
                                 &output_line_parameter_b,
                                 output_drift_time,
                                 output_drift_distance,
                                 output_residual,
                                 output_good_hit_flag);
      #ifdef GETNOISERATE_DEBUG
      // debugging
      if (good_track_flag) {
        cout << "Paused. Enter 'e' to exit or other key to continue." << endl;
        if (getchar() == 'e') {
          cout << "Program is stopped by user." << endl;
          return 1;
        }
      }
      #endif

      for (Int_t trigger_id = 0; trigger_id < event_trigger_length; 
           trigger_id++) {
        if (event_trigger[trigger_id][0] == 4) {
          total_trigger_count++;
        }
      }

      for (Int_t signal_id = 0; signal_id < event_signal_length; signal_id++) {
        if (output_good_hit_flag[signal_id] == -1) {
          if (event_signal[signal_id][0] == 4
              && event_signal[signal_id][1] != 1) {
            current_time = (event_signal[signal_id][3] +
                            event_signal[signal_id][4] / 128.0 ) * 
                            25.0;
            if (current_time > 0 && current_time < 70) {
              noise_count[event_signal[signal_id][1]][event_signal[signal_id][2]]++;
              total_noise_count++;
            }

            p_noise_data->Fill(current_event_id,
                               event_signal[signal_id][1],
                               event_signal[signal_id][2],
                               event_signal[signal_id][3],
                               event_signal[signal_id][4]);
          }
        } else if (output_good_hit_flag[signal_id] != -1) {
          if (event_signal[signal_id][0] == 4
              && event_signal[signal_id][1] != 1) {
            p_residual_hist->Fill(output_residual[signal_id] * 1000);
            p_residual_hist_vs_radius[output_good_hit_flag[signal_id]]->Fill(output_residual[signal_id] * 1000);
            p_radius_hist->Fill(output_drift_distance[signal_id]);
          }
        }
      }

      event_trigger_length = 0;
      event_signal_length = 0;
      event_length = 0;

      current_event_id = event_id;
    }
  }

  // calculate noise rate
  Double_t channel_noise_rate = 0;
  Double_t test_noise_rate;
  for (Int_t tdc_id = 0; tdc_id < MAX_TDC_QUANTITY; tdc_id++) {
    for (Int_t ch_id = 0; ch_id < MAX_TDC_CHANNEL_QUANTITY; ch_id++) {
      channel_noise_rate = noise_count[tdc_id][ch_id] / (0.00000007 * total_trigger_count);
      if (tdc_id >= 8 && tdc_id <= 11) {
        cout << "tdc#" << tdc_id << " ch#" << ch_id << " noise_rate: " << channel_noise_rate << endl;
      }
      test_noise_rate += channel_noise_rate;
    }
  }

  Double_t noise_rate = total_noise_count / (0.00000007 * total_trigger_count);
  cout << "noise counts: " << total_noise_count << endl;
  cout << "trigger counts: " << total_trigger_count << endl;
  cout << "total noise rate is: " << noise_rate << " Hz" << endl;
  cout << "test noise rate is: " << test_noise_rate << " Hz" << endl;

  p_residual_hist->Draw();
  p_output_rootfile->Write();

  // resume batch mode setup
  if (batch_mode_flag == false) {
    gROOT->SetBatch(kFALSE); 
  }

  return 0;
}