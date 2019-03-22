/*******************************************************************************
  file name: DecodeRawData.cxx
  author: Zhe Yang
  created: 01/25/2019
  last modified: 02/14/2019

  description:
  -Decode .raw data from HPTDC and save data to ntuple.
  -Use "root -l -b" command to run the code when debugging is not needed

  remark:
  -Learned basic decode method from Shuzhou Zhang, redeveloped and added new 
  function for new HPTDC data format.

*******************************************************************************/

#include "DecodeRawData.h"

#include <stdio.h>
#include <iostream>
#include <bitset>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h> 
#include <sys/types.h> 

#include "TFile.h"
#include "TDirectory.h"
#include "TNtuple.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "CheckEvent.cxx"
#include "CollectOneEvent.cxx"
#include "GetEventTrigger.cxx"
#include "GetHitInfo.h"
#include "GetHitLayerColumn.h"
#include "GetTrackFast.h"
#include "MdtSlewCorFuncHardcoded.cxx"

//#define DECODER_DEBUG // comment this when debugging information is not needed
//#define SET_MAXWORDS // comment this if you want to decode the whole data words
//#define SAVE_TRACKS_OUT_OF_ROOT // comment this to save decoding time, 
                                // otherwise time-consuming for picture-saving 
                                // is high

using namespace std;

double calibrated_p(const double &adc);
double ADC_Correction(double adc, vector<Double_t> v_LUT);

int DecodeRawData() {
  Bool_t batch_mode_flag;
  if (gROOT->IsBatch()) {
    batch_mode_flag = true;
  } else {
    batch_mode_flag = false;
  }
  gROOT->SetBatch(kTRUE); // stop canvas from popping up

  // open data file
  char filename[256];
  char input_filename[256];
  cout << "Please enter the input file name to analyse:" << endl;
  cin >> filename;
  strcpy(input_filename, "./");
  strcat(input_filename, filename);
  ifstream data_in_flow;
  data_in_flow.open(input_filename);
  while (!data_in_flow.is_open()) {
    cout << "Unable to open input file!" << endl;
    cout << "Please enter the input file name to analyse, enter 'NULL' to exit:" << endl;
    cin >> filename;
    if (strcmp(filename, "NULL") == 0) {
      cout << "NULL input file, exiting program." << endl;
      return 1;
    }
    strcpy(input_filename, "./");
    strcat(input_filename, filename);
    data_in_flow.open(input_filename);
  }
  data_in_flow.seekg(0, data_in_flow.end);
  unsigned int data_in_flow_length = data_in_flow.tellg();
  data_in_flow.seekg(0, data_in_flow.beg);

  // prepare output file
  char output_directoryname[256];
  strcpy(output_directoryname, filename);
  strcat(output_directoryname, ".dir");
  if (mkdir(output_directoryname, 0777) == -1) {
    cerr << strerror(errno) << endl; 
  }
  chdir(output_directoryname);
  char output_filename[256];
  strcpy(output_filename, filename);
  strcat(output_filename, ".out");
  ofstream data_out_flow;
  data_out_flow.open(output_filename);
  if (!data_out_flow.is_open()){
    cout << "Unable to open output data file!" << endl;
    return 1;
  }

  // prepare output Ntuple
  char output_root_filename[200];
  strcpy(output_root_filename, output_filename);
  strcat(output_root_filename, ".root");
  TFile *p_output_rootfile = new TFile(output_root_filename, "RECREATE");
  const char *ntuple_varlist = "type:tdc:data0:data1:data2:signal_flag";
    // for header_type = 2, data0 = 0, data1 = EventID, data2 = BunchID
    // for heaser_type = 4/5, data0 = channel, data1 = coarse, data2 = fine
  TNtuple *p_HPTDC_data = new TNtuple("HPTDCData", "HPTDCData", ntuple_varlist);

  // prepare output histogram
  TH1F *p_leading_time = new TH1F("leading time spectrum", 
                                  "leading time spectrum", 100, 0, 1000);
  TH1F *p_trailing_time = new TH1F("trailing time spectrum", 
                                   "trailing time spectrum", 100, 0, 1000);
  TH1F *p_hits_distribution[MAX_TUBE_LAYER];
  char histogram_name[256];
  for (Int_t layer_id = 0; layer_id != MAX_TUBE_LAYER; layer_id++) {
    sprintf(histogram_name, "layer_%d_hits_distribution", layer_id);
    p_hits_distribution[layer_id] = new TH1F(histogram_name, histogram_name, 54,
                                             0, 54);
  }
  TH1F *p_tdc_time[MAX_TDC_QUANTITY][MAX_TDC_CHANNEL_QUANTITY];
  TH1F *p_tdc_time_original[MAX_TDC_QUANTITY][MAX_TDC_CHANNEL_QUANTITY];
  TH1F *p_tdc_time_corrected[MAX_TDC_QUANTITY][MAX_TDC_CHANNEL_QUANTITY];
  TH1F *p_adc_time[MAX_TDC_QUANTITY][MAX_TDC_CHANNEL_QUANTITY];
  TH1F *p_tdc_tdc_time[MAX_TDC_QUANTITY];
  TH1F *p_tdc_tdc_time_original[MAX_TDC_QUANTITY];
  TH1F *p_tdc_tdc_time_corrected[MAX_TDC_QUANTITY];
  TH1F *p_tdc_adc_time[MAX_TDC_QUANTITY];
  TH1F *p_tdc_channel[MAX_TDC_QUANTITY];
  TDirectory *tdc_directory[MAX_TDC_QUANTITY];
  char directory_name[256];
  for (Int_t tdc_id = 0; tdc_id != MAX_TDC_QUANTITY; tdc_id++) {
    sprintf(directory_name, "TDC_%02d_of_%02d_Time_Spectrum", tdc_id,
            MAX_TDC_QUANTITY);
    tdc_directory[tdc_id] = p_output_rootfile->mkdir(directory_name);
    tdc_directory[tdc_id]->cd();
    if (mkdir(directory_name, 0777) == -1) {
      cerr << strerror(errno) << endl; 
    }

    sprintf(histogram_name, "tdc_%d_tdc_time_spectrum", tdc_id);
    p_tdc_tdc_time[tdc_id] = new TH1F(histogram_name, histogram_name, 300, -500,
                                      500);
    sprintf(histogram_name, "tdc_%d_tdc_time_spectrum_original", tdc_id);
    p_tdc_tdc_time_original[tdc_id] = new TH1F(histogram_name, histogram_name, 
                                               300, -500, 500);
    sprintf(histogram_name, "tdc_%d_tdc_time_spectrum_corrected", tdc_id);
    p_tdc_tdc_time_corrected[tdc_id] = new TH1F(histogram_name, histogram_name, 
                                               300, -500, 500);
    sprintf(histogram_name, "tdc_%d_adc_time_spectrum", tdc_id);
    p_tdc_adc_time[tdc_id] = new TH1F(histogram_name, histogram_name, 300, -500,
                                      500);
    sprintf(histogram_name, "tdc_%d_channel_distribution", tdc_id);
    p_tdc_channel[tdc_id] = new TH1F(histogram_name, histogram_name, 24, 0, 24);
    for (Int_t channel_id = 0; channel_id != MAX_TDC_CHANNEL_QUANTITY; 
         channel_id++) {
      sprintf(histogram_name, "tdc_%d_channel_%d_tdc_time_spectrum", tdc_id, 
              channel_id);
      p_tdc_time[tdc_id][channel_id] = new TH1F(histogram_name, histogram_name, 
                                                450, -500, 1000);
      sprintf(histogram_name, "tdc_%d_channel_%d_tdc_time_spectrum_original", 
              tdc_id, channel_id);
      p_tdc_time_original[tdc_id][channel_id] = new TH1F(histogram_name,
                                                         histogram_name, 450,
                                                        -500, 1000);
      sprintf(histogram_name, "tdc_%d_channel_%d_tdc_time_spectrum_corrected", 
              tdc_id, channel_id);
      p_tdc_time_corrected[tdc_id][channel_id] = new TH1F(histogram_name,
                                                         histogram_name, 450,
                                                        -500, 1000);
      sprintf(histogram_name, "tdc_%d_channel_%d_adc_time_spectrum", tdc_id, 
              channel_id);
      p_adc_time[tdc_id][channel_id] = new TH1F(histogram_name, histogram_name, 
                                                300, -500, 500);
    }
  }

    // prepare base of output track display
    const Double_t layer_distance = 13.0769836;
    const Double_t column_distance = 15.1;
    const Double_t radius = 7.5;

    TDirectory *event_track[10];
    event_track[0] = p_output_rootfile->mkdir("event_tracks_record_example0");
    if (mkdir("event_tracks_record_example0", 0777) == -1) {
        cerr << strerror(errno) << endl; 
    }
    char track_group_name[128];
    strcpy(track_group_name, "event_tracks_record_example0");
    TCanvas *track_base = new TCanvas("track base", "track base", 1200, 480);
    track_base->cd();
    Double_t center_x, center_y;
    Double_t track_corner_x[2] = {0, 800};
    Double_t track_corner_y[2] = {0, 320};
    TGraph * track_baseline = new TGraph(2, track_corner_x, track_corner_y);
    track_baseline->SetTitle("event");
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

    // prepare adc correction table
    vector<Double_t> v_LUT;
    for (Int_t i = 0; i != 400; i++) {
      v_LUT.push_back(109./calibrated_p(i));
    }

    // execute selected operation and save data to output file
    unsigned int word, id, header_type;
    unsigned int trigger_cache_length = 0;
    unsigned int trigger_cache[128][5];
    unsigned int signal_cache_length = 0;
    unsigned int signal_cache[128][5];
    unsigned int new_data_line[5];
    bool signal_flag[128];
    bool got_new_event = false;
    unsigned long total_loop = 0;
    unsigned long total_tracks = 0;
    int current_track_group = 0;
    unsigned int total_events;
    TCanvas *hit_display_cache;
    unsigned int event_count = 0;
    unsigned long current_event_id = 0;
    unsigned long event_id_group = 0;
    bool new_trigger_flag = false;
    bool new_signal_flag = false;
    bool plot_flag = true;
    cout << "Processing..." << endl;

    #ifdef SET_MAXWORDS
    while (data_in_flow.read((char *) &word, sizeof(word))
           && total_loop < MAX_WORD_QUANTITY) {
    #else
    while (data_in_flow.read((char *) &word, sizeof(word))) {
    #endif
    bitset<32> data_cache;
    bitset<4> header; // for header type
    data_cache = word;
    header = word >> 28; // get the four bits header of this word
    header_type = static_cast<unsigned int>((header.to_ulong()));

    unsigned int _type, _event_id, _bunch_id,  _tdc, _channel, _width, _coarse, _fine ,_leading_time;
    Double_t time_in_ns;
    bitset<4> type;
    bitset<4> tdc;
    bitset<12> event_id;
    bitset<12> bunch_id;
    bitset<5> channel;
    bitset<12> coarse;
    bitset<7> fine;

    unsigned int leading_mean_coarse = 0;
    unsigned int trailing_mean_coarse = 0;
    bool pass_event_check = false;
    unsigned int trigger_data_line[5];
    bool got_trigger = false;
    if (header_type == 2 || header_type == 3 || header_type == 0) {
      type = word >> 28;
      tdc = word >> 24;
      event_id = word >> 12;
      bunch_id = word;
      _type = static_cast<unsigned int>((type.to_ulong()));
      _tdc = static_cast<unsigned int>((tdc.to_ulong()));
      _event_id = static_cast<unsigned int>((event_id.to_ulong()));
      _bunch_id = static_cast<unsigned int>((bunch_id.to_ulong()));

      // analysis data if one event is over
      got_new_event = (current_event_id != _event_id);
      if (got_new_event) {
        event_count++;
        if (current_event_id == _event_id + 4096 * event_id_group + 4095) {
          event_id_group++;
        }
        current_event_id = _event_id + 4096 * event_id_group;
        leading_mean_coarse = 0;
        trailing_mean_coarse = 0;
        pass_event_check = false;
        got_trigger = false;

        #ifdef DECODER_DEBUG
        cout << "trigger: " << endl;
        for (int trigger_id = 0; trigger_id < trigger_cache_length; 
             trigger_id++) {
          cout << trigger_id << " : ";
          for (int i = 0; i < 5; i++) {
            cout << trigger_cache[trigger_id][i] << " ";
            }
          cout << endl;
        }
        cout << "signal: " << endl;
        for (int signal_id = 0; signal_id < signal_cache_length; signal_id++) {
          cout << signal_id << " : ";
          for (int i = 0; i < 5; i++) {
            cout << signal_cache[signal_id][i] << " ";
          }
          cout << endl;
        }
        cout << "Paused. Enter 'e' to exit or other key to continue." << endl;
        if (getchar() == 'e') {
          cout << "Program is stopped by user." << endl;
          return 1;
        }
        #endif

        // get trigger id, calculate tdc time & time difference
        // between leading & trailing edge
        int selected_trigger_id = 128;
        long current_coarse_difference = 2048;
        if (trigger_cache_length != 0) {
          long temp_coarse = 0;
          long total_coarse = 0;
          long temp_mean_coarse = 0;
          int trigger_quantity = 0;
          for (int trigger_id = 0; trigger_id < trigger_cache_length;
               trigger_id++) {
            if (trigger_cache[trigger_id][0] == 4) {
              total_coarse += trigger_cache[trigger_id][3];
              trigger_quantity++;
            }
          }
          if (trigger_quantity == 0) {
            temp_mean_coarse = 0;
          } else {
            temp_mean_coarse = total_coarse / trigger_quantity;
          }
          for (int trigger_id = 0; trigger_id < trigger_cache_length;
               trigger_id++) {
            if (trigger_cache[trigger_id][0] == 4) {
              temp_coarse = trigger_cache[trigger_id][3];
              if (abs(temp_coarse - temp_mean_coarse) < 
                  current_coarse_difference) {
                current_coarse_difference = abs(temp_coarse - temp_mean_coarse);
                selected_trigger_id = trigger_id;
                // selected_trigger_id = 0;
                got_trigger = true;
              }
              double tdc_time = (trigger_cache[trigger_id][3] + trigger_cache[trigger_id][4] / 128.0) * 25.0;
              p_tdc_time[trigger_cache[trigger_id][1]][trigger_cache[trigger_id][2]]->Fill(tdc_time);
              p_tdc_tdc_time[trigger_cache[trigger_id][1]]->Fill(tdc_time);
              p_tdc_channel[trigger_cache[trigger_id][1]]->Fill(trigger_cache[trigger_id][2]);

              double tdc_time2;
              for (int trigger_id2 = 0; trigger_id2 < trigger_cache_length; trigger_id2++) {
                if (trigger_cache[trigger_id2][0] == 5) {
                  if (trigger_cache[trigger_id][1] == trigger_cache[trigger_id2][1] && trigger_cache[trigger_id][2] == trigger_cache[trigger_id2][2]) {
                    tdc_time2 = (trigger_cache[trigger_id2][3] + trigger_cache[trigger_id2][4] / 128.0) * 25.0;
                    Double_t adc_time = tdc_time2 - tdc_time;
                    p_adc_time[trigger_cache[trigger_id][1]][trigger_cache[trigger_id][2]]->Fill(adc_time);
                    p_tdc_adc_time[trigger_cache[trigger_id][1]]->Fill(adc_time);
                    break;
                  }
                }
              }
            }
                        
            #ifdef DECODER_DEBUG
            cout << " current_coarse_difference= " << current_coarse_difference;
            cout << "  selected_trigger_id=" << selected_trigger_id << endl;
            #endif
          }

          trigger_data_line[0] = trigger_cache[selected_trigger_id][0];
          trigger_data_line[1] = trigger_cache[selected_trigger_id][1];
          trigger_data_line[2] = trigger_cache[selected_trigger_id][2];
          trigger_data_line[3] = trigger_cache[selected_trigger_id][3];
          trigger_data_line[4] = trigger_cache[selected_trigger_id][4];

          // calculate tdc time using leading edge
          int hit_layer, hit_column;
          for (int signal_id = 0; signal_id < signal_cache_length;
               signal_id++) {
            if (signal_cache[signal_id][0] == 4) {
              double current_signal_time = (signal_cache[signal_id][3] + signal_cache[signal_id][4]/ 128.0 ) * 25.0;
              double current_trigger_time = (trigger_data_line[3] + trigger_data_line[4]/ 128.0 ) * 25.0;
              double current_tdc_time = current_signal_time - current_trigger_time;
              if (got_trigger) {
                p_tdc_time[signal_cache[signal_id][1]][signal_cache[signal_id][2]]->Fill(current_tdc_time);
                p_tdc_time_original[signal_cache[signal_id][1]][signal_cache[signal_id][2]]->Fill(current_signal_time);
                p_tdc_tdc_time[signal_cache[signal_id][1]]->Fill(current_tdc_time);
                p_tdc_tdc_time_original[signal_cache[signal_id][1]]->Fill(current_signal_time);
                p_tdc_channel[signal_cache[signal_id][1]]->Fill(signal_cache[signal_id][2]);
              }
              GetHitLayerColumn(signal_cache[signal_id][1], signal_cache[signal_id][2], &hit_layer, &hit_column);
              p_hits_distribution[hit_layer]->Fill(hit_column);

              // calculate time difference between leading & trailing edge
              double adc_time = 0;
              bool got_adc = false;
              for (int signal_id2 = 0; signal_id2 < signal_cache_length; 
                   signal_id2++) {
                if ((signal_cache[signal_id2][0] == 5) && 
                    (signal_cache[signal_id2][1] == signal_cache[signal_id][1]) && (signal_cache[signal_id2][2] == signal_cache[signal_id][2])) {
                  adc_time = (signal_cache[signal_id2][3] + signal_cache[signal_id2][4] / 128.0 ) * 25.0 - (signal_cache[signal_id][3] + signal_cache[signal_id][4] / 128.0 ) * 25.0;

                  p_adc_time[signal_cache[signal_id][1]][signal_cache[signal_id][2]]->Fill(adc_time);
                  p_tdc_adc_time[signal_cache[signal_id][1]]->Fill(adc_time);
                  got_adc = true;
                }
              }

              if (got_trigger && got_adc) {
                double tdc_time_corrected = current_tdc_time - ADC_Correction(adc_time, v_LUT);
                p_tdc_time_corrected[signal_cache[signal_id][1]][signal_cache[signal_id][2]]->Fill(tdc_time_corrected);
                p_tdc_tdc_time_corrected[signal_cache[signal_id][1]]->Fill(tdc_time_corrected);
              }

              if (got_adc && adc_time > 70) {
                signal_flag[signal_id] = true;
              }
              else {
                signal_flag[signal_id] = false;
              }
            }
          }

          // make plot for event display
          pass_event_check = CheckEvent(signal_cache_length, signal_cache, &leading_mean_coarse, &trailing_mean_coarse);
          if (pass_event_check) {
            int temp_track_group;
            unsigned int tdc_id, channel_id;
            double hit_x, hit_y;
            temp_track_group = total_tracks / 100;
            if (current_track_group < temp_track_group) {
              if (temp_track_group < 10) {
                sprintf(track_group_name, "event_tracks_record_example%d", 
                        temp_track_group);
                event_track[temp_track_group] = p_output_rootfile->mkdir(track_group_name);

                #ifdef SAVE_TRACKS_OUT_OF_ROOT
                  if (mkdir(track_group_name, 0777) == -1) {
                    cerr << strerror(errno) << endl; 
                  }
                #endif
              }
              else {
                plot_flag = false;
              }
              current_track_group = temp_track_group;
            }
            if (plot_flag) {
              event_track[temp_track_group]->cd();
              track_base->cd();
              for (int signal_id = 0; signal_id < signal_cache_length;
                   signal_id++) {
                if (signal_cache[signal_id][0] == 4) {
                  tdc_id = signal_cache[signal_id][1];
                  channel_id = signal_cache[signal_id][2];
                  GetHitInfo(tdc_id, channel_id, &hit_x, &hit_y);
                  hit_model[signal_id] = new TEllipse(hit_x, hit_y, radius, radius);
                  hit_model[signal_id]->SetFillColor(kRed);
                  hit_model[signal_id]->Draw();
                }
              }

              char canvas_name[256];
              char canvas_output_name[256];
              cout << "entry" << p_HPTDC_data->GetEntries() << endl;
              cout << "curret_event_id:" << current_event_id << " _event_id:" << _event_id << " event_id_group:" << event_id_group << endl;
              if (current_event_id == 0) {
                cout << "id error" << endl;
                getchar();
              }
              sprintf(canvas_name, "selected_event_id_%lu",
                      current_event_id - 1); // already move to next event's header
              strcpy(canvas_output_name, canvas_name);
              strcat(canvas_output_name, ".png");
              hit_display_cache = new TCanvas(canvas_name,canvas_name, 1200, 
                                              480);
              hit_display_cache->Divide(1, 1);
              hit_display_cache->cd(1);
              track_baseline->SetNameTitle(canvas_name, canvas_name);
              track_base->DrawClonePad();

              if (true) {
                event_track[temp_track_group]->WriteTObject(hit_display_cache);
                // cout << track_group_name << endl;
                #ifdef SAVE_TRACKS_OUT_OF_ROOT
                  chdir(track_group_name);
                    hit_display_cache->SaveAs(canvas_output_name);
                    chdir("..");
                #endif
                total_tracks++;
              }
              #ifdef DECODER_DEBUG
                cout << "Paused. Enter 'e' to exit or other key to continue." << endl;
                if (getchar() == 'e') {
                  cout << "Program is stopped by user." << endl;
                  return 1;
                }
              #endif
              delete hit_display_cache;
              for (int signal_id = 0; signal_id < signal_cache_length;
                   signal_id++) {
                if (signal_cache[signal_id][0] == 4) {
                  delete hit_model[signal_id];
                }
              }
            }
          }
        }

        #ifdef DECODER_DEBUG
          for (int signal_id = 0; signal_id < signal_cache_length; signal_id++) {
            for (int i = 0; i < 5; i++) {
              cout << signal_cache[signal_id][i] << " " ;
            }
          }

          cout << "Paused. Enter 'e' to exit or other key to continue." << endl;
          if (getchar() == 'e') {
          cout << "Program is stopped by user." << endl;
            return 1;
          }
        #endif


        for (int trigger_id = 0; trigger_id < trigger_cache_length;
             trigger_id++) {
          p_HPTDC_data->Fill(trigger_cache[trigger_id][0],
                             trigger_cache[trigger_id][1],
                             trigger_cache[trigger_id][2],
                             trigger_cache[trigger_id][3],
                             trigger_cache[trigger_id][4],
                             false);
          for (int i = 0; i < 5; i++) {
            trigger_cache[trigger_id][i] = 0;
          }
        }
        for (int signal_id = 0; signal_id < signal_cache_length; signal_id++) {
          p_HPTDC_data->Fill(signal_cache[signal_id][0],
                             signal_cache[signal_id][1],
                             signal_cache[signal_id][2],
                             signal_cache[signal_id][3],
                             signal_cache[signal_id][4],
                             signal_flag[signal_id]);
          for (int i = 0; i < 5; i++) {
            signal_cache[signal_id][i] = 0;
          }
        }
        
        trigger_cache_length = 0;
        signal_cache_length = 0;
      }
            
      #ifdef DECODER_DEBUG
        cout << endl; 
        cout << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
        cout << endl << "header= " << header_type << endl;
        cout << data_cache << endl;
        cout << _type << " || " << _tdc << " || " << _event_id << " || " << _bunch_id << endl << endl;
      #endif

      // save data to file if got data
        data_out_flow << _type << " " << _tdc << " " << 0 << " " << current_event_id << " " <<  _bunch_id << endl;

      // fill ntuple
      p_HPTDC_data->Fill(_type, _tdc, 0, current_event_id, _bunch_id, false);

    } else if (header_type == 4 || header_type == 5) {
      type = word >> 28;
      tdc = word >> 24;
      channel = word >> 19;
      coarse = word >> 7;
      fine = word;
      _type = static_cast<unsigned int>((type.to_ulong()));
      _tdc = static_cast<unsigned int>((tdc.to_ulong()));
      _channel = static_cast<unsigned int>((channel.to_ulong()));
      _coarse = static_cast<unsigned int>((coarse.to_ulong()));
      _fine = static_cast<unsigned int>((fine.to_ulong()));

      // collect data
      new_data_line[0] = _type;
      new_data_line[1] = _tdc;
      new_data_line[2] = _channel;
      new_data_line[3] = _coarse;
      new_data_line[4] = _fine;
      if (trigger_cache_length > 128 || signal_cache_length > 128) {
        cout << "warning: data overflow! Set to Max" << endl;
        trigger_cache_length = 128;
        signal_cache_length = 128;
      }

      /*CollectOneEvent(&trigger_cache_length, trigger_cache,
                        &signal_cache_length, signal_cache, new_data_line);*/

      if (new_data_line[1] == 1) {
        trigger_cache[trigger_cache_length][0] = new_data_line[0];
        trigger_cache[trigger_cache_length][1] = new_data_line[1];
        trigger_cache[trigger_cache_length][2] = new_data_line[2];
        trigger_cache[trigger_cache_length][3] = new_data_line[3];
        trigger_cache[trigger_cache_length][4] = new_data_line[4];
        trigger_cache_length++;
      }
      if (new_data_line[1] != 1) {
        signal_cache[signal_cache_length][0] = new_data_line[0];
        signal_cache[signal_cache_length][1] = new_data_line[1];
        signal_cache[signal_cache_length][2] = new_data_line[2];
        signal_cache[signal_cache_length][3] = new_data_line[3];
        signal_cache[signal_cache_length][4] = new_data_line[4];
        signal_cache_length++;
      }

      #ifdef DECODER_DEBUG
        cout << "header= " << header_type << endl;
        cout << data_cache << endl;
        cout << _type << " || " << _tdc << " || " << _channel << " || " << _coarse << " || " << _fine << endl;
        // cout << "trigger_cache_length: " << trigger_cache_length << endl;
        // cout << "signal_cache_length: " << signal_cache_length << endl;
        //getchar();
      #endif

      // save data to file if got data
      data_out_flow << _type << " " << _tdc << " " << _channel << " " <<  _coarse << " " << _fine << endl;

      // fill ntuple
      //p_HPTDC_data->Fill(_type, _tdc, _channel, _coarse, _fine);

      // fill histogram
      time_in_ns = (_coarse + _fine / 128.0 ) * 25.0;
      if (header_type == 4) {
        p_leading_time->Fill(time_in_ns);
      }
      if (header_type == 5) {
        p_trailing_time->Fill(time_in_ns);
      }
    } else {
      #ifdef DECODER_DEBUG
        cout << endl << "header= " << header_type << endl;
        cout << data_cache << endl;
      #endif
        //cout << endl << "++++++++++++++++++++++++++++++++++++++++++++" << endl << "++++++++++++++++++++++++++++++++++++++++++++" << endl << "unknown type" << endl << "++++++++++++++++++++++++++++++++++++++++++++" << endl << "++++++++++++++++++++++++++++++++++++++++++++" << endl;
    }
    total_loop++;
    //cout << "loop: " << total_loop << endl;
    #ifdef SET_MAXWORDS
      if (total_loop % 10000 == 0) {
        cout << "Decoding progress: "
             << total_loop * 4.0 * 100 / (MAX_WORD_QUANTITY * 4.0) << " %"
             << endl;
      }
    #else
      if (total_loop % 10000 == 0) {
        cout << "Decoding progress: "
             << total_loop * 4.0 * 100 / data_in_flow_length << " %"
             << endl;
      }
    #endif
    }

  cout << "Decoding completed !" << endl;

  // plot the time spectrum for leading and trailing edge
  cout << "Making plots... " << endl;
  cout << endl << "Plotting total leading edge spectrum... " << endl;
  TCanvas *p_leading_canvas = new TCanvas("leading", "leading", 0, 0, 800,
                                          450);
  p_leading_canvas->cd();
  p_leading_time->GetXaxis()->SetTitle("time/ns");
  p_leading_time->Draw();

  cout << endl << "Plotting total trailing edge spectrum... " << endl;
  TCanvas *p_trailing_canvas = new TCanvas("trailing", "trailing", 0, 450,
                                           800, 450);
  p_trailing_canvas->cd();
  p_trailing_time->GetXaxis()->SetTitle("time/ns");
  p_trailing_time->Draw();

  cout << endl << "Plotting hits distribution... " << endl;
  TCanvas *p_hits_canvas[MAX_TUBE_LAYER];
  char canvas_name[256];
  for (Int_t layer_id = 0; layer_id != MAX_TUBE_LAYER; layer_id++) {
    sprintf(canvas_name, "layer_%d_hits_distribution", layer_id);
    p_hits_canvas[layer_id] = new TCanvas(canvas_name, canvas_name);
    p_hits_canvas[layer_id]->cd();
    p_hits_distribution[layer_id]->Draw();
  }

  p_output_rootfile->Write();

  // export data to output directory
  #ifdef SAVE_TRACKS_OUT_OF_ROOT
    TCanvas *p_output_canvas = new TCanvas("", "");
    p_output_canvas->cd();
    for (Int_t tdc_id = 0; tdc_id != MAX_TDC_QUANTITY; tdc_id++) {
      sprintf(directory_name, "TDC_%02d_of_%02d_Time_Spectrum", tdc_id,
              MAX_TDC_QUANTITY);
      chdir(directory_name);

      p_tdc_tdc_time[tdc_id]->Draw();
      sprintf(output_filename, "tdc_%d_tdc_time_spectrum.png", tdc_id);
      p_output_canvas->SaveAs(output_filename);

      p_tdc_tdc_time_original[tdc_id]->Draw();
      sprintf(output_filename, "tdc_%d_tdc_time_spectrum_original.png", tdc_id);
      p_output_canvas->SaveAs(output_filename);

      p_tdc_tdc_time_corrected[tdc_id]->Draw();
      sprintf(output_filename, "tdc_%d_tdc_time_spectrum_corrected.png", tdc_id);
      p_output_canvas->SaveAs(output_filename);

      p_tdc_adc_time[tdc_id]->Draw();
      sprintf(output_filename, "tdc_%d__adc_time_spectrum.png", tdc_id);
      p_output_canvas->SaveAs(output_filename);

      p_tdc_channel[tdc_id]->Draw();
      sprintf(output_filename, "tdc_%d__channel_distribution.png", tdc_id);
      p_output_canvas->SaveAs(output_filename);
        
      for (Int_t channel_id = 0; channel_id != MAX_TDC_CHANNEL_QUANTITY;
           channel_id++) {
        p_tdc_time[tdc_id][channel_id]->Draw();
        sprintf(output_filename, "tdc_%d__channel_%d__tdc_time_spectrum.png", 
                tdc_id, channel_id);
        p_output_canvas->SaveAs(output_filename);

        p_tdc_time_original[tdc_id][channel_id]->Draw();
        sprintf(output_filename, 
                "tdc_%d__channel_%d__tdc_time_spectrum_original.png", 
                tdc_id, channel_id);
        p_output_canvas->SaveAs(output_filename);

        p_tdc_time_corrected[tdc_id][channel_id]->Draw();
        sprintf(output_filename, 
                "tdc_%d__channel_%d__tdc_time_spectrum_corrected.png", 
                tdc_id, channel_id);
        p_output_canvas->SaveAs(output_filename);
            
        p_adc_time[tdc_id][channel_id]->Draw();
        sprintf(output_filename, "tdc_%d__channel_%d__adc_time_spectrum.png",
                tdc_id, channel_id);
        p_output_canvas->SaveAs(output_filename);
      }
      chdir("..");
    }
  #endif
    
  delete p_output_rootfile;
  delete track_base;
  delete p_leading_canvas;
  delete p_trailing_canvas;

  strcpy(output_filename, filename);
  strcat(output_filename, ".out");
  cout << endl;
  cout << "Text data saved to: " << output_filename << endl;
  cout << "NTuple data saved to: " << output_root_filename << endl;
  cout << endl << "Work done." << endl;

  if (batch_mode_flag == false) {
  gROOT->SetBatch(kFALSE); // resume setup
  }
  return 0;
}

double calibrated_p(const double &w) {
  return std::exp(1.11925e+00 + 2.08708e-02 * w);
}

double ADC_Correction(double w, vector<Double_t> v_LUT) {
  if( w > 400.0 || w < 0 ) return 0;
  return v_LUT[static_cast<int>(w)];
}